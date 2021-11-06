PROGRAM RANS_MPI
    USE MPI 
    USE RANS_lib 
    USE readin 
    USE omp_lib

    IMPLICIT NONE 

    REAL(kind = 8) :: rho, nu, U_in  ! Flow properties 

    ! Domain
    INTEGER :: n, nhx, nhy 
    REAL(KIND = 8) :: x, y, hx, hy
    REAL(KIND = 8), ALLOCATABLE :: x_tec(:), y_tec(:)
    REAL(KIND = 8), ALLOCATABLE :: meshY(:)
    REAL(KIND = 8) :: dt 

    ! Velocity field 
    REAL(KIND = 8), ALLOCATABLE :: U(:,:),  V(:,:), Unew(:,:), Vnew(:,:), Uout(:,:), Vout(:,:), Pout(:,:)

    ! Pressure field
    REAL(KIND = 8), ALLOCATABLE :: P(:,:), Pc(:,:), trial(:,:), trial_lag(:,:)
    REAL(KIND = 8) :: prs_coeff

    ! Residual/Divergence field 
    REAL(KIND = 8) :: div_sum
    REAL(KIND = 8), ALLOCATABLE :: residual(:,:), div(:,:), div_gl(:,:)

    ! Turbulence model (eddy viscosity)
    REAL(KIND = 8), ALLOCATABLE :: delY_u(:,:), delYY_u(:,:), delYX_u(:,:),ml(:,:), tmp(:,:)

    ! Iteration/relaxation factor
    INTEGER :: n_count, p_count, i, j, k, ii, jj, kk
    REAL(KIND = 8) :: w, relx_pc

    ! Error/ Steady state criteria
    REAL(KIND = 8) :: U_change, U_change_max, resid_pc,resid_pc_max 

    ! Domain Decomposition
    INTEGER :: no_node 
    INTEGER :: ll, lh, lh_prs, ll_prs ! Local Low, Local High, local pressure low/high  
    INTEGER, ALLOCATABLE :: lglel(:) ! Local to GLoblal ELement

    ! Solver timer
    REAL :: t_prsslv1, t_prsslv2
    DOUBLE precision :: t_solve1, t_solve2, t2solve

    ! Data IO
    INTEGER :: UNIT 
    CHARACTER(LEN = 80) :: filename 

    ! MPI cart decompose
    INTEGER :: dims(0:1)

    ! MPI initialise 
    LOGICAL :: reorder, periods(0:1)
    INTEGER :: cart_comm
    INTEGER, PARAMETER :: DOWNc = 0
    INTEGER, PARAMETER :: UPc = 1
    INTEGER, PARAMETER :: LEFTc = 2
    INTEGER, PARAMETER :: RIGHTc = 3
    INTEGER :: neighbours_ranks(0:3)
    INTEGER :: ierr, pid, np, tag 
    INTEGER(KIND = 4) :: stat(MPI_STATUS_SIZE, 8), istat(MPI_STATUS_SIZE)
    INTEGER :: left, right, req(8)
    INTEGER, ALLOCATABLE :: displacement(:)

    tag = 1
    CALL MPI_Init(ierr)
    CALL MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)


    ! INITIALISE value 
    call valInit()
    call init_mpi(np,neighbours_ranks,cart_comm)

    IF ((np .GT. nhx) .OR. (np .GT. nhy)) then 
        IF(PID == 0) write(*,*) 'Too little cells for too many PIDs'
        call MPI_FINALIZE(ierr)
        return 
    end IF 

    left = neighbours_ranks(2)
    right = neighbours_ranks(3)

    CALL mxlen(nhx,nhy,y,meshY,ml)

    t_solve1 = MPI_Wtime()

    DO WHILE (U_change .GE. U_change_max)
        n_count = n_count + 1


        call eddy_visc(nhx,nhy,Unew,delY_u,ml,delYY_u,delYX_u,ll,lh,left,right,tag,cart_comm,req)
  
        ! U - Velocity solve
        DO i = 2,nhy - 1
            DO j = ll,lh 
                Unew(i,j)= (1-w)*U(i,j)+w*(U(i,j)-dt*(P(i+1,j)-P(i-1,j))/(2*hx) &
                +nu*dt*(1/(hx*hx)*(U(i+1,j)-2.*U(i,j)+U(i-1,j)) & 
                +1/(hy*hy)*(U(i,j+1)-2.*U(i,j)+U(i,j-1))) &
                -dt*U(i,j)/(hx)*(U(i,j)-U(i-1,j)) &
                -dt*V(i,j)/(2*hy)*(U(i,j+1)-U(i,j-1))+dt*delYY_u(i,j))
            end DO 
        end DO

        ! Update Boundary Condtion
        Unew(nhy,:)=Unew(nhy-1,:);
        Unew(1,2:nhx-1)=-Unew(2,2:nhx-1)+2*U_in;
        Unew(2:nhy,1)=-Unew(2:nhy,2);
        Unew(2:nhy,nhx)=-Unew(2:nhy,nhx-1);
        U_change= MAXVAL(ABS((Unew(:,ll:lh)-U(:,ll:lh))/dt)) ! MPI_Allreduce?

        CALL MPI_Allreduce(U_change, U_change, 1, MPI_DOUBLE_PRECISION, MPI_SUM, cart_comm, ierr)

        ! MPI send and recv
		CALL MPI_SEND_RECV(Unew,nhy,nhx,ll,lh,left,right,tag,cart_comm,req)

        CALL MPI_waitall(4,req,stat,ierr)

        ! V - Velocity solve
        DO i = 2,nhy - 1
            DO j = ll,lh
                Vnew(i,j)= (1-w)*V(i,j)+w*(V(i,j)-dt*(P(i,j+1)-P(i,j-1))/(2*hy) &
                +nu*dt*(1/(hx*hx)*(V(i+1,j)-2.*V(i,j)+V(i-1,j)) &
                +1/(hy*hy)*(V(i,j+1)-2.*V(i,j)+V(i,j-1))) &
                -dt*U(i,j)/(hx)*(V(i,j)-V(i-1,j)) &
                -dt*V(i,j)/(2*hy)*(V(i,j+1)-V(i,j-1))+dt*delYX_u(i,j))
            end DO
        end DO

        ! Update Boundary Condition 
        Vnew(nhy,:)=Vnew(nhy-1,:);
        Vnew(1,:)=-Vnew(2,:);
        Vnew(:,1)=-Vnew(:,2);
        Vnew(:,nhx)=-Vnew(:,nhx-1);

        ! MPI send and recv
		CALL MPI_SEND_RECV(Vnew,nhy,nhx,ll,lh,left,right,tag,cart_comm,req)

        CALL MPI_waitall(4,req,stat,ierr)

        ! Divergence solve
        DO i = 2,nhy - 1
            DO j = ll,lh 
                div(i,j)=(Unew(i+1,j)-Unew(i-1,j))/(2*hx) & 
                +(Vnew(i,j+1)-Vnew(i,j-1))/(2*hy);
            end DO 
        end DO

        ! MPI send and recv
		CALL MPI_SEND_RECV(div,nhy,nhx,ll,lh,left,right,tag,cart_comm,req)     

        CALL MPI_waitall(4,req,stat,ierr)


        div_sum = sum(abs(div))/((nhx-2)*(nhy-2))

        resid_pc = 1
        p_count = 0 
        Pc = 0.0
        trial = 0.0
        trial_lag = 0.0
        
        DO WHILE ((resid_pc .GT. resid_pc_max) .AND. (p_count .LT. 100)) 
            p_count = p_count + 1
            DO i = 2,nhy - 1
                DO j = ll,lh 
                    trial(i,j) = trial_lag(i,j)*(1-relx_pc)+relx_pc*0.25*((div(i,j)/dt - &
                        (trial_lag(i-1,j) + trial_lag(i+1,j))/hx**2 - &
                        (trial_lag(i,j-1)+trial_lag(i,j+1))/hy**2)*prs_coeff)
                end DO
            end DO

            resid_pc = sum(abs(trial(:,ll:lh)-abs(trial_lag(:,ll:lh)))/((nhx-2)*(nhy-2)))
            CALL MPI_Allreduce(resid_pc, resid_pc, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

            trial_lag = trial
            call MPI_SEND_RECV(trial,nhy,nhx,ll,lh,left,right,tag,cart_comm,req)
    
            CALL MPI_waitall(4,req,stat,ierr)

        end DO


        Pc = prs_coeff*trial+trial

        Pc(1,:) = Pc(2,:)
        Pc(:,1) = Pc(:,2)
        Pc(:,nhx) = Pc(:,nhx-1)
        Pc(nhy,:) = 0 

        ! Update Pressure
        P = P + Pc 

        DO i = 2,nhy - 1
            DO j = ll,lh
                Unew(i,j)= Unew(i,j)-dt*(Pc(i+1,j)-Pc(i-1,j))/(2*hx)
                Vnew(i,j)= Vnew(i,j)-dt*(Pc(i,j+1)-Pc(i,j-1))/(2*hy)
            end DO
        end DO

        ! Velocity Boundary Condition
        Unew(nhy,:)=Unew(nhy-1,:);
        Unew(1,2:nhx-1)=-Unew(2,2:nhx-1)+2*U_in;
        Unew(2:nhy,1)=-Unew(2:nhy,2);
        Unew(2:nhy,nhx)=-Unew(2:nhy,nhx-1);
        
        Vnew(nhy,:)=Vnew(nhy-1,:);
        Vnew(1,:)=-Vnew(2,:);
        Vnew(:,1)=-Vnew(:,2);
        Vnew(:,nhx)=-Vnew(:,nhx-1);
        
        U=Unew;
        V=Vnew;

        ! MPI send and recv
		CALL MPI_SEND_RECV(Unew,nhy,nhx,ll,lh,left,right,tag,cart_comm,req)
    
        CALL MPI_waitall(4,req,stat,ierr)
    
	    ! MPI send and recv
		CALL MPI_SEND_RECV(Vnew,nhy,nhx,ll,lh,left,right,tag,cart_comm,req)
    
        CALL MPI_waitall(4,req,stat,ierr)  

        if (PID == 0) then 
                WRITE(*,'(a40)') '-----------------------------------------'
                WRITE(*,'(a20,i20.1)') '# Current Iteration = ', n_count
                WRITE(*,'(a20,E20.6)') '# Current Delta U = ', U_change
                WRITE(*,'(a20,i20.6)') '# Pressure count = ', p_count
        end if 

    end DO

    t_solve2 = MPI_Wtime()
    t2solve = t_solve2 - t_solve1

    ! Collect data to write to disc
    IF (PID == 0) then 
        Uout(:,ll:lh) = U(:,ll:lh)
        DO kk = 1, np - 1
            call Decompose(nhx,kk,np,no_node,lglel)
            ll_prs = lglel(1) 
            lh_prs = lglel(size(lglel))
            if (kk == np - 1) lh_prs = lh_prs - 1
            CALL MPI_Recv(Uout(:,ll_prs:lh_prs),(lh_prs-ll_prs+1)*nhy, MPI_DOUBLE_PRECISION,kk,tag,cart_comm,istat,ierr)
        end DO
    ELSE 
        CALL MPI_Send(U(:,ll:lh),(lh-ll+1)*nhy,MPI_DOUBLE_PRECISION,0,tag,cart_comm,ierr)
    end IF

    call MPI_Bcast(Uout,nhx*nhy,MPI_DOUBLE_PRECISION,0,cart_comm,ierr)

    CALL MPI_FINALIZE(ierr)

    CALL write_to_screen(PID,nhx,nhy,n_count,t2solve)

    ! Data IO
    filename = 'Uout.dat'
    UNIT = 1001 

    OPEN(UNIT,FILE = filename, FORM = 'FORMATTED', STATUS = 'REPLACE')
        WRITE(UNIT,*) Uout 
    CLOSE(UNIT)

    filename = 'Vout.dat'
    UNIT = 1001 

    OPEN(UNIT,FILE = filename, FORM = 'FORMATTED', STATUS = 'REPLACE')
        WRITE(UNIT,*) V 
    CLOSE(UNIT)

    filename = 'Pout.dat'
    UNIT = 1001 

    OPEN(UNIT,FILE = filename, FORM = 'FORMATTED', STATUS = 'REPLACE')
        WRITE(UNIT,*) P 
    CLOSE(UNIT)

Contains
SUBROUTINE AllocateMemory()
	IMPLICIT NONE
	ALLOCATE(U(nhy,nhx))
	ALLOCATE(Unew(nhy,nhx))
	ALLOCATE(Uout(nhy,nhx))
	ALLOCATE(V(nhy,nhx))
	ALLOCATE(Vnew(nhy,nhx))
	ALLOCATE(Vout(nhy,nhx))
	ALLOCATE(P(nhy,nhx))
	ALLOCATE(Pc(nhy,nhx))
	ALLOCATE(Pout(nhy,nhx))
	ALLOCATE(residual(nhy,nhx))
	ALLOCATE(div(nhy,nhx))
	ALLOCATE(div_gl(nhy,nhx))
	ALLOCATE(meshY(nhy))
	ALLOCATE(delY_u(nhy,nhx))
	ALLOCATE(delYY_u(nhy,nhx))
	ALLOCATE(delYX_u(nhy,nhx))
	ALLOCATE(ml(nhy,nhx))
	ALLOCATE(tmp(nhy,nhx))

    ALLOCATE(trial(nhx,nhy))
    ALLOCATE(trial_lag(nhx,nhy))
END SUBROUTINE AllocateMemory

SUBROUTINE valInit()
	IMPLICIT NONE
	rho = 1
    nu = 0.001
    U_in = 60.0

    ! n = 4 
    ! x = 8 ! Domain length
    ! y = 0.1
    ! hx = x/2**n 
    ! hy = y/2**n 
    ! nhx = NINT(x/hx) + 2
    !nhy = NINT(y/hy) + 2
    dt = 0.000001 ! may work for bigger dt 

    call userinp(x,y,hx,hy,nhx,nhy,nu,rho,U_change_max,resid_pc_max)

    call Decompose(nhx,PID,np,no_node,lglel)

    ll = lglel(1) 
    lh = lglel(size(lglel))

    if (PID == 0) then 
        ll = 2
    end if

    if (PID == np - 1) lh = lh - 1
	CALL AllocateMemory()
	CALL linspace(-hy/2,0.1+hy/2,hy,meshY)

    U = 0.0
    U(2:nhy - 1,1:2) = U_in 
    Unew = U
    V = 0.0
    Vnew = V 
    residual = 0.0 
    div = 0.0
    delY_u = 0.0
    delYY_u = 0.0
    delYX_u = 0.0
    P = 0.0
    P(:,1) = 1.0
    Pc = P 

    U_change = 1
    U_change_max = 1e-4
    resid_pc = 1
    resid_pc_max = 1e-4

    w = 1.9
    relx_pc = 1.5
    prs_coeff = -2/hx**2 - 2/hy**2
    prs_coeff = 1/prs_coeff

    trial(:,:) = 0
    trial_lag(:,:) = 0 
END SUBROUTINE valInit

SUBROUTINE eddy_visc(nhx,nhy,Unew,delY_u,ml,delYY_u,delYX_u,ll,lh,left,right,tag,cart_comm,req)
	
	IMPLICIT NONE
	INTEGER :: nhx, nhy 
	REAL(KIND = 8),DIMENSION(nhy,nhx), INTENT(IN) :: Unew,ml
	REAL(KIND = 8),DIMENSION(nhy,nhx), INTENT(INOUT) :: delYY_u,delYX_u,delY_u
	INTEGER :: i,j
	INTEGER :: ll, lh
	INTEGER :: cart_comm,tag
	INTEGER :: left, right,ierr, req(8)
	
	
         DO i = 2,nhy - 1
             DO j = ll,lh
             delY_u(i,j) = (Unew(i+1,j)-Unew(i-1,j))/(2*hy)
             end DO 
         end DO
         delY_u(nhy,:) = 0
         delY_u(1,2:nhx-1) = 0
         delY_u(2:nhy,1) = -delY_u(2:nhy,2)
         delY_u(2:nhy,nhx) = -delY_u(2:nhy,nhx-1)
		! MPI send and recv
		CALL MPI_SEND_RECV(delY_u,nhy,nhx,ll,lh,left,right,tag,cart_comm,req)

        CALL MPI_waitall(4,req,stat,ierr)

        tmp = ml**2*ABS(delY_u(:,ll:lh))*delY_u(:,ll:lh)
        DO i = 2,nhy - 1
            DO j = ll,lh 
                tmp(i,j) = ml(i,j)**2*ABS(delY_u(i,j))*delY_u(i,j)
            end DO
        end DO
        ! MPI send and recv
		CALL MPI_SEND_RECV(tmp,nhy,nhx,ll,lh,left,right,tag,cart_comm,req)

        CALL MPI_waitall(4,req,stat,ierr)

        DO i = 2,nhy - 1
            DO j = ll,lh
                delYY_u(i,j) = (tmp(i+1,j)-tmp(i-1,j))/(2*hy);
            end DO
        end DO

		! MPI send and recv
		CALL MPI_SEND_RECV(delYY_u,nhy,nhx,ll,lh,left,right,tag,cart_comm,req)

        CALL MPI_waitall(4,req,stat,ierr)

        DO i = 2,nhy - 1
            DO j = ll,lh
                delYX_u(i,j) = (tmp(i,j+1)-tmp(i,j-1))/(2*hx);
            end DO
        end DO

		! MPI send and recv
		CALL MPI_SEND_RECV(delYX_u,nhy,nhx,ll,lh,left,right,tag,cart_comm,req)

        CALL MPI_waitall(4,req,stat,ierr)


END SUBROUTINE eddy_visc

END PROGRAM RANS_MPI 