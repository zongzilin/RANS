PROGRAM RANS_MPI
    USE MPI 
    USE RANS_lib 

    IMPLICIT NONE 

    REAL(kind = 8) :: rho, nu, U_in  ! Flow properties 

    ! Domain
    INTEGER :: n, nhx, nhy 
    REAL(KIND = 8) :: x, y, hx, hy
    REAL(KIND = 8), ALLOCATABLE :: meshY(:)
    REAL(KIND = 8) :: dt 

    ! Velocity field 
    REAL(KIND = 8), ALLOCATABLE :: U(:,:),  V(:,:), Unew(:,:), Vnew(:,:)

    ! Pressure field
    REAL(KIND = 8), ALLOCATABLE :: P(:,:), Pc(:,:)

    ! Residual/Divergence field 
    REAL(KIND = 8) :: div_sum
    REAL(KIND = 8), ALLOCATABLE :: residual(:,:), div(:,:)

    ! Turbulence model (eddy viscosity)
    REAL(KIND = 8), ALLOCATABLE :: delY_u(:,:), delYY_u(:,:), delYX_u(:,:),ml(:,:), tmp(:,:)

    ! Iteration/relaxation factor
    INTEGER :: n_count, p_count, i, j, k, ii, jj 
    REAL(KIND = 8) :: w, relx_pc

    ! Error/ Steady state criteria
    REAL(KIND = 8) :: U_change, U_change_max, resid_pc,resid_pc_max 

    ! Domain Decomposition
    INTEGER :: no_node 
    INTEGER :: ll, lh ! Local Low, Local High 
    INTEGER, ALLOCATABLE :: lglel(:) ! Local to GLoblal ELement

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

    tag = 1
    CALL MPI_Init(ierr)
    CALL MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)


    ! INITIALISE value 
    ! ###########################  USE SUBROUTINE !!!!!! ###########################!
    rho = 1
    nu = 0.001
    U_in = 60.0

    n = 4 
    x = 8 ! Domain length
    y = 0.1
    hx = x/2**n 
    hy = y/2**n 
    nhx = NINT(x/hx) + 2
    nhy = NINT(y/hy) + 2
    dt = 0.00001 ! may work for bigger dt 

    call Decompose(nhx,PID,np,no_node,lglel)

    ll = lglel(1) 
    lh = lglel(size(lglel))

    if (PID == 0) then 
        ll = 2
    end if

    if (PID == np - 1) lh = lh - 1

    ALLOCATE(U(nhy,nhx))
    ALLOCATE(Unew(nhy,nhx))
    ALLOCATE(V(nhy,nhx))
    ALLOCATE(Vnew(nhy,nhx))
    ALLOCATE(P(nhy,nhx))
    ALLOCATE(Pc(nhy,nhx))
    ALLOCATE(residual(nhy,nhx))
    ALLOCATE(div(nhy,nhx))
    ALLOCATE(meshY(nhy))
    ALLOCATE(delY_u(nhy,nhx))
    ALLOCATE(delYY_u(nhy,nhx))
    ALLOCATE(delYX_u(nhy,nhx))
    ALLOCATE(ml(nhy,nhx))
    ALLOCATE(tmp(nhy,nhx))

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

    !################### TURN THIS INTO SUBROUNTINE call it mpicart #################################!
    dims = (/0,0/)
    CALL MPI_Dims_create(np, 1, dims, ierr)

    ! Creates communicator on cartesian topology
    reorder = .TRUE.
    periods = (/.FALSE., .FALSE./)
    CALL MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, reorder, cart_comm, ierr)

    CALL MPI_Cart_shift(cart_comm, 1, 1, neighbours_ranks(LEFTc), neighbours_ranks(RIGHTc), ierr)

    !write(*,*) 'this is PID', PID, 'with right neighour', neighbours_ranks(3), ' and left neighour ',neighbours_ranks(2)

    !##############################################################################################!

    CALL mxlen(nhx,nhy,y,meshY,ml) ! Mixing length uniform for all PID 

    !DO WHILE (U_change .GE. U_change_max)
    DO ii = 1,1
         n_count = n_count + 1

    !     ! ############## MADE THIS SECTION A SUBROUNTINE (CALL IT eddy_visc) ############################!
         DO i = 2,nhy - 1
             DO j = ll,lh
             delY_u(i,j) = (Unew(i+1,j)-Unew(i-1,j))/(2*hy)
             end DO 
         end DO
         delY_u(nhy,:) = 0
         delY_u(1,2:nhx-1) = 0
         delY_u(2:nhy,1) = -delY_u(2:nhy,2)
         delY_u(2:nhy,nhx) = -delY_u(2:nhy,nhx-1)

        call sendTo(PID,np,left,right) ! SHOULD CHANGE TO MPI CART SHIFT

        left = neighbours_ranks(2)
        right = neighbours_ranks(3)

        CALL MPI_iRecv(delY_u(1:nhy,ll-1),nhy,MPI_DOUBLE_PRECISION,left,tag,MPI_COMM_WORLD,req(2),ierr)
        CALL MPI_iSend(delY_u(1:nhy,ll),nhy,MPI_DOUBLE_PRECISION,left,tag,MPI_COMM_WORLD,req(1),ierr)
        CALL MPI_iRecv(delY_u(1:nhy,lh+1),nhy,MPI_DOUBLE_PRECISION,right,tag,MPI_COMM_WORLD,req(4),ierr)
        CALL MPI_iSend(delY_u(1:nhy,lh),nhy,MPI_DOUBLE_PRECISION,right,tag,MPI_COMM_WORLD,req(3),ierr)

        CALL MPI_waitall(4,req,stat,ierr)

        tmp = ml**2*ABS(delY_u)*delY_u

        DO i = 2,nhy - 1
            DO j = ll,lh
                delYY_u(i,j) = (tmp(i+1,j)-tmp(i-1,j))/(2*hy);
            end DO
        end DO

        CALL MPI_iRecv(delYY_u(1:nhy,ll-1),nhy,MPI_DOUBLE_PRECISION,left,tag,MPI_COMM_WORLD,req(2),ierr)
        CALL MPI_iSend(delYY_u(1:nhy,ll),nhy,MPI_DOUBLE_PRECISION,left,tag,MPI_COMM_WORLD,req(1),ierr)
        CALL MPI_iRecv(delYY_u(1:nhy,lh+1),nhy,MPI_DOUBLE_PRECISION,right,tag,MPI_COMM_WORLD,req(4),ierr)
        CALL MPI_iSend(delYY_u(1:nhy,lh),nhy,MPI_DOUBLE_PRECISION,right,tag,MPI_COMM_WORLD,req(3),ierr)

        CALL MPI_waitall(4,req,stat,ierr)

  

        DO i = 2,nhy - 1
            DO j = ll,lh
                delYX_u(i,j) = (tmp(i,j+1)-tmp(i,j-1))/(2*hx);
            end DO
        end DO

        CALL MPI_iRecv(delYX_u(1:nhy,ll-1),nhy,MPI_DOUBLE_PRECISION,left,tag,MPI_COMM_WORLD,req(2),ierr)
        CALL MPI_iSend(delYX_u(1:nhy,ll),nhy,MPI_DOUBLE_PRECISION,left,tag,MPI_COMM_WORLD,req(1),ierr)
        CALL MPI_iRecv(delYX_u(1:nhy,lh+1),nhy,MPI_DOUBLE_PRECISION,right,tag,MPI_COMM_WORLD,req(4),ierr)
        CALL MPI_iSend(delYX_u(1:nhy,lh),nhy,MPI_DOUBLE_PRECISION,right,tag,MPI_COMM_WORLD,req(3),ierr)

        CALL MPI_waitall(4,req,stat,ierr)

    !     !############################################################################################!

    !     ! U - Velocity solve
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
        U_change= MAXVAL(ABS((Unew-U)/dt)) ! MPI_Allreduce?

        CALL MPI_iRecv(Unew(1:nhy,ll-1),nhy,MPI_DOUBLE_PRECISION,left,tag,MPI_COMM_WORLD,req(2),ierr)
        CALL MPI_iSend(Unew(1:nhy,ll),nhy,MPI_DOUBLE_PRECISION,left,tag,MPI_COMM_WORLD,req(1),ierr)
        CALL MPI_iRecv(Unew(1:nhy,lh+1),nhy,MPI_DOUBLE_PRECISION,right,tag,MPI_COMM_WORLD,req(4),ierr)
        CALL MPI_iSend(Unew(1:nhy,lh),nhy,MPI_DOUBLE_PRECISION,right,tag,MPI_COMM_WORLD,req(3),ierr)

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

        CALL MPI_iRecv(Vnew(1:nhy,ll-1),nhy,MPI_DOUBLE_PRECISION,left,tag,MPI_COMM_WORLD,req(2),ierr)
        CALL MPI_iSend(Vnew(1:nhy,ll),nhy,MPI_DOUBLE_PRECISION,left,tag,MPI_COMM_WORLD,req(1),ierr)
        CALL MPI_iRecv(Vnew(1:nhy,lh+1),nhy,MPI_DOUBLE_PRECISION,right,tag,MPI_COMM_WORLD,req(4),ierr)
        CALL MPI_iSend(Vnew(1:nhy,lh),nhy,MPI_DOUBLE_PRECISION,right,tag,MPI_COMM_WORLD,req(3),ierr)

        CALL MPI_waitall(4,req,stat,ierr)

        ! Divergence solve
        DO i = 2,nhy - 1
            DO j = ll,lh 
                div(i,j)=(Unew(i+1,j)-Unew(i-1,j))/(2*hx) & 
                +(Vnew(i,j+1)-Vnew(i,j-1))/(2*hy);
            end DO 
        end DO

        CALL MPI_iRecv(div(1:nhy,ll-1),nhy,MPI_DOUBLE_PRECISION,left,tag,MPI_COMM_WORLD,req(2),ierr)
        CALL MPI_iSend(div(1:nhy,ll),nhy,MPI_DOUBLE_PRECISION,left,tag,MPI_COMM_WORLD,req(1),ierr)
        CALL MPI_iRecv(div(1:nhy,lh+1),nhy,MPI_DOUBLE_PRECISION,right,tag,MPI_COMM_WORLD,req(4),ierr)
        CALL MPI_iSend(div(1:nhy,lh),nhy,MPI_DOUBLE_PRECISION,right,tag,MPI_COMM_WORLD,req(3),ierr)

        CALL MPI_waitall(4,req,stat,ierr)

        if(PID == 0) write(*,*) div(:,2)

        div_sum = sum(abs(div))/((nhx-2)*(nhy-2))


        resid_pc = 1
        p_count = 0 
        Pc = 0.0

        !DO WHILE ((resid_pc .GT. resid_pc_max) .AND. (p_count .LT. 100)) 
        DO jj = 1,1
            p_count = p_count + 1
            DO i = 2,nhy - 1
                DO j = ll,lh
                    residual(i,j) = -1/(hx*hx)*(Pc(i+1,j)-2.*Pc(i,j)+Pc(i-1,j)) &
                                    -1/(hy*hy)*(Pc(i,j+1)-2.*Pc(i,j)+Pc(i,j-1)) &
                                    +(1./dt)*div(i,j)
                end DO
            end DO

            CALL MPI_iRecv(div(1:nhy,ll-1),nhy,MPI_DOUBLE_PRECISION,left,tag,MPI_COMM_WORLD,req(2),ierr)
            CALL MPI_iSend(div(1:nhy,ll),nhy,MPI_DOUBLE_PRECISION,left,tag,MPI_COMM_WORLD,req(1),ierr)
            CALL MPI_iRecv(div(1:nhy,lh+1),nhy,MPI_DOUBLE_PRECISION,right,tag,MPI_COMM_WORLD,req(4),ierr)
            CALL MPI_iSend(div(1:nhy,lh),nhy,MPI_DOUBLE_PRECISION,right,tag,MPI_COMM_WORLD,req(3),ierr)
    
            CALL MPI_waitall(4,req,stat,ierr)

            
        
            Pc(1,:)=Pc(2,:);
            Pc(:,1)=Pc(:,2);
            Pc(:,nhx)=Pc(:,nhx-1);
            Pc(nhy,:)=0;
            
            resid_pc = sum(abs(residual))/((nhx-2)*(nhy-2))

            !if (PID == 1) write(*,*) residual(:,ll)

        end DO


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

        ! if (PID == 0) then 
        !         WRITE(*,'(a40)') '-----------------------------------------'
        !         WRITE(*,'(a20,i20.1)') '# Current Iteration = ', n_count 
        !         WRITE(*,'(a20,E20.6)') '# Current Delta U = ', U_change
        ! end if 
        

    end DO

    CALL MPI_FINALIZE(ierr)

    ! ! Data IO
    ! filename = 'Uout.dat'
    ! UNIT = 1001 

    ! OPEN(UNIT,FILE = filename, FORM = 'FORMATTED', STATUS = 'REPLACE')
    !     WRITE(UNIT,*) U 
    ! CLOSE(UNIT)

    ! filename = 'Vout.dat'
    ! UNIT = 1001 

    ! OPEN(UNIT,FILE = filename, FORM = 'FORMATTED', STATUS = 'REPLACE')
    !     WRITE(UNIT,*) V 
    ! CLOSE(UNIT)

    ! filename = 'Pout.dat'
    ! UNIT = 1001 

    ! OPEN(UNIT,FILE = filename, FORM = 'FORMATTED', STATUS = 'REPLACE')
    !     WRITE(UNIT,*) P 
    ! CLOSE(UNIT)



END PROGRAM RANS_MPI 
