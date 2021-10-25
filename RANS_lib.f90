MODULE RANS_lib 
    CONTAINS 

    SUBROUTINE int_linspace(first,last,dx,x)
        USE iso_fortran_env

        ! Numerical Precision
        INTEGER ,PARAMETER :: DBL = SELECTED_REAL_KIND(p=13)

        ! Counter 
        INTEGER :: nx 

        INTEGER, INTENT(IN) :: first,last,dx 
        INTEGER,ALLOCATABLE, INTENT(OUT) :: x(:)

        nx = (last-first)/dx + 1

        ALLOCATE(x(nx))
        x(1) = first 

        DO i = 2,nx 
            x(i) = x(i-1) + dx
        end DO 
    
    END SUBROUTINE int_linspace

    SUBROUTINE linspace(first,last,dx,x)
        USE iso_fortran_env

        ! Numerical Precision
        INTEGER ,PARAMETER :: DBL = SELECTED_REAL_KIND(p=13)

        ! Counter 
        INTEGER :: nx 

        REAL(KIND = REAL64), INTENT(IN) :: first,last,dx 
        REAL(KIND = REAL64),ALLOCATABLE, INTENT(OUT) :: x(:)

        nx = NINT(ABS((last-first)/dx + 1))

        ALLOCATE(x(nx))
        x(1) = first 

        DO i = 2,nx 
            x(i) = x(i-1) + dx
        end DO 

    END SUBROUTINE linspace 

    SUBROUTINE Decompose(nx,PID,np,no_node,tmp)
        IMPLICIT NONE 

        INTEGER, INTENT(IN) :: nx, np, PID 
        INTEGER, INTENT(OUT) :: no_node
        INTEGER, ALLOCATABLE :: tmp(:)

        INTEGER :: rem 

        no_node = nx/np

        IF (PID .EQ. np - 1) then 
            rem = nx - no_node*PID 
            no_node = rem
        end IF

        IF (PID .NE. np-1) then 
            CALL int_linspace(PID*no_node+1, PID*no_node+no_node,1,tmp)
        Else
            Call int_linspace(nx-no_node+1,nx ,1,tmp)
        end IF

    END SUBROUTINE 

    SUBROUTINE sendTo(PID,np,left,right) !######################## LEGACY #######################!
        ! Check current PID location and neighbouring PID
        USE MPI
        IMPLICIT NONE 

        INTEGER, INTENT(IN) :: PID, np 
        INTEGER, INTENT(OUT) :: left, right 

        IF (PID == 0) then
            left = MPI_PROC_NULL ! IF boundary PID, end isend/irecv immediately
            right = PID + 1
        ELSE IF (PID == np - 1) then 
            left = PID - 1
            right = MPI_PROC_NULL ! IF boundary PID, end isend/irecv immediately
        ELSE
            left = PID - 1
            right = PID + 1
        end IF 

    end SUBROUTINE

    SUBROUTINE mxlen(nhx,nhy,h,y,ml_out)
        ! Calculate Mixing length 
        IMPLICIT none 

        REAL(kind = 8), ALLOCATABLE, INTENT(IN) :: y(:)
        REAL(kind = 8), INTENT(IN) :: h
        INTEGER , INTENT(IN) :: nhx,nhy 
        REAL(kind = 8), ALLOCATABLE, INTENT(OUT) :: ml_out(:,:)

        REAL(kind = 8), ALLOCATABLE :: ml(:)
        REAL(kind = 8) :: del 

        del = h/2 

        ALLOCATE(ml(nhy))
        ALLOCATE(ml_out(nhy,nhx))

        ml = y 
        ml(nhy/2 + 1: nhy) = abs(y(nhy/2 + 1: nhy)-h)
        ml = del*(0.14-0.08*(1-ml/del)**2-0.06*(1-ml/del)**4);

        ml(nhy) = ml(1)

        ml_out = spread(ml,2,nhx)

    end SUBROUTINE

    SUBROUTINE init_mpi(np,neighbours_ranks,cart_comm)
        USE MPI 

        IMPLICIT NONE 

        INTEGER, INTENT(IN) :: np 
        INTEGER, INTENT(OUT) :: neighbours_ranks(0:3),cart_comm

        LOGICAL :: reorder, periods(0:1)
        INTEGER :: dims(0:1), ierr 
        INTEGER, PARAMETER :: DOWNc = 0
        INTEGER, PARAMETER :: UPc = 1
        INTEGER, PARAMETER :: LEFTc = 2
        INTEGER, PARAMETER :: RIGHTc = 3

        dims = (/0,0/)
        CALL MPI_Dims_create(np, 1, dims, ierr)
    
        ! Creates communicator on cartesian topology
        reorder = .TRUE.
        periods = (/.FALSE., .FALSE./)
        CALL MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, reorder, cart_comm, ierr)
    
        CALL MPI_Cart_shift(cart_comm, 1, 1, neighbours_ranks(LEFTc), neighbours_ranks(RIGHTc), ierr)

    end subroutine

end MODULE 