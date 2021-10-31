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

    
    SUBROUTINE MPI_SEND_RECV(Vector,nhy,nhx,ll,lh,left,right,tag,cart_comm,req)
		USE MPI
		IMPLICIT NONE
		
		REAL(KIND = 8), dimension(nhy,nhx) :: Vector
		INTEGER :: nhx, nhy 
		INTEGER :: ll, lh
		INTEGER :: cart_comm,tag
		INTEGER :: left, right,ierr, req(8)

		CALL MPI_iRecv(Vector(1:nhy,ll-1),nhy,MPI_DOUBLE_PRECISION,left,tag,cart_comm,req(2),ierr)
        CALL MPI_iSend(Vector(1:nhy,ll),nhy,MPI_DOUBLE_PRECISION,left,tag,cart_comm,req(1),ierr)
        CALL MPI_iRecv(Vector(1:nhy,lh+1),nhy,MPI_DOUBLE_PRECISION,right,tag,cart_comm,req(4),ierr)
        CALL MPI_iSend(Vector(1:nhy,lh),nhy,MPI_DOUBLE_PRECISION,right,tag,cart_comm,req(3),ierr)

	END SUBROUTINE 


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

    SUBROUTINE write_to_screen(PID,nhx,nhy,n_count,t2solve)

        IMPLICIT NONE 

        INTEGER, INTENT(IN) :: PID, nhy, nhx, n_count
        REAL(KIND = 8), INTENT(IN) :: t2solve

        if (PID == 0) then 
            WRITE(*,'(a40)') '-------RUN SUCCESSFUL, DYING......-----------'
            WRITE(*,'(a20,i20.1)') '# Number of Grids = ', nhx*nhy
            WRITE(*,'(a20,i20.1)') '# Iteration to converge = ', n_count
            !WRITE(*,'(a20,E20.6)') ' # Time in pressure solver = ', t_prsslv2 - t_prsslv1
            WRITE(*,'(a20,f20.6)') '# Time to converge = ', t2solve
        end if 
    end subroutine

    subroutine tecplot_2D ( iunit, nx, ny, x, y, T )

        IMPLICIT NONE
      
        Integer ( kind = 4 ) iunit,nx,ny,i,j,ierr
        Real ( kind = 8 ) T(nx,ny)
        Real ( kind = 8 ) x(nx)
        Real ( kind = 8 ) y(ny)
      
        Character(80), parameter ::  file_name = 'TecPlot2D.tec'
      
        open ( unit = iunit, file = file_name, form = 'formatted', access = 'sequential', status = 'replace', iostat = ierr )
      
        if ( ierr /= 0 ) then
          write ( *, '(a)' ) '  Error opening file : tecplot_2D '
          stop
        end if
         
        write ( iunit, '(a)' ) 'Title="' // trim ( 'Temperature Data' ) // '"'
        write ( iunit, '(a)' ) 'Variables=' // trim ( '"X","Y","T"' )
      
        write ( iunit, '(a)' ) ' '
        write ( iunit, '(a,i6,a,i6,a)' ) 'Zone I=', ny, ', J=', nx, ', F=POINT'
       
        do i = 1, nx
          do j = 1, ny
            write ( iunit, '(2f10.3,g15.6)' ) x(i), y(j), T(i,j)
          end do
        end do
        
        close ( unit = iunit )
      
      end subroutine tecplot_2D

end MODULE 

