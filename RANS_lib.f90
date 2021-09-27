MODULE RANS_lib 
    CONTAINS 

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

    SUBROUTINE mxlen(nhy,h,y,ml)
        ! Calculate Mixing length (eddy viscousity)
        IMPLICIT none 

        REAL(kind = 8), ALLOCATABLE, INTENT(IN) :: y(:)
        REAL(kind = 8), INTENT(IN) :: h
        INTEGER , INTENT(IN) :: nhy 
        REAL(kind = 8), ALLOCATABLE, INTENT(OUT) :: ml(:)

        REAL(kind = 8) :: del 

        del = h/2 

        ALLOCATE(ml(nhy))

        ml = y 
        ml(nhy/2 + 1: nhy) = abs(y(nhy/2 + 1: nhy)-h)
        ml = del*(0.14-0.08*(1-ml/del)**2-0.06*(1-ml/del)**4);



    end SUBROUTINE

end MODULE 