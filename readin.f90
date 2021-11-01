module readin 
    contains

    subroutine userinp(x,y,nhx,nhy,nu,density,U_change_max,resid_pc_max)
        implicit none

        ! IO variables
        character(128) :: buff = '' 
        character(128) :: buff1 
        integer :: stat , stat_inp
        integer :: line_no = 0
        logical :: end_of_file = .FALSE. 
    
        ! Mesh variables
        REAL(KIND = 8), INTENT(OUT) :: x, y 
        INTEGER, INTENT(OUT) :: nhx, nhy 
    
        ! Physical variables
        REAL(KIND = 8), INTENT(OUT) :: nu, density
    
        ! Tolerance
        REAL(KIND = 8), INTENT(OUT) :: U_change_max, resid_pc_max
    
        ! Turbulent model
        character(128) :: turb_model
    
        integer :: fff,i  
    
        open(1001,file = 'user.inp')
        fff = 1 
    
        ! Read nx and ny 
        do i = 1,21
            read(1001,'(a)',iostat = stat_inp) buff 
    
            line_no = line_no + 1   
    
            ! Convert tabs to white spaces  
            call StripSpaces(buff)
    
            
           ! Read mesh
            if (line_no == 2) then
                buff1 = buff(6:10)
                write(*,*) buff1
                read(buff1,*,iostat=stat) x 
            end if
            if (line_no == 3) then
                buff1 = buff(6:10)
                read(buff1,*,iostat=stat) y
            end if        
            if (line_no == 4) then
                buff1 = buff(6:10)
                read(buff1,*,iostat=stat) nhx 
            end if
            if(line_no == 5) then 
                buff1 = buff(6:10)
                read(buff1,*,iostat=stat) nhy
            end if 
    
            ! Read physical properties
            if(line_no == 8) then 
                buff1 = buff(5:10)
                write(*,*) buff1 
                read(buff1,*,iostat=stat) nu 
            endif 
            if(line_no == 9) then
                buff1 = buff(10:20)
                read(buff1,*,iostat=stat) density
            endif 
    
            ! Tolerances
            if(line_no == 12) then
                buff1 = buff(13:20)
                read(buff1,*,iostat=stat) U_change_max
            endif 
            if(line_no == 13) then
                buff1 = buff(20:30)
                read(buff1,*,iostat=stat) resid_pc_max
            endif 
    
            ! Read turbulent model 
            if (line_no == 19) then 
                buff1 = buff(1:20)
                read(buff1,*,iostat=stat) turb_model
            endif             
    
            !write(*,*) buff, line_no
    
            ! if (stat_inp.lt.0) then
            !     end_of_file = .true.
            !     goto 9999
            ! endif 
        end do 
    
        write(*,*) 'x length: ' , x 
        write(*,*) 'y length: ', y
        write(*,*) 'number of x grid: ', nhx 
        write(*,*) 'number of y grid: ', nhy
        write(*,*) 'viscousity :', nu 
        write(*,*) 'Density: ', density
        write(*,*) 'U tolerance: ', U_change_max
        write(*,*) 'Pressure tolerance: ', resid_pc_max

    end subroutine

    subroutine StripSpaces(string)
        character(len=*) :: string
        integer :: stringLen 
        integer :: last, actual

        stringLen = len (string)
        last = 1
        actual = 1

        do while (actual < stringLen)
            if (string(last:last) == ' ') then
                actual = actual + 1
                string(last:last) = string(actual:actual)
                string(actual:actual) = ' '
            else
                last = last + 1
                if (actual < last) &
                    actual = last
            endif
        end do

        end subroutine
end module