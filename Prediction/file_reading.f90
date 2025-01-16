        subroutine file_reading(atms,x,y,z,n,fname)
        real, dimension(n) :: x,y,z
        character(len=2),dimension(n) :: atms
        character(len=50) :: fname
        
        open(55,file=fname)
        do i=1,n
        read(55,*) atms(i),x(i),y(i),z(i)
        end do
        close(55)

        end

