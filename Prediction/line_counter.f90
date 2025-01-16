        subroutine line_counter(fname,n)
        integer io,n
        character(len=50) :: fname
        open(45,file=fname,status='unknown',action='read')
        n=0

        do
          read(45,*,iostat=io)
          if(io/=0) exit
          n=n+1
        end do
        rewind(45)
        close(45)

        end
