        subroutine vecdecision(ix,filename,nAtoms,ix1)
        implicit none
        character(len=50) :: filename
        real,dimension(nAtoms) ::x,y,z
        real,dimension(3) ::temp
        character(len=2),dimension(nAtoms) :: atms
        integer nAtoms,i,j,p
        integer,dimension(3) :: ix
        integer ix1,loc

        real rij,rij1,rij2,maxv


        call file_reading(atms,x,y,z,nAtoms,filename)
        ! 1--> 1,2 ; 2--> 1,3 ; 3--> 2,3
        p=0
        maxv=0.0
        do i=1,3
           do j=i+1,3
           p=p+1
           rij=sqrt((x(ix(i))-x(ix(j)))**2 + &
                     (y(ix(i))-y(ix(j)))**2 + &
                     (z(ix(i))-z(ix(j)))**2)
           temp(p)=rij
           if(temp(p)>maxv) then
              maxv=temp(p)
              loc=p
           end if

           end do
       end do

       if(loc==1) then
          rij1=sqrt((x(ix(1))-x(ix(3)))**2 + &
      &             (y(ix(1))-y(ix(3)))**2 + &
      &             (z(ix(1))-z(ix(3)))**2)

          rij2=sqrt((x(ix(2))-x(ix(3)))**2 + &
       &            (y(ix(2))-y(ix(3)))**2 + &
       &            (z(ix(2))-z(ix(3)))**2)
           if(rij1>rij2) then
              ix1=1
           else
              ix1=2
           end if

        else if(loc==2) then

         rij1=sqrt((x(ix(1))-x(ix(2)))**2 + &
        &           (y(ix(1))-y(ix(2)))**2 + &
        &           (z(ix(1))-z(ix(2)))**2)

          rij2=sqrt((x(ix(2))-x(ix(3)))**2 + &
        &           (y(ix(2))-y(ix(3)))**2 + &
        &           (z(ix(2))-z(ix(3)))**2)
           if(rij1>rij2) then
              ix1=1
           else
              ix1=3
           end if

         else  if(loc==3) then

          rij1=sqrt((x(ix(1))-x(ix(2)))**2 + &
        &           (y(ix(1))-y(ix(2)))**2 + &
        &           (z(ix(1))-z(ix(2)))**2)

          rij2=sqrt((x(ix(1))-x(ix(3)))**2 + &
       &            (y(ix(1))-y(ix(3)))**2 + &
       &            (z(ix(1))-z(ix(3)))**2)
           if(rij1>rij2) then
              ix1=2
           else
              ix1=3
           end if
      end if

        ix1=ix(ix1)
        
!        write(*,*) ix1,"between",ix(1),ix(2),ix(3)
        end
