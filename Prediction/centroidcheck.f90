        subroutine centroidcheck(a,filename,nAtoms,smc,flag)
        implicit none
        real,dimension(nAtoms) :: x,y,z
        character(len=2),dimension(nAtoms) :: atms
        integer, dimension(smc) :: a
        real, dimension(smc) :: dist
        integer flag,i,j,nAtoms,smc,p,q
        real rij,xsum,ysum,zsum,maxRc,minRc,Rc
        character(len=50) :: filename
        integer f
       
        call file_reading(atms,x,y,z,nAtoms,filename)

        xsum=(sum(x(a(1:smc))))/dble(smc)
        ysum=(sum(y(a(1:smc))))/dble(smc)
        zsum=(sum(z(a(1:smc))))/dble(smc)

!       f=0
!       write(454,*) nAtoms+3
!       write(454,*) " " 
!       do j=1,nAtoms
!          f=0
!           do i=1,smc
!           if(j==a(i)) f=1
!           end do
!           if(f==0)  write(454,*) atms(j),x(j),y(j),z(j)
!           if(f==1)  write(454,*) "Li",x(j),y(j),z(j)
!       end do
!       write(454,*) "X",xsum,ysum,zsum

        do i=1,smc
        dist(i)=sqrt((x(a(i))-xsum)**2+(y(a(i))-ysum)**2+ &
     & (z(a(i))-zsum)**2)
        end do

!        Rc=minval(dist)
         Rc=1.2
!        Rc=2.0
!        write(*,*) "Cutoff from centroid",Rc

        flag=0
        do i=1,nAtoms
           
           rij=sqrt((x(i)-xsum)**2+(y(i)-ysum)**2+(z(i)-zsum)**2)
           if(rij<=Rc) then
           flag=1
           end if
        end do

!        if(flag==1) write(*,*) "Failed centroid",a(1:smc)
        end

