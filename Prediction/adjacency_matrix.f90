        subroutine adjacency_matrix(x,y,z,atms,nAtoms,con_mat)
        real,dimension(nAtoms) :: x,y,z
        character(len=2), dimension(nAtoms) :: atms
        integer,dimension(nAtoms,nAtoms) :: con_mat
        integer nAtoms
        integer i,j
        real rij

!       con_mat(1:nAtoms,1:nAtoms) = 0

        do i=1,nAtoms
           do j=1,nAtoms
           con_mat(i,j)=0
           end do
        end do

        do i=1,nAtoms
           do j=i+1,nAtoms
           rij=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
           if(rij<=3.8) then
           con_mat(i,j)=1
           con_mat(j,i)=1
           end if
           end do
        end do


        end


        
