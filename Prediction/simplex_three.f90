        subroutine simplex_three(nAtoms,con_mat,ncount,filename)
        implicit none
        integer,dimension(nAtoms,nAtoms) :: con_mat
        integer,dimension(3) :: a
        integer,dimension(3,3) :: new_mat

        integer nAtoms,ncount
                
        integer i,j,k
        integer p,q,smc
        integer flag,scol,flagm,flagp

        character(len=50) :: filename

        ncount=0
        smc=3

        do i=1,nAtoms
           do j=i+1,nAtoms
              do k=j+1,nAtoms
                a(1)=i
                a(2)=j
                a(3)=k
              
                do p=1,3
                   do q=1,3
                   new_mat(p,q)=con_mat(a(p),a(q))
                   end do
                end do

                flag=0
                do p=1,3
                scol=sum(new_mat(p,1:3))
                if(scol.ne.2) flag=1
                end do

                if(flag==0) then
                call vectorcheck(a,filename,nAtoms,flagm,smc)
                if(flagm==0) then
                   ncount=ncount+1
                   write(57,*) i,j,k,0,0,0
                else
                write(58,*) i-1,j-1,k-1
                end if

                end if

              end do  
           end do
         end do

        end
