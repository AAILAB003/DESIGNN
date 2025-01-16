        subroutine simplex_six(nAtoms,con_mat,ncount,filename)

        implicit none
        integer,dimension(nAtoms,nAtoms) :: con_mat
        integer,dimension(6) :: a
        integer,dimension(3) :: scl,ix,tempa
        integer,dimension(6,6) :: new_mat
        integer,dimension(6,6) :: upd_mat
        character(len=50) :: filename

        integer nAtoms,ncount
                
        integer i,j,k,l,m,n
        integer p,q,t,r
        integer flag,scol
        
        integer flag1,flagn,flag3,flagm,flagp
        integer,dimension(6) :: flag2
        integer smc
        ncount=0
        smc=6
        do i=1,nAtoms
           do j=i+1,nAtoms
              do k=j+1,nAtoms
                 do l=k+1,nAtoms
                    do m=l+1,nAtoms
                       do n=m+1,nAtoms
                a(1)=i
                a(2)=j
                a(3)=k
                a(4)=l
                a(5)=m
                a(6)=n
              
                do p=1,6
                   do q=1,6
                   new_mat(p,q)=con_mat(a(p),a(q))
                   end do
                end do

                flag=0
                do p=1,6
                scol=sum(new_mat(p,1:6))
                if(scol.ne.2) flag=1
                end do

                if(flag==0) then                 !if-1
                   flag1=0
                   flagn=0
                   do p=1,6                 !do-1
                  ix(1)=p
                   t=1
                
                      do q=1,6
                       if(new_mat(p,q).ne.0) then  !if-2
                        t=t+1
                        ix(t)=q
                      end if                  !close if-2
                      end do

                      do q=1,3
                         do r=1,3
                         upd_mat(q,r)=new_mat(ix(q),ix(r))
                         end do
                      end do

                      do q=1,3
                      scl(q)=sum(upd_mat(q,1:3))
                      end do

                      if((scl(1)==2.and.scl(2)==1.and.scl(3)==1).or. &
        &               (scl(1)==1.and.scl(2)==2.and.scl(3)==1).or.   &
        &                (scl(1)==1.and.scl(2)==1.and.scl(3)==2)) then
                         flag1=flag1+1
                      end if

                      if(flag1==6) then
                      call centroidcheck(a,filename,nAtoms,smc,flagn)
                        if(flagn==0) then
                        call vectorcheck(a,filename,nAtoms, & 
                                         flagm,smc)
                        if(flagm==0) then
                   write(57,*) i,j,k,l,m,n
                   ncount=ncount+1
                   else
                   write(58,*) i-1,j-1,k-1,l-1,m-1,n-1
                   end if
                   end if
                   end if

                      end do !close do-1
              end if  !close-if-1

                

                  end do
                end do
               end do
              end do  
           end do
         end do

        end
