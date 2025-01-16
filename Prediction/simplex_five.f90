        subroutine simplex_five(nAtoms,con_mat,ncount,filename)

        implicit none
        integer,dimension(nAtoms,nAtoms) :: con_mat
        integer,dimension(5) :: a,flag2
        integer,dimension(3) :: scl,ix,tempa
        integer,dimension(5,5) :: new_mat
        integer,dimension(5,5) :: upd_mat
        character(len=50) :: filename

        integer nAtoms,ncount
                
        integer i,j,k,l,m,n
        integer p,q,t,r
        integer flag,scol
        
        integer flag1,flag3
        integer flagn,smc,flagm,flagp

        ncount=0
        smc=5
        do i=1,nAtoms
           do j=i+1,nAtoms
              do k=j+1,nAtoms
                 do l=k+1,nAtoms
                    do m=l+1,nAtoms
                a(1)=i
                a(2)=j
                a(3)=k
                a(4)=l
                a(5)=m
              
                do p=1,5
                   do q=1,5
                   new_mat(p,q)=con_mat(a(p),a(q))
                   end do
                end do

                flag=0
                do p=1,5
                scol=sum(new_mat(p,1:5))
                if(scol.ne.2) flag=1
                end do

                if(flag==0) then                 !if-1
                        flag1=0
                        do p=1,5                 !do-1
                      ix(1)=p
                      t=1
                      do q=1,5
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
        &                (scl(1)==1.and.scl(2)==2.and.scl(3)==1).or.   &
        &                (scl(1)==1.and.scl(2)==1.and.scl(3)==2)) then
                         flag1=flag1+1
                      end if

                     if(flag1==5) then
                       call centroidcheck(a,filename,nAtoms,smc,flagn)
                          if(flagn==0) then
                        call vectorcheck(a,filename,nAtoms, &
       &                                 flagm,smc)
                        if(flagm==0) then
                          write(57,*) i,j,k,l,m,0
                          ncount=ncount+1
                          else
                          write(58,*) i-1,j-1,k-1,l-1,m-1
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

        end
