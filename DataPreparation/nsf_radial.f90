       subroutine nsf_radial(eta,Rs,Rc,x,y,z,atms,nAtoms,&
     & G,F,pair_wt,edgno)

        implicit none

        real, dimension(nAtoms) :: G,G_i        
        real, dimension(nAtoms) :: x,y,z
        real, dimension(edgno) :: pair_wt

        integer, dimension(nAtoms) :: nelem,F
        character(len=2),dimension(nAtoms) :: atms
        real Rc,eta,rij,Rs,fc,pi,Gf
        integer nAtoms,i,j,findex,cjj,ctr
        integer edgno

        pi=3.142
        ctr=0
        open(145,file='indexes-pair')
        open(146,file='sizes-pair')


        do i=1,nAtoms       ! 1-do
           G(i)=0.0
           !nelem(i)=0.0
           if(F(i)==1) then   ! 1-if

           do j=1,nAtoms      ! 2-do

           rij=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
           if (rij>0.0.and.rij<=Rc) then   !2-if
           write(145,*) i-1,j-1
           fc=0.5*(cos((pi*pi*rij)/(180*Rc))+1.0)
           ctr=ctr+1
           pair_wt(ctr)=(exp(-eta*(rij-Rs)**2)*fc)
           !nelem(i)=nelem(i)+1
           else
           fc=0.0
           end if                          !2-end if
           G(i)=G(i)+(exp(-eta*(rij-Rs)**2)*fc) 
           end do   !2-enddo

           !if(nelem(i)>0) then   !3-if
           !G_i(i)=G(i)/dble(nelem(i))
           !else
           !G_i(i)=0.0
           !end if               !3-endif

           G_i(i)=G(i)

           else if(F(i)==0) then  !1-elseif
 
           do j=1,nAtoms      
              if(F(j)==1) then
                 rij=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
                 if (rij>0.0.and.rij<=Rc) then
                 fc=0.5*(cos((pi*pi*rij)/(180*Rc))+1.0)
                 write(145,*) i-1,j-1
                 ctr=ctr+1
                 pair_wt(ctr)=(exp(-eta*(rij-Rs)**2)*fc)
                 !nelem(i)=nelem(i)+1
                 else
                 fc=0.0
                 end if
                 G(i)=G(i)+(exp(-eta*(rij-Rs)**2)*fc)
              end if
           end do

           !if(nelem(i)>0) then
           !G_i(i)=G(i)/dble(nelem(i))
           !else
           !G_i(i)=0.0
           !end if

           G_i(i)=G(i)

          end if
        end do   

        do i=1,nAtoms
        G(i)=G_i(i)
        end do

        write(146,*) ctr
        close(145)
        close(146)

       end


