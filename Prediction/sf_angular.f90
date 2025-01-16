        subroutine sf_angular(atms,x,y,z,nAtoms, &
     &  eta,Rs,zeta,l,Rc,G2,F)

        real, dimension(nAtoms)   :: x,y,z,G2
        integer,dimension(nAtoms) :: F
        character(len=2),dimension(nAtoms) :: atms

        integer i,j,k
        real Gtemp,G
        real vec1,vecx1,vecy1,vecz1
        real vec2,vecx2,vecy2,vecz2
        real dot_pdkx, dot_pdky, dot_pdkz, dot_pdk
        real cos_ijk, ang_rad, deg_theta

        real zeta,l,Rc,Rs,eta,pi,rij,rik,rjk,fc_ij,fc_ik,fc_jk

        integer nfaces, nAtoms

        pi=3.142

        do i=1,nAtoms
           Gtemp=0.0

           do j=1,nAtoms
           vecx1 = x(i)-x(j)
           vecy1 = y(i)-y(j)
           vecz1 = z(i)-z(j)

           vec1  = vecx1**2 + vecy1**2 + vecz1**2
           rij=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
           
           do k=1,nAtoms
           if(F(i)==1.or.F(j)==1.or.F(k)==1) then

           rik=sqrt((x(i)-x(k))**2+(y(i)-y(k))**2+(z(i)-z(k))**2)
           rjk=sqrt((x(j)-x(k))**2+(y(j)-y(k))**2+(z(j)-z(k))**2)

           if(rij<=Rc.and.rik<=Rc.and.rjk<=Rc.and. &
     &       (i.ne.j.and.i.ne.k)) then
           
           fc_ij=0.5*(cos((pi*pi*rij)/(180.0*Rc))+1.0)
           fc_ik=0.5*(cos((pi*pi*rik)/(180.0*Rc))+1.0)
           fc_jk=0.5*(cos((pi*pi*rjk)/(180.0*Rc))+1.0)

!           write(*,*) i,j,k,fc_ij,fc_ik,fc_jk

           vecx2 = x(i)-x(k)
           vecy2 = y(i)-y(k)
           vecz2 = z(i)-z(k)
           vec2  = vecx2**2 + vecy2**2 + vecz2**2

           dot_pdkx= vecx1 * vecx2
           dot_pdky= vecy1 * vecy2
           dot_pdkz= vecz1 * vecz2
           dot_pdk=dot_pdkx+dot_pdky+dot_pdkz

           cos_ijk= dot_pdk/(sqrt(vec1 * vec2))
           ang_rad=acos(cos_ijk)
           deg_theta= (ang_rad*180)/3.142

           Gtemp=Gtemp+(((1.0+(l*cos_ijk))**zeta)              &
     &                   *(exp(-eta*((rij**2)+(rik**2)+(rjk**2)))) &
     &                   *(fc_ij * fc_ik * fc_jk ))

           G=(((1.0+(l*cos_ijk))**zeta)              &
     &                   *(exp(-eta*((rij**2)+(rik**2)+(rjk**2)))) &
     &                   *(fc_ij * fc_ik * fc_jk ))
           write(38,*) deg_theta,G

            end if
            end if
            end do
          end do
          G2(i)=2**(1-zeta)*Gtemp
          end do

          write(38,*) 

        
          end 

