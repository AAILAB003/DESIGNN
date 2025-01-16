       subroutine vectorcheck(a,filename,nAtoms,flagm,smc)
       implicit none
       integer, dimension(smc) :: a
       character(len=50) :: filename
       real, dimension(nAtoms) :: x,y,z
       character(len=2),dimension(nAtoms) :: atms

       real, dimension(smc) :: dist
       integer, dimension(smc) :: ainx
       integer flagm,nAtoms,i,j,p,flagp,smc

       real vecx1,vecy1,vecz1
       real vecx2,vecy2,vecz2
       real vecx,vecy,vecz
       real nvecx,nvecy,nvecz
       real vecl,factor,rij
       real vrotx,vroty,vrotz

       real rij1,rij2,tmpd,Rc
       real av,bv,cv,arad
       integer tmpi
       real projxa1,projya1,projza1
       real projxa2,projya2,projza2
       real projxb1,projyb1,projzb1
       real projxb2,projyb2,projzb2
       real minxa1,minya1,minza1
       real minxa2,minya2,minza2
       real minxb1,minyb1,minzb1
       real minxb2,minyb2,minzb2
       real nminx1,nminy1,nminz1
       real nminx2,nminy2,nminz2
       real comx,comy,comz
       integer f
       

       call file_reading(atms,x,y,z,nAtoms,filename)

       comx=sum(x(1:nAtoms))/dble(nAtoms)
       comy=sum(y(1:nAtoms))/dble(nAtoms)
       comz=sum(z(1:nAtoms))/dble(nAtoms)
       
        f=0
       write(454,*) nAtoms+2+2
        write(454,*) " "
        do j=1,nAtoms
           f=0
            do i=1,smc
            if(j==a(i)) f=1
            end do
            if(f==0)  write(454,*) atms(j),x(j),y(j),z(j)
            if(f==1)  write(454,*) "Li",x(j),y(j),z(j)
        end do

       
       !--------------------------------------------------------------!
       !               Found the centroid of a simplex                !
       !--------------------------------------------------------------!

       vrotx=sum(x(a(1:smc)))/dble(smc)
       vroty=sum(y(a(1:smc)))/dble(smc)
       vrotz=sum(z(a(1:smc)))/dble(smc)

        write(454,*) "X",vrotx,vroty,vrotz
       !--------------------------------------------------------------!
       !                     Ordering of vectors                      !
       !--------------------------------------------------------------!

       ! Determine the distance between centroid and index

       do i=1,smc
       ainx(i) = a(i)
       dist(i) = sqrt((x(a(i))-vrotx)**2 + (y(a(i))-vroty)**2 +      &
     &                (z(a(i))-vrotz)**2)
       end do

       ! Arranging in the ascending order of distances

       do i=1,smc
          do j=i+1,smc
          if(dist(i)>dist(j)) then

          tmpd=dist(j)
          dist(j)=dist(i)
          dist(i)=tmpd

          tmpi=ainx(j)
          ainx(j)=ainx(i)
          ainx(i)=tmpi

          end if
          end do
      end do

      !---------------------------------------------------------------!
      !       Finding four projected points for occ/vac check         !
      !---------------------------------------------------------------!

      factor=2.8

      ! For the vectors made of shortest lengths
       vecx1=vrotx-x(ainx(1))
       vecy1=vroty-y(ainx(1))
       vecz1=vrotz-z(ainx(1))

       vecx2=vrotx-x(ainx(2))
       vecy2=vroty-y(ainx(2))
       vecz2=vrotz-z(ainx(2))

       av=sqrt((vecx1**2)+(vecy1**2)+(vecz1**2))
       bv=sqrt((vecx2**2)+(vecy2**2)+(vecz2**2))
       cv=sqrt((x(ainx(1))-x(ainx(2)))**2 + (y(ainx(1))-y(ainx(2)))**2+&
     &         (z(ainx(1))-z(ainx(2)))**2)

       arad=(av**2) + (bv**2) - (cv**2)

       if(arad<0.01) then
       vecx2=vrotx-x(ainx(3))
       vecy2=vroty-y(ainx(3))
       vecz2=vrotz-z(ainx(3))
       end if

       vecx=(vecy1*vecz2)-(vecy2*vecz1)
       vecy=-((vecx1*vecz2)-(vecz1*vecx2))
       vecz=(vecx1*vecy2)-(vecy1*vecx2)

       projxa1=vrotx+vecx
       projya1=vroty+vecy
       projza1=vrotz+vecz

       projxa2=vrotx-vecx
       projya2=vroty-vecy
       projza2=vrotz-vecz

       vecl=sqrt((projxa1-vrotx)**2 + (projya1-vroty)**2 +      &
     &           (projza1-vrotz)**2)

       minxa1=factor*((projxa1-vrotx)/vecl) + vrotx
       minya1=factor*((projya1-vroty)/vecl) + vroty
       minza1=factor*((projza1-vrotz)/vecl) + vrotz

       minxa2=factor*((projxa2-vrotx)/vecl) + vrotx
       minya2=factor*((projya2-vroty)/vecl) + vroty
       minza2=factor*((projza2-vrotz)/vecl) + vrotz

      ! For the vectors made of longest lengths

       vecx1=vrotx-x(ainx(smc))
       vecy1=vroty-y(ainx(smc))
       vecz1=vrotz-z(ainx(smc))

       vecx2=vrotx-x(ainx(smc-1))
       vecy2=vroty-y(ainx(smc-1))
       vecz2=vrotz-z(ainx(smc-1))

       av=sqrt((vecx1**2)+(vecy1**2)+(vecz1**2))
       bv=sqrt((vecx2**2)+(vecy2**2)+(vecz2**2))
       cv=sqrt((x(ainx(smc))-x(ainx(smc-1)))**2 + &
     &         (y(ainx(smc))-y(ainx(smc-1)))**2 +&
     &         (z(ainx(smc))-z(ainx(smc-1)))**2)

       arad=((av**2) + (bv**2) - (cv**2))
       if(arad<0.01) then
       vecx2=vrotx-x(ainx(smc-2))
       vecy2=vroty-y(ainx(smc-2))
       vecz2=vrotz-z(ainx(smc-2))
       end if

       vecx=(vecy1*vecz2)-(vecy2*vecz1)
       vecy=-((vecx1*vecz2)-(vecz1*vecx2))
       vecz=(vecx1*vecy2)-(vecy1*vecx2)

       projxb1=vrotx+vecx
       projyb1=vroty+vecy
       projzb1=vrotz+vecz

       projxb2=vrotx-vecx
       projyb2=vroty-vecy
       projzb2=vrotz-vecz

       vecl=sqrt((projxb1-vrotx)**2 + (projyb1-vroty)**2 +      &
     &           (projzb1-vrotz)**2)

       minxb1=factor*((projxb1-vrotx)/vecl) + vrotx
       minyb1=factor*((projyb1-vroty)/vecl) + vroty
       minzb1=factor*((projzb1-vrotz)/vecl) + vrotz

       minxb2=factor*((projxb2-vrotx)/vecl) + vrotx
       minyb2=factor*((projyb2-vroty)/vecl) + vroty
       minzb2=factor*((projzb2-vrotz)/vecl) + vrotz

       !--------------------------------------------------------------!
       ! Obtained four projected points :                             !
       !                                                              !
       !(i)   minxa1,minya1,minza1                                    !
       !(ii)  minxa2,minya2,minza2                                    !
       !(iii) minxb1,minyb1,minzb1                                    !
       !(iv)  minxb2,minyb2,minzb2                                    !
       !                                                              !
       !----------------------END of PROJECTIONS----------------------!

       
       ! Finding the two sets of vectors to average

       !-------Set One----->

       rij1=sqrt((minxa1-minxb1)**2+(minya1-minyb1)**2+ &
     & (minza1-minzb1)**2)
       rij2=sqrt((minxa1-minxb2)**2+(minya1-minyb2)**2+ &
     & (minza1-minzb2)**2)

       if(rij1<=rij2) then

       nminx1=(minxa1+minxb1)/2.0
       nminy1=(minya1+minyb1)/2.0
       nminz1=(minza1+minzb1)/2.0
       
       else if(rij2<=rij1) then

       nminx1=(minxa1+minxb2)/2.0
       nminy1=(minya1+minyb2)/2.0
       nminz1=(minza1+minzb2)/2.0

       end if

       !-------Set Two----->

        rij1=sqrt((minxa2-minxb1)**2+(minya2-minyb1)**2+ &
     & (minza2-minzb1)**2)
       rij2=sqrt((minxa2-minxb2)**2+(minya2-minyb2)**2+ &
     & (minza2-minzb2)**2)

       if(rij1<=rij2) then

       nminx2=(minxa2+minxb1)/2.0
       nminy2=(minya2+minyb1)/2.0
       nminz2=(minza2+minzb1)/2.0

       else if(rij2<=rij1) then

       nminx2=(minxa2+minxb2)/2.0
       nminy2=(minya2+minyb2)/2.0
       nminz2=(minza2+minzb2)/2.0

       end if
 
       write(454,*) "H",nminx1,nminy1,nminz1
       write(454,*) "H",nminx2,nminy2,nminz2

       !--------------------------------------------------------------!
       !--------------------------------------------------------------!
       !        Checking on both sides of the face for vacancy        !
       !--------------------------------------------------------------!
       
       !``````````````````````````````````````````````````````````````!
       !        if flagm or flagp is 0,  it implies a vacancy         !
       !``````````````````````````````````````````````````````````````!

!        Rc=sqrt((factor**2)+(dist(smc)**2))-factor
        Rc=sqrt((factor**2)+(dist(smc)**2))-1.2
!       Rc=sqrt(2*(factor**2))

        rij1=sqrt((comx-nminx1)**2+(comy-nminy1)**2+(comz-nminz1)**2)
        rij2=sqrt((comx-nminx2)**2+(comy-nminy2)**2+(comz-nminz2)**2)

       ! BOTH POINTS NEED NOT BE CHECKED

       if(rij1>=rij2) then
       flagm=0
       do i=1,nAtoms
       f=0
          do j=1,smc
          if(i==a(j)) f=1
          end do
          if(f==0) then
       rij=sqrt((nminx1-x(i))**2+(nminy1-y(i))**2+(nminz1-z(i))**2)
          if(rij<=Rc) flagm=1
          end if
       end do
       write(454,*) "He",nminx1,nminy1,nminz1

       else if(rij2>rij1) then
       flagm=0
       do i=1,nAtoms
       f=0
          do j=1,smc
          if(i==a(j)) f=1
          end do
          if(f==0) then
       rij=sqrt((nminx2-x(i))**2+(nminy2-y(i))**2+(nminz2-z(i))**2)
          if(rij<=Rc) flagm=1
        end if
       end do

       write(454,*) "He",nminx2,nminy2,nminz2

       end if

       end

