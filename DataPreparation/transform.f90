        program transform
        implicit none
        real,allocatable,dimension(:) :: x,y,z,Gi,Gj,Gk
        integer,allocatable,dimension(:) :: Ai,Aj,Ak,Am
        character(len=2),allocatable,dimension(:) :: atms
        integer i,tempind,j,k,t
        real cx,cy,cz,rt 
        real vecx12,vecy12,vecz12,vecx13,vecy13,vecz13
        real vec12,vec13,dot,theta,temp
        real, dimension(3) :: r,r_a
        real, dimension(2) :: r_t
        integer, dimension(3) :: indx
        integer indx_1,indx_2,f_indx,s_indx
        real min_r_a,s_bl,r_data
        integer nfaces,nAtoms

        real p1,p2,p3
        real r1,r2,r3
        real q1,q2,q3
        real detmc
        real,allocatable,dimension(:,:) :: mt
        real,allocatable,dimension(:,:) :: m,mmcinv
        real,allocatable,dimension(:,:) :: mmcinvmt
        real,allocatable,dimension(:,:) :: mc, mcinv
        real,allocatable,dimension(:)   :: pq,pr
        real,allocatable,dimension(:,:) :: top, topp

        real a,b,c,theta_check
        character(len=50) :: fname, gname



        fname='special_esp'
        fname=trim(fname)
        gname='faces'
        gname=trim(gname)

        call line_counter(fname,nAtoms)
        call line_counter(gname,nfaces)

        allocate(atms(nAtoms),x(nAtoms),y(nAtoms),z(nAtoms))
        allocate(Ai(nfaces),Aj(nfaces),Ak(nfaces),Am(nfaces))
        allocate(Gi(nfaces),Gj(nfaces),Gk(nfaces))

        open(41,file='special_esp')
        open(42,file='faces')
        open(118,file='Structure_Descriptor')
        open(119,file='R-Property')
        open(121,file='T-Property')

        do i=1,nAtoms
        read(41,*) atms(i),x(i),y(i),z(i)
        end do

        do i=1,nfaces
        read(42,*) Ai(i),Aj(i),Ak(i),Am(i)
        end do

        do t=1,nfaces

        if(Am(t).ne.0) then
        allocate(mt(2,3),m(3,2),mmcinv(3,2))
        allocate(mmcinvmt(3,3),mc(2,2), mcinv(2,2))
        allocate(pq(3),pr(3))
        allocate(top(3,1),topp(3,1))


        !-------------------------------------------------------------!
        !    PROJECTION OF THE MINIMA VECTOR ON THE PLANE OF 3 ATOMS  !
        !-------------------------------------------------------------!
        ! Orienting according to the origin  w.r.t the central atom P !
        !-------------------------------------------------------------!

        r(1)=sqrt((x(Ai(t))-x(Am(t)))**2 + (y(Ai(t))-y(Am(t)))**2 + &
     &            (z(Ai(t))-z(Am(t)))**2)
        indx(1) = Ai(t)

        r(2)=sqrt((x(Aj(t))-x(Am(t)))**2 + (y(Aj(t))-y(Am(t)))**2 + &
     &           (z(Aj(t))-z(Am(t)))**2)
        indx(2) = Aj(t)

        r(3)=sqrt((x(Ak(t))-x(Am(t)))**2 + (y(Ak(t))-y(Am(t)))**2 + &
     &           (z(Ak(t))-z(Am(t)))**2)
        indx(3) = Ak(t)
      
        do j=1,3
           do k=j+1,3
           if(r(j)>=r(k)) then

             temp=r(j)
             r(j)=r(k)
             r(k)=temp

             tempind=indx(j)
             indx(j)=indx(k)
             indx(k)=tempind

           end if
           end do
        end do

        p1=x(indx(1))
        p2=y(indx(1))
        p3=z(indx(1))

        q1=x(indx(2))
        q2=y(indx(2))
        q3=z(indx(2))

        r1=x(indx(3))
        r2=y(indx(3))
        r3=z(indx(3))

        top(1,1)=x(Am(t))
        top(2,1)=y(Am(t))
        top(3,1)=z(Am(t))

        cx=(p1+q1+r1)/3.0
        cy=(p2+q2+r2)/3.0
        cz=(p3+q3+r3)/3.0

        q1=q1-p1
        q2=q2-p2
        q3=q3-p3
      
        r1=r1-p1
        r2=r2-p2
        r3=r3-p3
      
        top(1,1)= top(1,1)-p1
        top(2,1)= top(2,1)-p2
        top(3,1)= top(3,1)-p3
      
        p1=p1-p1
        p2=p2-p2
        p3=p3-p3
     
       ! write(*,*) p1,p2,p3
       ! write(*,*) q1,q2,q3
       ! write(*,*) r1,r2,r3
       ! read(*,*)

        pq(1)=q1-p1
        pq(2)=q2-p2
        pq(3)=q3-p3
      
        pr(1)=r1-p1
        pr(2)=r2-p2
        pr(3)=r3-p3

        m(1,1)=pq(1)
        m(2,1)=pq(2)
        m(3,1)=pq(3)
   
        m(1,2)=pr(1)
        m(2,2)=pr(2)
        m(3,2)=pr(3)
   
   
        mt(1,1)=m(1,1)
        mt(1,2)=m(2,1)
        mt(2,1)=m(1,2)
        mt(2,2)=m(2,2)
        mt(1,3)=m(3,1)
        mt(2,3)=m(3,2)

        do i=1,2
            do j=1,2
            mc(i,j)=0.0
                 do k=1,3
                   mc(i,j)= mc(i,j) + mt(i,k) * m(k,j)
                 end do
           end do
        end do

        detmc= mc(1,1) * mc(2,2) - mc(1,2) * mc(2,1)
      

        mcinv(1,1)=(1.0/detmc)* mc(2,2)
        mcinv(2,2)=(1.0/detmc)* mc(1,1)
        mcinv(1,2)=-(1.0/detmc)* mc(2,1)
        mcinv(2,1)=-(1.0/detmc)* mc(1,2)
      
        do i=1,3
          do j=1,2
              mmcinv(i,j)=0.0
              do k=1,2
                mmcinv(i,j) = mmcinv(i,j) + m(i,k) * mcinv(k,j)
              end do
          end do
        end do
      
        do i=1,3
           do j=1,3
               mmcinvmt(i,j)=0.0
               do k=1,2
                 mmcinvmt(i,j) = mmcinvmt(i,j) + mmcinv(i,k) * mt(k,j)
               end do
           end do
         end do
      
         do i=1,3
             do j=1,1
                 topp(i,j)=0.0
                 do k=1,3
                   topp(i,j) = topp(i,j) + mmcinvmt(i,k) * top(k,j)
                 end do
             end do
           end do
        
           !Minima returned to actual orientiation in space on 3atom
           !plane

           topp(1,1)=topp(1,1)+x(indx(1))
           topp(2,1)=topp(2,1)+y(indx(1))
           topp(3,1)=topp(3,1)+z(indx(1))

        ! Find the shortest bond length between metal atoms
           
           r_a(1)= sqrt((x(Ai(t))-x(Aj(t)))**2 +        &
     &                  (y(Ai(t))-y(Aj(t)))**2 +        &
     &                  (z(Ai(t))-z(Aj(t)))**2)
                   
           r_a(2)= sqrt((x(Ai(t))-x(Ak(t)))**2 +        &
     &                  (y(Ai(t))-y(Ak(t)))**2 +        &
     &                  (z(Ai(t))-z(Ak(t)))**2)

           r_a(3)= sqrt((x(Aj(t))-x(Ak(t)))**2 +        &
     &                  (y(Aj(t))-y(Ak(t)))**2 +        &
     &                  (z(Aj(t))-z(Ak(t)))**2)

          
           min_r_a=minval(r_a)

           if(r_a(1)==minval(r_a)) then
             indx_1=Ai(t)      
             indx_2=Aj(t)
           else if(r_a(2)==minval(r_a)) then
             indx_1=Ai(t)
             indx_2=Ak(t)      
           else if(r_a(3)==minval(r_a)) then        
             indx_1=Aj(t)
             indx_2=Ak(t)
           end if
          
            

       !----------------------------------------------!
       ! find the distances between the plane minima  !
       ! and all the three atoms                      !
       !----------------------------------------------!

       r(1)=sqrt((topp(1,1)-x(Ai(t)))**2 +    &    !P atom
     &           (topp(2,1)-y(Ai(t)))**2 +    &
     &           (topp(3,1)-z(Ai(t)))**2)
       indx(1)=Ai(t)

       r(2)=sqrt((topp(1,1)-x(Aj(t)))**2 +    &    !Q atom
     &           (topp(2,1)-y(Aj(t)))**2 +    &
     &           (topp(3,1)-z(Aj(t)))**2)
       indx(2)=Aj(t)
       r(3)=sqrt((topp(1,1)-x(Ak(t)))**2 +    &    !R atom
     &           (topp(2,1)-y(Ak(t)))**2 +    &
     &           (topp(3,1)-z(Ak(t)))**2)
       indx(3)=Ak(t)



        ! Find the distance between the minima on plane and the shortest
        ! bond length between atoms obtained above

           r_t(1)=sqrt((topp(1,1)-x(indx_1))**2 +    &    
     &                 (topp(2,1)-y(indx_1))**2 +    &
     &                 (topp(3,1)-z(indx_1))**2)

           r_t(2)=sqrt((topp(1,1)-x(indx_2))**2 +    &    
     &                 (topp(2,1)-y(indx_2))**2 +    &
     &                 (topp(3,1)-z(indx_2))**2)

        
     !   write(116,*) minval(r_t),(Gi(t)+Gj(t)+Gk(t))/minval(r_t)
     !   write(116,*) minval(r_t),(Gi(t)+Gj(t)+Gk(t))/minval(r_a)
     !   write(118,'(f10.5,A)',advance='no') &
     ! &             (Gi(t)+Gj(t)+Gk(t))/minval(r_t),","

!        write(118,'(f10.1,A)',advance='no') &
!     &             (Gi(t)+Gj(t)+Gk(t))/minval(r_a),","
        r_data=minval(r_t)

           if(r_t(1)<r_t(2)) then
             f_indx=indx_1
             s_indx=indx_2
             a=r_t(1)
             c=r_t(2)
           else if(r_t(2)<r_t(1)) then
             f_indx=indx_2
             s_indx=indx_1
             a=r_t(2)
             c=r_t(1)
           end if

        ! find the component

!         vecx12=topp(1,1)-x(indx(1))
!         vecy12=topp(2,1)-y(indx(1))
!         vecz12=topp(3,1)-z(indx(1))

        vecx12=x(f_indx)-topp(1,1)
        vecy12=y(f_indx)-topp(2,1)
        vecz12=z(f_indx)-topp(3,1)

        vec12=sqrt(vecx12**2+vecy12**2+vecz12**2)

        vecx13=x(f_indx)-x(s_indx)
        vecy13=y(f_indx)-y(s_indx)
        vecz13=z(f_indx)-z(s_indx)

        vec13=sqrt(vecx13**2+vecy13**2+vecz13**2)
        b=vec13
        dot=((vecx12*vecx13)+(vecy12*vecy13)+(vecz12*vecz13))
        theta=acos(dot/(vec12*vec13)) * 57.2958
        theta_check=acos((a**2+b**2-c**2)/(2*a*b)) * 57.2958

        write(119,*) minval(r_t)
        write(121,*) theta_check

         deallocate(mt,m,mmcinv)
         deallocate(mmcinvmt,mc, mcinv)
         deallocate(pq,pr)
         deallocate(top,topp)

        end if
        end do
        close(119)
        close(121)

        call system("cat R-Property >> r-value")
        call system("cat T-Property >> t-value")
        stop
        end

        



