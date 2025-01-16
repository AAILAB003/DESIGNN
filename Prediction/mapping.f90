        subroutine mapping(anAtoms,nAtoms,xa,ya,za,x,y,z, &
     &             atms,aatms,shell_file,full_file)
!       program designn
        implicit none
        real,dimension(nAtoms)     :: x,y,z
        real,dimension(anAtoms)    :: xa,ya,za
        character(len=2),dimension(nAtoms) :: atms
        character(len=2),dimension(anAtoms) ::aatms
        character(len=50)   :: full_file,shell_file
        character(len=50)   :: tempfile,iflname
        integer nAtoms,anAtoms

        integer i,j,k
        real, allocatable,dimension(:)     :: G
        integer, allocatable,dimension(:)  :: G_atms
        real :: fc,eta,Rc,Rs,rij,pi

        real rik,rjk
        integer nfaces
        integer, allocatable,dimension(:)  :: Ai,Aj,Ak,Bk
        real, allocatable,dimension(:)     :: Gface

        real, dimension(3)                 :: r
        integer, dimension(3,3)            :: indx
        real temp
        integer temp1,temp2,temp3

        integer p
        real,dimension(2)                  :: Gi
        integer,dimension(3)               :: ind
        integer atom_i,atom_j,atom_k

        real, allocatable, dimension(:)    :: rp, tp, tp_r
        real ax,ay,az,bx,by,bz,cx,cy,cz
        real a,b,c,alpha_r,alpha_d
        real kd,kx,ky,kz
        real vx,vy,vz,vx_1,vy_1,vz_1
        real kvx,kvy,kvz
        real d_kv,kx_d_kv,ky_d_kv,kz_d_kv
        real, allocatable, dimension(:) :: vx_rot,vy_rot,vz_rot
        real vecx1,vecy1,vecz1,vecx2,vecy2,vecz2
        real proj_x,proj_y,proj_z
        real dist,factor,check
        integer flag
        real, allocatable, dimension(:) :: min_x,min_y,min_z
        real, allocatable, dimension(:) :: new_min_x,new_min_y,new_min_z
        real, allocatable, dimension(:) :: new_vx_rot,new_vy_rot 
        real, allocatable, dimension(:) :: new_vz_rot
        real temp_count
        real temp_minx,temp_miny,temp_minz
        real temp_vxrot,temp_vyrot,temp_vzrot 

        integer nflag

        real fcut
        integer ct
        real, allocatable, dimension(:) :: arij
        integer pos1,pos2
        pi=3.142
        allocate(G(nAtoms),G_atms(nAtoms))
        
        open(1992,file='ESP_Framework_Predicted.xyz')
!       open(39,file='data.txt')
        factor=2.5

        !=============================================================!
        !        File preparation for ESP Framework Deduction         !
        !=============================================================!

!        write(*,*) "INSIDE DESIGNN MODULE--->"
!        write(*,*) "Number of atoms in full:",anAtoms
!        write(*,*) "Number of atoms in shell:",nAtoms

        !---------Only for the faces on the metal shell---------------!

!         allocate(arij((nAtoms*(nAtoms-1))/2))
!        write(*,*) (nAtoms*(nAtoms-1))/2
!        ct=0
!         do i=1,nAtoms
!            do j=i+1,nAtoms
!               rij=sqrt((x(i)-x(j))**2 + (y(i)-y(j))**2 +(z(i)-z(j))**2)
!               ct=ct+1
!               arij(ct)=rij
!            end do
!        end do

!        fcut=(minval(arij)+maxval(arij))/2
!        write(*,*) "THE FACE CUTOFF is:",fcut
        
        !-------------------------------------------------------------!
        !---Determining the face-index list upto six atoms------------!
        !-------------------------------------------------------------!

        call faceindexfinder(shell_file,anAtoms)
        iflname='Face-Indexes'  
        call system("cp Face-Indexes face-index.txt")
        call feature_calculator(anAtoms,full_file,iflname)

        !-------------------------------------------------------------!
        !    Determining the shortest bond length in the 3-atom face  !
        !-------------------------------------------------------------!
        !       Rearrange the face list for Rodriguez Rotation        !
        !-------------------------------------------------------------!


        tempfile='faces'
        call line_counter(tempfile,nfaces)
        allocate(Ai(nfaces),Aj(nfaces),Ak(nfaces))

        open(38,file='faces')
        do i=1,nfaces
        read(38,*) Ai(i),Aj(i),Ak(i)
        end do
        close(38)
        
        do i=1,nfaces

           r(1)=sqrt((x(Ai(i))-x(Aj(i)))**2 +   &
     &               (y(Ai(i))-y(Aj(i)))**2 +   &
     &               (z(Ai(i))-z(Aj(i)))**2)
           indx(1,1)=Ai(i)
           indx(1,2)=Aj(i)

           r(2)=sqrt((x(Ai(i))-x(Ak(i)))**2 +   &
     &               (y(Ai(i))-y(Ak(i)))**2 +   &
     &               (z(Ai(i))-z(Ak(i)))**2)                       

           indx(2,1)=Ai(i)
           indx(2,2)=Ak(i)

           r(3)=sqrt((x(Aj(i))-x(Ak(i)))**2 +   &
     &               (y(Aj(i))-y(Ak(i)))**2 +   &
     &               (z(Aj(i))-z(Ak(i)))**2)                       

           indx(3,1)=Aj(i)
           indx(3,2)=Ak(i)


          do j=1,3
           do k=j+1,3
           if(r(j)>=r(k)) then

             temp=r(j)
             r(j)=r(k)
             r(k)=temp

             temp1=indx(j,1)
             indx(j,1)=indx(k,1)
             indx(k,1)=temp1

             temp2=indx(j,2)
             indx(j,2)=indx(k,2)
             indx(k,2)=temp2
        
           end if
           end do
           end do 
        
          do j=1,2
             do k=1,2
          if(indx(1,j)==indx(2,k)) then
          ind(1)=indx(1,j)
          pos1=j
          pos2=k
          end if
          end do
          end do

          do j=1,2
          if(pos1.ne.j) ind(2)=indx(1,j)
          if(pos2.ne.j) ind(3)=indx(2,j)
          end do

           atom_i=ind(1)
           atom_j=ind(2)
           atom_k=ind(3)
       
      Ai(i)=atom_i
      Aj(i)=atom_j
      Ak(i)=atom_k
        end do

!        close(39)

!        write(*,*)"===================================================="
!        write(*,*)"    PREPARATION OF PREDICTION FILES ---> DONE       "
!        write(*,*)"===================================================="

        !-----------------------------------------------------------------!
        ! (a) Deep-Network training of three atom faces for ESP Framework !
        ! (b) Prediction of the three atom systems in file                ! 
        !-----------------------------------------------------------------!

!        write(*,*)"===================================================="
!        write(*,*)" DEEP LEARNING MODEL FOR TRAINING R & THETA : START "
!        write(*,*)"===================================================="

         call system("python3.9 prediction-model-DESIGNN-R.py")
         call system("python3.9 prediction-model-DESIGNN-T.py")

         call system('paste -d " " predicted_r predicted_t > pred_file')
         call system("sed -i 's/\[/ /g' pred_file")
         call system("sed -i 's/tensor/ /g' pred_file")
         call system("sed -i 's/grad_fn=<AddmmBackward0>//g' pred_file")
         call system("sed -i 's/]/ /g' pred_file")
         call system("sed -i 's/(/ /g' pred_file")
         call system("sed -i 's/)/ /g' pred_file")
         call system("sed -i 's/,/ /g' pred_file")


!        write(*,*)"===================================================="
!        write(*,*)" END OF DEEP LEARNING TRAINNG AND PREDICTION        "
!        write(*,*)"===================================================="
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
        !     POST-PROCESSING of the output file for final framework      !
        !=================================================================!
        !        RODRIGUEZ FORMULA FOR ROTATION IN AN ARBITRARY AXIS      !
        !=================================================================!
        !        vrot = vcos(t) + (k x v)sin(t) + k(k.v)(1-cos(t))        !
        !-----------------------------------------------------------------!        
        
        allocate(rp(nfaces),tp(nfaces),tp_r(nfaces))
        allocate(min_x(nfaces),min_y(nfaces),min_z(nfaces))
        allocate(vx_rot(nfaces),vy_rot(nfaces),vz_rot(nfaces))

        allocate(new_min_x(nfaces),new_min_y(nfaces),new_min_z(nfaces))
        allocate(new_vx_rot(nfaces),new_vy_rot(nfaces))
        allocate(new_vz_rot(nfaces))

        open(40,file='pred_file')
        do i=1,nfaces
        read(40,*) rp(i),tp(i)
        end do
        close(40)

        do i=1,nfaces

        write(350,*) "4"
        write(350,*) "Face no index:",Ai(i),Aj(i),Ak(i)
        write(350,*) "Mg",x(Ai(i)),y(Ai(i)),z(Ai(i))
        write(350,*) "Mg",x(Aj(i)),y(Aj(i)),z(Aj(i))
        write(350,*) "Mg",x(Ak(i)),y(Ak(i)),z(Ak(i))

        ax=x(Ai(i))-x(Aj(i))
        ay=y(Ai(i))-y(Aj(i))
        az=z(Ai(i))-z(Aj(i))

        bx=x(Ai(i))-x(Ak(i))
        by=y(Ai(i))-y(Ak(i))
        bz=z(Ai(i))-z(Ak(i))

        cx=x(Aj(i))-x(Ak(i))
        cy=y(Aj(i))-y(Ak(i))
        cz=z(Aj(i))-z(Ak(i))
 
        a=sqrt( ax**2 + ay**2 + az**2)
        b=sqrt( bx**2 + by**2 + bz**2)
        c=sqrt( cx**2 + cy**2 + cz**2)

        !-----------------------------------------------------------------!
        !        Evaluation of individual terms in the Formulation        !
        !-----------------------------------------------------------------!                

        alpha_r=acos(((a**2)+(b**2)-(c**2))/(2*a*b))
        alpha_d=acos(((a**2)+(b**2)-(c**2))/(2*a*b)) * 57.2958
        kd = a*b*sin(alpha_r)
        kx= (  (ay*bz - az*by) )/kd
        ky= ( -(ax*bz - az*bx) )/kd
        kz= (  (ax*by - ay*bx) )/kd

        vx=(ax/a)*rp(i) 
        vy=(ay/a)*rp(i) 
        vz=(az/a)*rp(i) 


        ! Evaluation of vcos(t)           ===>
        
        tp_r(i) = tp(i)*(pi/180.0)
       
        vx_1=vx*cos(tp_r(i))
        vy_1=vy*cos(tp_r(i))
        vz_1=vz*cos(tp_r(i))


        ! Evaluation of (k x v)sin(t)     ===>

        kvx= (  (ky*vz - kz*vy))  * sin(tp_r(i))
        kvy= ( -(kx*vz - kz*vx))  * sin(tp_r(i))
        kvz= (  (kx*vy - ky*vx))  * sin(tp_r(i))

        ! Evaluation of k(k.v)(1-cos(t))  ===>
        d_kv = (kx*vx + ky*vy + kz*vz)*(1.0-cos(tp_r(i)))
        kx_d_kv = kx * d_kv
        ky_d_kv = ky * d_kv
        kz_d_kv = kz * d_kv
       
        ! Addition of terms               ===> 
        vx_rot(i) =  (vx_1 +  kvx + kx_d_kv)
        vy_rot(i) =  (vy_1 +  kvy + ky_d_kv)
        vz_rot(i) =  (vz_1 +  kvz + kz_d_kv)
        write(350,*) "Li",vx_rot(i),vy_rot(i),vz_rot(i)

        vx_rot(i) = x(Ai(i)) - (vx_1 +  kvx + kx_d_kv)
        vy_rot(i) = y(Ai(i)) - (vy_1 +  kvy + ky_d_kv)
        vz_rot(i) = z(Ai(i)) - (vz_1 +  kvz + kz_d_kv)
        write(350,*) "Be",vx_rot(i),vy_rot(i),vz_rot(i) 
        
        !----------------------------------------------------------!
        !         2D-3D Projection of the Predicted Point          !
        !----------------------------------------------------------!

        vecx1=vx_rot(i)-x(Ai(i))
        vecy1=vy_rot(i)-y(Ai(i))
        vecz1=vz_rot(i)-z(Ai(i))

        vecx2=vx_rot(i)-x(Aj(i))
        vecy2=vy_rot(i)-y(Aj(i))
        vecz2=vz_rot(i)-z(Aj(i))

        proj_x= vx_rot(i) + (  (vecy1*vecz2 - vecz1*vecy2))
        proj_y= vy_rot(i) + ( -(vecx1*vecz2 - vecz1*vecx2))
        proj_z= vz_rot(i) + (  (vecx1*vecy2 - vecy1*vecx2))

        dist=sqrt((proj_x-vx_rot(i))**2 + (proj_y-vy_rot(i))**2 +      &
     &            (proj_z-vz_rot(i))**2)
       
        min_x(i)=factor*((proj_x-vx_rot(i))/dist) + vx_rot(i)
        min_y(i)=factor*((proj_y-vy_rot(i))/dist) + vy_rot(i)
        min_z(i)=factor*((proj_z-vz_rot(i))/dist) + vz_rot(i)
        write(350,*) "X",min_x(i),min_y(i),min_z(i)
        flag=0
        do j=1,anAtoms
        check=sqrt((xa(j)-min_x(i))**2+(ya(j)-min_y(i))**2+     &
     &             (za(j)-min_z(i))**2)
        if(check<2.0) then
          flag=1
        end if
        end do

        if(flag==1) then
        
        vecx1=vx_rot(i)-x(Aj(i))
        vecy1=vy_rot(i)-y(Aj(i))
        vecz1=vz_rot(i)-z(Aj(i))

        vecx2=vx_rot(i)-x(Ai(i))
        vecy2=vy_rot(i)-y(Ai(i))
        vecz2=vz_rot(i)-z(Ai(i))

        proj_x= vx_rot(i) + (  (vecy1*vecz2 - vecz1*vecy2))
        proj_y= vy_rot(i) + ( -(vecx1*vecz2 - vecz1*vecx2))
        proj_z= vz_rot(i) + (  (vecx1*vecy2 - vecy1*vecx2))

        dist=sqrt((proj_x-vx_rot(i))**2 + (proj_y-vy_rot(i))**2 +      &
     &            (proj_z-vz_rot(i))**2)

        min_x(i)=factor*((proj_x-vx_rot(i))/dist) + vx_rot(i)
        min_y(i)=factor*((proj_y-vy_rot(i))/dist) + vy_rot(i)
        min_z(i)=factor*((proj_z-vz_rot(i))/dist) + vz_rot(i)
        end if
        
        write(350,*)  "X",min_x(i),min_y(i),min_z(i)
        end do

        ! if points are still there remove
        !----------------------------------------------------------------!
        !           Averaging the closely predicted ESP Minima           !
        !----------------------------------------------------------------!           

!        write(*,*) "Removing average points started"
        nflag=1
        do while(nflag>=1) 
            p=0
            i=1
            do while(i<=nfaces) 
            if(min_x(i).ne.1000.0.and.min_y(i).ne.1000.0.and.          &
     &         min_z(i).ne.1000.0) then
               write(200,*) "Combination of actual face no.",i
        temp_count=1.0
        temp_minx=min_x(i)
        temp_miny=min_y(i)
        temp_minz=min_z(i)

        temp_vxrot=vx_rot(i)
        temp_vyrot=vy_rot(i)
        temp_vzrot=vz_rot(i)

           do j=i+1,nfaces

           rij=sqrt((min_x(i)-min_x(j))**2 + (min_y(i)-min_y(j))**2 + &
     &              (min_z(i)-min_z(j))**2)
           if(rij<=1.5) then
           write(200,*) j
           temp_count=temp_count+1.0
       
           temp_minx=temp_minx+min_x(j)
           temp_miny=temp_miny+min_y(j)
           temp_minz=temp_minz+min_z(j)

           temp_vxrot=temp_vxrot+vx_rot(j)
           temp_vyrot=temp_vyrot+vy_rot(j)
           temp_vzrot=temp_vzrot+vz_rot(j)

           min_x(j)=1000.00
           min_y(j)=1000.00
           min_z(j)=1000.00

           end if
           end do
           p=p+1
           i=i+1
           new_min_x(p) = temp_minx/temp_count
           new_min_y(p) = temp_miny/temp_count
           new_min_z(p) = temp_minz/temp_count

           new_vx_rot(p) = temp_vxrot/temp_count
           new_vy_rot(p) = temp_vyrot/temp_count
           new_vz_rot(p) = temp_vzrot/temp_count


           else
           i=i+1
           end if
           end do


           do i=1,p
           min_x(i)=new_min_x(i)
           min_y(i)=new_min_y(i)
           min_z(i)=new_min_z(i)

           vx_rot(i)=new_vx_rot(i)
           vy_rot(i)=new_vy_rot(i)
           vz_rot(i)=new_vz_rot(i)

           end do
           nfaces=p

           nflag=0
        do i=1,p
           do j=i+1,p
           rij=sqrt((min_x(i)-min_x(j))**2 + (min_y(i)-min_y(j))**2 + &
     &              (min_z(i)-min_z(j))**2)
           if(rij<=1.5) then
           nflag=1
           end if
           end do
        end do

       end do! flag do while 


        do i=1,anAtoms
        write(1992,*) aatms(i),xa(i),ya(i),za(i)
        end do

           do j=1,p
           flag=0
              do i=1,anAtoms
              rij=sqrt((min_x(j)-xa(i))**2 +   &
     &                 (min_y(j)-ya(i))**2 +   &
     &                 (min_z(j)-za(i))**2)
              if(rij<1.5) then        
              flag=1
              end if
              end do
           if(flag.ne.1) then
             write(1992,*) "X",min_x(j),min_y(j),min_z(j)
             write(198,*) "draw arrow {",vx_rot(j),vy_rot(j),vz_rot(j),&
     &                     "}"," {",min_x(j),min_y(j),min_z(j),"}"

           end if
            end do
        
        write(198,*) " "
        close(1992)


        end
