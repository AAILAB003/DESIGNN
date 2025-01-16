        program topofaceassociate
        implicit none
        real, allocatable,dimension(:) :: x,y,z,d_tp,G
        real, allocatable,dimension(:,:) :: Gf
        integer,allocatable,dimension(:) :: d_tp_idx
        integer,allocatable,dimension(:,:) :: MX, MXtemp

        character(len=2),allocatable,dimension(:) :: atms
        integer nAtoms,i,j,ctx,nMetal,k,index_j,ii,nlines
        real rij
        character(len=50) :: filename,flname,iflname
        integer flag
        character(len=150),allocatable,dimension(:) :: paths


        write(*,*) "Enter the list of filename with path:"
        read(*,*) filename

        call line_counter(filename,nlines)
        write(*,*) "Found lines:" ,nlines
        allocate(paths(nlines))

        open(56,file=filename)
        do i=1,nlines
        read(56,'(A)') paths(i)
        end do
        close(56)


        !-------------------------------------------------------------!
        !    Reading of the variables etc for the structure file      !
        !-------------------------------------------------------------!

        do ii=1,nlines
        filename=trim(paths(ii))
        write(*,*) filename
        call line_counter(filename,nAtoms)
        allocate(x(nAtoms),y(nAtoms),z(nAtoms),atms(nAtoms),G(nAtoms))

        !write(*,*) "The number of atoms found:",nAtoms
        open(55,file=filename)
        do i=1,nAtoms
        read(55,*) atms(i),x(i),y(i),z(i)
        end do

        !-------------------------------------------------------------!
        !       Determine the number ESP Topography Minima            !
        !-------------------------------------------------------------!


        ctx=0
        do i=1,nAtoms
        if(atms(i)=='X') ctx=ctx+1
        end do
        
        !-------------------------------------------------------------!
        ! PREPARATION OF FILE FOR ATOM FEATURES FOR THE SPECIFIC FACE !
        !-------------------------------------------------------------!
        
         nMetal=nAtoms-ctx
         open(37,file='faces') ! USED ONLY TO GET THE TOPOGRAPHY PROJ.
         allocate(d_tp(nMetal),d_tp_idx(nMetal))
         allocate(MX(ctx,6),Mxtemp(1,6))

         open(78,file='Coordinates')
         do j=ctx+1,nAtoms
         write(78,*) atms(j),x(j),y(j),z(j)
         end do
         close(78)

         open(47,file='Face-Inde')

         MX(1:ctx,1:6)=0
         MXtemp(1,1:6)=0
        
         do i=1,ctx
             k=0

             !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
             !------TOPOGRAPHY PROJECTION THROUGH TWO VECTORS-------!
             !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

             do j=ctx+1,nAtoms
             k=k+1
             rij=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
             d_tp(k)=rij
             d_tp_idx(k)=k+ctx
             end do
             call asc_order_dist(d_tp,d_tp_idx,nMetal,MXtemp)

             write(37,*) d_tp_idx(1),d_tp_idx(2),d_tp_idx(3),i

            do k=ctx+1,nAtoms
            flag=0
            
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
            !------DETERMINING THE TYPE OF FACE ENGAGED WITH "X"---!
            ! AND THE FEATURE REPRESENTATION WITH INDEXES OF ATOMS !
            !------------------------------------------------------!

            do j=1,6
             if(Mxtemp(1,j).ne.0.and.MXtemp(1,j)==k) then
               index_j=Mxtemp(1,j)
               flag=1
             end if
             end do ! end of j
            
             if(flag==1) then
                write(47,'(i0,A)',advance='no') 1," "
             else 
                write(47,'(i0,A)',advance='no')  0," "
             end if
             end do ! end of k
             write(47,*) 

         end do

         close(37)
         close(38)
         close(47)

         !---------------------------------------------------------!
         !----REPRESENTATION OF  "X" AS R-VALUE AND THETA-VALUE----!
         !DETERMINATION OF THE INDEXES INVOLVED FOR EACH STRUCTURE !
         !---------------------------------------------------------!
        

         flname='Coordinates'
         iflname='Face-Inde'
         call feature_calculator(nMetal,flname,iflname)

         open(57,file='special_esp')
         do i=1,nAtoms
         write(57,*) atms(i),x(i),y(i),z(i)
         end do
         close(57)

         call system("./transformation.out")
         do i=1,ctx
         call pair_indices(x,y,z,atms,nAtoms,ctx)
         end do
        
         deallocate(x,y,z,atms,G,d_tp,d_tp_idx,Mx,Mxtemp)
         call system("cat Face-Inde >> face-index.txt")

         end do  !! end of ii
         stop
         end
