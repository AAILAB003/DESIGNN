        subroutine faceindexfinder(filename,tAtoms)
        implicit none
        character(len=50) :: filename ! the shell file
        integer nAtoms
        !double precision, allocatable,dimension(:) :: x,y,z
        real, allocatable,dimension(:) :: x,y,z
        character(len=2),allocatable,dimension(:) :: atms

        integer, allocatable,dimension(:,:) :: con_mat
        integer,allocatable,dimension(:,:) :: facindx
        integer i,j,k,l,m,p
        integer simth_nct
        integer tAtoms,nfac,tcount
        !take in the total number of atoms in the full system

        integer ix1,ix2,ix3,nix1,nix2,flag,tfac
        real, allocatable,dimension(:,:) :: distmat
        integer, allocatable,dimension(:) :: tempindx
        integer, allocatable,dimension(:,:) :: indxnoter
        character(len=50) :: fname,flname
        real minrij,minrik,minrjk,rij


       !write(*,*) "Received filename:",filename,tAtoms
         
        open(47,file='Face-Indexes')
        open(38,file='faces')

        call line_counter(filename,nAtoms)
        allocate(x(nAtoms),y(nAtoms),z(nAtoms),atms(nAtoms))
        call file_reading(atms,x,y,z,nAtoms,filename)


        allocate(con_mat(nAtoms,nAtoms))
        
        do i=1,nAtoms
           con_mat(i,1:nAtoms)=0
        end do

        call adjacency_matrix(x,y,z,atms,nAtoms,con_mat)

        
        tcount=0
        call simplex_three(nAtoms,con_mat,simth_nct,filename)
        tcount=tcount+simth_nct
        call simplex_four(nAtoms,con_mat,simth_nct,filename)
        tcount=tcount+simth_nct
        call simplex_five(nAtoms,con_mat,simth_nct,filename)
        tcount=tcount+simth_nct
        call simplex_six(nAtoms,con_mat,simth_nct,filename)
        tcount=tcount+simth_nct
        call system("cp fort.57 Face-List")
        nfac=tcount

        allocate(facindx(nfac,6))
        facindx(1:nfac,1:6)=0

        fname='Face-List'
        open(480,file=fname)
        do i=1,nfac
        read(480,*) (facindx(i,j),j=1,6)
        end do
        close(480)

        !tAtoms= total atoms in the full_system

        do i=1,nfac        ! i-do
        do j=1,tAtoms   ! j-do
           flag=0
           do k=1,6     ! k-do
           if(facindx(i,k)==j) then ! 1-if
             flag=1
           end if               ! 1-endif
           end do ! k-enddo

           if(flag==1) then     ! 2-if
            write(47,'(i0,A)',advance='no') 1," "
           else
            write(47,'(i0,A)',advance='no') 0," "
           end if               ! 2-endif
           end do ! end of j
           write(47,*)
           end do ! end of i

        close(47)

        !finding the three atoms which are to be used for projection
        
        flname='Face-Indexes'
        call line_counter(flname,tfac)
        allocate(indxnoter(tfac,tAtoms))
        
        open(47,file='Face-Indexes')
        do i=1,nfac
        read(47,*) (indxnoter(i,j),j=1,tAtoms)
        end do
        close(47)

        do i=1,nfac
           p=0
           allocate(tempindx(6))
           tempindx(1:6) =0
           do j=1,tAtoms
              if(indxnoter(i,j)==1)then
               p=p+1
               tempindx(p)=j
              end if 
           end do
          
          !write(*,*) "Found atoms:",p
          allocate(distmat(p,p))
          distmat(1:p,1:p)=0.0
          do j=1,p
             ix1=tempindx(j)
             do k=j+1,p
             ix2=tempindx(k)
             rij=sqrt((x(ix1)-x(ix2))**2+(y(ix1)-y(ix2))**2+ & 
     &       (z(ix1)-z(ix2))**2)         
             distmat(j,k)=rij
             distmat(k,j)=rij
             end do
          end do

          minrij=10.0
          do j=1,p
             do k=j+1,p
               if(distmat(j,k)<minrij) then
                 ix1=j
                 ix2=k
                 minrij=distmat(j,k)
               end if
             end do
          end do

          minrjk=10.0
          do j=1,p
          if(distmat(ix1,j)<minrjk.and.j.ne.ix2.and.ix1.ne.j) then
             nix1=j
             minrjk=distmat(ix1,j)
          end if
          end do

          minrjk=10.0
          do j=1,p
          if(distmat(ix2,j)<minrjk.and.j.ne.ix1.and.ix2.ne.j) then
             nix2=j
             minrjk=distmat(ix2,j)
          end if
          end do

          if(distmat(ix1,nix1) < distmat(ix2,nix2)) then
             ix3=nix1
          else
             ix3=nix2
          end if
          

          write(38,'(i0,i5,i5)') tempindx(ix1),tempindx(ix2),& 
     & tempindx(ix3)
          deallocate(distmat,tempindx)

        end do
        close(38)

        end

