        subroutine feature_calculator(nAtoms,flname,iflname)
        implicit none
        character(len=50) :: flname,iflname
        real,allocatable,dimension(:) :: x,y,z,G,Gang,pair_wt
        character(len=2),allocatable,dimension(:) :: atms
        integer,allocatable,dimension(:,:) :: FIndex
        integer,allocatable,dimension(:) :: F
        real,allocatable,dimension(:,:) :: Gf,Gfang,Fpair
        integer ii,jj,i,j,nAtoms,nrows,edno

        character(len=50) :: paramfile
        integer n_rsf,param,flag
        real, allocatable,dimension(:) :: eta,Rs,Rc,eta1,Rs1,Rc1

        integer n_asf,p
        real, allocatable,dimension(:) :: zeta,lamb

        !--------------------------------------------------!
        !         R.Symmetry Function Parameters           !
        !--------------------------------------------------!

        open(199,file='atomic-features.txt')
        open(200,file='Edge-features.txt')
        paramfile='RSF_Parameter_file'
        call line_counter(paramfile,n_rsf)
!        write(*,*) "Number of lines found :",n_rsf
        allocate(eta(n_rsf),Rs(n_rsf),Rc(n_rsf))

        open(56,file=paramfile)
        do i=1,n_rsf
        read(56,*) eta(i),Rs(i),Rc(i)
        end do
        close(56)

        paramfile='ASF_Parameter_file'
        call line_counter(paramfile,n_asf)
!        write(*,*) "Number of lines found :",n_asf
        allocate(eta1(n_asf),Rs1(n_asf),zeta(n_asf),& 
      & lamb(n_asf),Rc1(n_asf))

        open(56,file=paramfile)
        do i=1,n_asf
        read(56,*) eta1(i),Rs1(i),zeta(i),lamb(i),Rc1(i)
        end do
        close(56)

        !--------------------------------------------------!
        !                  Files Reading                   !
        !--------------------------------------------------!

        call line_counter(flname,nAtoms)
        allocate(x(nAtoms),y(nAtoms),z(nAtoms),atms(nAtoms))
!        write(*,*) "The number of atoms found:",nAtoms

        call file_reading(atms,x,y,z,nAtoms,flname)
        
        !--------------------------------------------------!
        !         Files Reading of Face Indexes            !
        !--------------------------------------------------!

        iflname='Face-Indexes'
        iflname=trim(iflname)

        call line_counter(iflname,nrows)
!        write(*,*) "FOUND LINES IN FACE-INDEXES:",nrows,nAtoms
        allocate(FIndex(nrows,nAtoms))
        allocate(F(nAtoms))

        open(58,file=iflname)
        do i=1,nrows
           read(58,*) (FIndex(i,j),j=1,nAtoms)
        end do
        close(58)
        !--------------------------------------------------!
        !         R.Symmetry Function Calculations         !
        !--------------------------------------------------!

        allocate(Gf(n_rsf,nAtoms),G(nAtoms))
        allocate(Gfang(n_asf,nAtoms),Gang(nAtoms))

        edno=nAtoms*nAtoms

        allocate(Fpair(n_rsf,edno))
        FPair(1:n_rsf,1:edno)=0.0

        open(445,file='nAtoms.txt')
        do ii=1,nrows
           write(445,*) nAtoms
           do jj=1,nAtoms
           F(jj)=Findex(ii,jj)
           end do

        edno=nAtoms*nAtoms

        allocate(pair_wt(edno))
        pair_wt(1:edno)=0.0

           do i=1,n_rsf
           call nsf_radial(eta(i),Rs(i),Rc(i),x,y,z,atms, &
     & nAtoms,G,F,pair_wt,edno)

                do j=1,nAtoms
                   Gf(i,j)=G(j)
                end do

                do j=1,edno
                FPair(i,j)=pair_wt(j)
                end do
           end do


           call system("cat indexes-pair >> Edge-Weight-Pairs")
           call system("cat sizes-pair >> Edge-Wt-Sizes-Pairs")
           do i=1,n_asf
           call sf_angular(atms,x,y,z,nAtoms,eta1(i),Rs1(i),zeta(i),& 
     & lamb(i),Rc1(i),Gang,F)
               do j=1,nAtoms
                   Gfang(i,j)=Gang(j)
               end do
           end do

          do j=1,nAtoms
          write(199,*) Gf(1:n_rsf,j),Gfang(1:n_asf,j)
          end do
          
          do j=1,edno
          if(sum(FPair(1:n_rsf,j))> 0.0) then
          write(200,*) FPair(1:n_rsf,j)
          end if
          end do
          deallocate(pair_wt)

          call pair_indices(x,y,z,atms,nAtoms)
          end do
        close(445)

        call system("mv fort.455 Pairs.txt")
        call system("mv fort.456 nPairs-list.txt")

        end
