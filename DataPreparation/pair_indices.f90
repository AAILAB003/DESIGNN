        subroutine pair_indices(x,y,z,atms,nAtoms,ctx)
        real,dimension(nAtoms) :: x,y,z
        real,dimension(nAtoms) :: xt,yt,zt
        character(len=2),dimension(nAtoms) :: atms,atmst
        integer nAtoms,ctx,i,j,nAtm,pcount
        real Rc,rij

        Rc=4.2
        j=0

        do i=ctx+1,nAtoms
        j=j+1
        xt(j)=x(i)
        yt(j)=y(i)
        zt(j)=z(i)
        atmst(j)=atms(i)
        end do

        nAtm=j
        pcount=0
        
        open(455,file='Pair-List')
        do i=1,nAtm
           do j=1,nAtm
           rij=sqrt((xt(i)-xt(j))**2+(yt(i)-yt(j))**2+(zt(i)-zt(j))**2)
           if(rij>0.0.and.rij<=Rc) then
           pcount=pcount+1
           write(455,'(i0,A,i0)') i-1," ",j-1 ! python-indexing 
          ! write(*,'(i0,A,i0)') i-1," ",j-1 ! python-indexing 
           end if
           end do
        end do

        open(456,file='nPairList-structure')
        write(456,*) pcount
        close(455)
        close(456)

        call system("cat Pair-List >> Pairs.txt")
        call system("cat nPairList-structure >> nPairs-list.txt")

        end





