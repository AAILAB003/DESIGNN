        subroutine pair_indices(xt,yt,zt,atms,nAtoms)
        real,dimension(nAtoms) :: xt,yt,zt
        character(len=2),dimension(nAtoms) :: atms,atmst
        integer nAtoms,i,j,nAtm,pcount
        real Rc,rij

        Rc=3.8
        do i=1,nAtoms
           do j=1,nAtoms
           rij=sqrt((xt(i)-xt(j))**2+(yt(i)-yt(j))**2+(zt(i)-zt(j))**2)
           if(rij>0.0.and.rij<=Rc) then
           pcount=pcount+1
           write(455,'(i0,A,i0)') i-1," ",j-1 ! python-indexing 
           end if
           end do
        end do
        write(456,*) pcount

        end





