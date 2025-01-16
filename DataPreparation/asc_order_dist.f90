        subroutine asc_order_dist(lt,li,n,MX)
        real, dimension(n) :: lt
        integer,dimension(1,6) :: MX

        integer, dimension(n) :: li
        integer i,j,p,q
        
        !-------------------------------------------------------------!
        ! Arrange the distances from small to large and interchange   !
        ! the indexes accordingly.                                    !
        !-------------------------------------------------------------!
        p=1
        q=6

        do i=1,n
           do j=i+1,n
           if(lt(i)>lt(j)) then

           tmplt=lt(j)
           lt(j)=lt(i)
           lt(i)=tmplt

           tmpli=li(j)
           li(j)=li(i)
           li(i)=tmpli

           end if
           end do
        end do

        do i=1,3
        MX(p,i)=li(i)
        end do

        i=3
        do while(i<=5)
           j=i+1
           if(abs(lt(i)-lt(j))<0.3) then
            MX(p,j)=li(j)
           i=i+1
           else
           i=6
           end if
        end do


           
        end
                   
