      subroutine pack1(indi,indj,ij)
use indices

!     this subroutine packs indices i and j into a single integer
!     note that i is STRICTLY <= j 
!
!
integer,intent(in)::indi,indj
integer,intent(inout)::ij
integer::i,j
         i=min0(indi,indj)
         j=max0(indi,indj)
!      ij=i+j*(j-1)/2
       ij=i+offset1(j)
       return
      end subroutine pack1

      
