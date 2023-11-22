      subroutine pack2(indi,indj,ij)
use indices
!
!     this subroutine packs indices i and j into a single integer
!     note that i is STRICTLY <= j 
!
!
      integer i,j,k,l,ij,nlen,indi,indj
!         if(indi.ne.indj)then
         i=min0(indi,indj)
         j=max0(indi,indj)
       ij=i+offset2(j)
       return
      end subroutine pack2

      
