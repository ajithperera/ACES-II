       subroutine equi(comp,ncount,comppov,e)

       implicit none
#include <aces.h>
#include <iuhf.com>
#include <mol.com>
#include <grid.com>
#include <kscore.com>

c        M_REAL
c     &  CA(totprim,naobasfn)

         integer ref, k,i,j,comp(natoms),comppov(natoms)

       integer ncount, e(natoms,natoms)

       i=1
       do j=1,ncount 
          ref=kscore(polist+j-1)
CSSS          write(*,*) 'ref=',ref
         do k=0,comppov(j)-1
          e(ref,ref+k)=comp(i)
          i=i+1
         end do
       end do
   
      return
      end
 
