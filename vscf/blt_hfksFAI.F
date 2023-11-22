      Subroutine blt_hfksFai(C,iuhf,ld3,Fock,scr1,ld2,ld1,
     &        scrtmp,CT,tfia,tfib,
     &        occvrt1,occvrt2,
     &        pop)

      implicit none
#include <aces.h>
#include <symm2.com>
#include <machsp.com>
       integer ld1,ld2,ld3,iuhf
       integer occvrt2,occvrt1
       integer pop(8,2)
       integer indx2

       M_REAL
     &       C((IUHF+1)*LD3),
     &       fock((iuhf+1)*ld1),
     &       scr1(ld2),
     &       scrtmp(ld2),
     &       CT(ld2),
     &       tfia(occvrt1),
     &       tfib(occvrt2)
      integer i,j,n,a,ind,ispin
      INDX2(I,J,N)=I+(J-1)*N

      do ispin=1,iuhf+1
         ind=1
        do i=1,nirrep
c         if (ispin .eq. 1 ) goto 102
         if(nbfirr(i) .eq.0) goto 101
         call expnd2(fock(((ispin-1)*ld1)+itriof(i)),scr1,
     &              nbfirr(i))
         call transp(C(((ispin-1)*ld3)+isqrof(i)),CT,
     &             nbfirr(i),nbfirr(i))
         call mxm(CT,nbfirr(i),scr1,nbfirr(i),scrtmp,nbfirr(i))
         call mxm(scrtmp,nbfirr(i),C(((ispin-1)*ld3)+isqrof(i)),
     &          nbfirr(i),scr1,nbfirr(i))
c       call kkk(nbfirr(i)*nbfirr(i),scr1)
c 102    continue
         do j=1,pop(i,ispin)
           do a=pop(i,ispin)+1,nbfirr(i)
            if(ispin .eq. 1) then
c                tfia(ind)=scr1(j+(a-1)*nbfirr(i))
               tfia(ind)=scr1(indx2(j,a,nbfirr(i)))
            ind=ind+1
           else
               tfib(ind)=scr1(indx2(j,a,nbfirr(i)))
               ind=ind+1
            end if
           end do 
         end do
101    continue
      end do
      end do
      call putrec(20,'JOBARC','HFKSFIAA',occvrt1*iintfp,
     &             tfia)
      if(iuhf .eq. 1) then
      call putrec(20,'JOBARC','HFKSFIAB',occvrt2*iintfp,
     &             tfib)
      end if
      return
      end   
