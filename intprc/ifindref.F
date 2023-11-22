      integer function ifindref(nbas,erank,efrom,eto,espin)
c
cjp this function finds to which reference configuration a given
cjp excitation gives rise to, is used to identify Heff offdiagonal
cjp matrix elements in focklist (singles), dijab (doubles)
cjp and is written in a general way to be used later for triples too.
c
      implicit none
#include "bwcc.com"
c
      integer nbas
      integer erank,efrom,eto,espin
      dimension efrom(*),eto(*),espin(*)
c
      integer i,j,k,result
      integer xocc
      dimension xocc(maxorb,2)
c
cjp check legitimacy of that excitation and prepare the target occupation pattern
      do i=1,2
      do j=1,nbas
        xocc(j,i)=iocc(j,iref,i)
      enddo
      enddo
      do i=1,erank
      xocc(efrom(i),espin(i)) = xocc(efrom(i),espin(i)) - 1
      xocc(eto(i),espin(i)) = xocc(eto(i),espin(i)) + 1
      enddo
      do i=1,2
      do j=1,nbas
         if(xocc(j,i) .lt. 0 .or. xocc(j,i) .gt. 1) then
              write(6,*) 'internal error in ifindref detected',iref,j,i
              stop 'internal error in ifindref detected'
         endif
      enddo
      enddo
c
cjp now search for xocc pattern between all references
      result=0
      do k=1,nref
      if(k.ne.iref) then
          do i=1,2
          do j=1,nbas
              if(xocc(j,i).ne.iocc(j,k,i)) goto 1
          enddo
          enddo
          result=k
          goto 2
1         continue
      endif
      enddo
2     continue
c
      if(result.eq.0) write(6,*) '@IFINDREF-I there are internal excita
     +tions not corresponding to any reference configuration. You are u
     +sing an incomplete model space.'
c
      ifindref=result
      return
      end
