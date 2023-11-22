
      subroutine resetfno(eval,nbas,file)
      implicit none
#include "aces2.h"
C     Common blocks
      integer iintln,iifltln,iintfp 
      common/machsp/iintln,iifltln,iintfp 
C     Input Variables
      integer nbas
      character*2 file
C     Pre-allocated Local Variables
      double precision eval(nbas)      
C     Local variables
      integer nirrep,nocco(2),vrtdrop(8,2),vrt(8,2),uhf,virtd,virtk,
     &   virt,irrep,ii,spin,ieval,dropocc
      character*8 recname(2)
      data recname /'SCFEVLA0','SCFEVLB0'/
C - - - - - - -- -- - - - - - - - - - - - - - - - - - - - - - - -
      uhf = 0
      if (iflags(11).ne.0) uhf = 1
      call aces_ja_init
      call getrec(20,'JOBARC','NOCCORB',2,nocco)
      call getrec(20,'JOBARC','COMPNIRR',1,nirrep)
      call getrec(20,'JOBARC','SYMPOPVA',nirrep,vrt(1,1))
      call getrec(20,'JOBARC','SYMPOPVB',nirrep,vrt(1,2))
      call getrec(20,'JOBARC','FNODROP ',16,vrtdrop)
      call getrec(20,'JOBARC','FNOFREEZ',1,dropocc)
      do spin=1,uhf+1
         call getrec(20,'JOBARC',recname(spin),nbas*iintfp,eval)
         if (file.eq.'DM') then
            ieval=nocco(spin)+dropocc
         else
            ieval=nocco(spin)
         endif
         do irrep=1,nirrep
            virt=vrt(irrep,spin)
            virtd=vrtdrop(irrep,spin)
            if (file.eq.'DM') then
               virtk=virt
            else
               virtk=virt-virtd
            endif
            ieval=ieval+virtk
            do ii=1,virtd
               eval(ieval+ii)=eval(ieval+ii)-999.0d0
            end do
            ieval=ieval+virtd
         end do
         call putrec(20,'JOBARC',recname(spin),nbas*iintfp,eval)
      end do
      call aces_ja_fin
      end
