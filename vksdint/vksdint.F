       subroutine  vksdint

       implicit none

#include <aces.h>
#include <flags.com>
#include <flags2.com>
#include <iuhf.com>
#include <mol.com>
#include <kscore.com>
#include <symm2_ks.com>
#include <timeinfo.com>
#include <int.com>
#include <intgrtflags.com>

      integer intnumradpts,max_angpts,islateron
      integer maxmem, i0tmp
      logical kshf

      call callstack_init('INTGRT')

      int_ks=.true.
     
      call timer(0)
c      call aces_init_rte
c      call aces_ja_init

      call initintgrt(kshf,intnumradpts,islateron)
      call getrec(1,'JOBARC','IFLAGS  ',100,iflags)
      call getrec(1,'JOBARC','IFLAGS2 ',500,iflags2)
      iuhf=0
      if (iflags(11).ge.1) iuhf=1
 
      call getrec(-1,'JOBARC','COMPNIRR',1,nirrep)
      call sb_com_symm2_ks

      maxmem=iflags(36)
      maxicore=500000

c      do memknown=0,1
      do memknown=-1,-1
         i0tmp = 1
         if (memknown.ne.0) then
           kscore(1) = 0
            do while ((kscore(1).eq.0).and.(maxmem.gt.0))
               call aces_malloc(maxmem,kscore,i0tmp)
               if (kscore(1).eq.0) maxmem = maxmem - 1024*1024
            end do
         end if
         call ks_mem_init(kscore(i0tmp),maxmem)
         call sb_com_mol

         call anlytc
         call setupgrdAG(max_angpts)

         idns=0
 900     continue

         call primvecint
         call numintAG(.true.)


         call sb_mem_fin


c     end do memknown=0,1
      end do

c      call aces_ja_fin
      call timer(1)
      write(*,9020) timetot

 9020 format(t3,'@CHECKOUT-I, Total execution time : ',f12.4,
     &       ' seconds.')


      call callstack_term
      Return
      end
 


