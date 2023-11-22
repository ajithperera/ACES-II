
      subroutine setderiv(listl1,listl1off,listl2,listl2rs,listr1,
     &                    listr1off,listr2,listr2rs,listgrl,listgtl,
     &                    listgrlof,listgtlof,listtmp,listtmpof,listt1,
     &                    listt1off,listz1,listz1off,listz2,uhf,nonhf)
      implicit none
C Common blocks
#include "flags.h"
#include "flags.com"
      logical ss,sd,ds,dd
      common/drvhbar/ss,sd,ds,dd
      integer irrepx
      common/statsym/irrepx
      logical gabcd
      common/abcd/gabcd
      integer maxexp,nreduce,ntol,nsizec
      common/extrap/maxexp,nreduce,ntol,nsizec
      integer irpdpd(8,22),isytyp(2,500),id(18)
      common/sympop/irpdpd,isytyp,id
      logical mbpt2,cc
      common/reftype/mbpt2,cc
C Input variables
      integer uhf
C Output variables
      logical nonhf
      integer listl1,listl1off,listl2,listl2rs,listr1,listr1off,listr2,
     &        listr2rs,listgrl,listgtl,listgrlof,listgtlof,listtmp,
     &        listtmpof,listt1,listt1off,listz1,listz1off,listz2
C External functions
      integer idsymsz
      external idsymsz
C Local variables
      integer spin,len

      ss=.true.
      sd=.true.
      ds=.true.
      dd=.true.
      mbpt2=.false.
      cc=.true.
      gabcd=iflags(h_IFLAGS_gamma_abcd).eq.1
      irrepx=1
      nonhf = (iflags(h_IFLAGS_nonhf).ne.0).or.
     &   (iflags(h_IFLAGS_reference).eq.2).or.
     &   ((iflags(h_IFLAGS_reference).eq.1).and.
     &   ((iflags(32).ne.0).or.(iflags(33).ne.0)))
C Lists
      listr1=290
      listl1=190
      listr1off=0
      listl1off=0
      listl2=144
      listr2=244
      listl2rs=134
      listr2rs=234
      listgrl=400
      listgtl=0
      listgrlof=2
      listgtlof=2
      listtmp=290
      listtmpof=2
      listt1=90
      listt1off=0
      listz1=390
      listz1off=0
      listz2=344
      maxexp = iflags(h_IOPPAR_estate_maxc)
      nreduce = min(12,maxexp-3)
      ntol = iflags(h_IFLAGS_estate_tol)
      len = 0
      do spin = 1,uhf+1
        len = len+irpdpd(1,8+spin)
      end do
      len = len + idsymsz(1,isytyp(1,46),isytyp(2,46))
      if (uhf .ne. 0) then
        len = len + idsymsz(1,isytyp(1,44),isytyp(2,44))
        len = len + idsymsz(1,isytyp(1,45),isytyp(2,45))
      endif
      nsizec = len

      return
      end
