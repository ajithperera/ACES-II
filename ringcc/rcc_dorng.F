      subroutine RCC_DORNG(Work,Maxcor,iUHF)
C
#include "flags.h" /* for h_IFLAGS_* and h_IFLAGS2_* */

      implicit none

      integer WORK(Maxcor),Maxcor,iUHF
      Logical MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD
C
      If (iUHF.eq.0) then
         If (DRCCD) THEN
            Call Rcc_rng(Work,Maxcor,3,iUHF)
         Elseif (RCCD) THEN
            Call Rcc_rng_4rhf_rccd(Work,Maxcor,3,iUHF)
         Endif 
      Else
         Call Rcc_rng(Work,Maxcor,1,iUHF)
         Call Rcc_rng(Work,Maxcor,2,iUHF)
         Call Rcc_rng(Work,Maxcor,3,iUHF)
      Endif
C
      Return
      End

