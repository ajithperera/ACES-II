      subroutine dre3en(iCore,iCoreDim,iUHF,iDoPPL)
C
#include "flags.h" /* for h_IFLAGS_* and h_IFLAGS2_* */

      implicit none

      integer iCore(*), iCoreDim, iUHF, iDoPPL, iflags,
     &        iflags2
C
      LOGICAL         MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC,NO_REDUNDANT
      LOGICAL       ROHF4,ITRFLG
      COMMON /SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC
      COMMON /FLAGS/ iFlags(100)
      COMMON /FLAGS2/ iFlags2(500)
      COMMON /ROHF/ ROHF4,ITRFLG

c ----------------------------------------------------------------------
C
      if (iDoPPL.eq.0) then
         call drlad(iCore,iCoreDim,iUHF,1)
         if (iFlags(93).eq.2) then
            call draolad(iCore,iCoreDim,iUHF,.false.,1,1,43,60,243,260)
         else
            call drlad(iCore,iCoreDim,iUHF,6)
         end if
      else
         call drlad(iCore,iCoreDim,iUHF,1)
      end if

#ifdef _DEBUG_LVL0
      Write(6,*) "After mnij+abef"
      call check_t2(icore,Icoredim,Iuhf)
#endif 
      if (iFlags(2).gt.9.and.SING1.and.(.not.QCISD).or.
     &    (ROHF4.and.ITRFLG)
     &   ) then
           call t12int2(iCore,iCoreDim,iUHF)
      end if
 
#ifdef _DEBUG_LVL0
      Write(6,*) "After t1int2"
      call check_t2(icore,Icoredim,Iuhf)
#endif 
C
         if (iUHF.eq.0) then
            call drrng(iCore,iCoreDim,3,iUHF)
         else
            call drrng(iCore,iCoreDim,1,iUHF)
            call drrng(iCore,iCoreDim,2,iUHF)
            call drrng(iCore,iCoreDim,3,iUHF)
         end if
#ifdef _DEBUG_LVLM
      Write(6,*) "After       mbej"
      call check_t2(icore,Icoredim,Iuhf)
#endif 
C
      return
c     end subroutine dre3en
      end

