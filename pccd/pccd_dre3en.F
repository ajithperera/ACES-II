      subroutine pccd_dre3en(iCore,iCoreDim,iUHF,iDoPPL)
C
#include "flags.h" /* for h_IFLAGS_* and h_IFLAGS2_* */

      implicit none

      integer iCore(*), iCoreDim, iUHF, iDoPPL, iflags,
     &        iflags2
      COMMON /FLAGS/ iFlags(100)
      COMMON /FLAGS2/ iFlags2(500)
C
      if (iDoPPL.eq.0) then
         call Pccd_drlad(iCore,iCoreDim,iUHF,1)
         if (iFlags(93).eq.2) then
            Call Pccd_draolad(iCore,iCoreDim,iUHF,.false.,1,1,
     &                        43,60,243,260)
         else
            Call Pccd_drlad(iCore,iCoreDim,iUHF,6)
         end if
      else
         Call Pccd_drlad(iCore,iCoreDim,iUHF,1)
      end if

#ifdef _DEBUG_LVLM
      call pccd_check_t2(icore,Icoredim,Iuhf)
#endif 
      Call Pccd_drrng(iCore,iCoreDim,3,iUHF)

#ifdef _DEBUG_LVLM
      call pccd_check_t2(icore,Icoredim,Iuhf)
#endif 
C
      return
c     end subroutine dre3en
      end

