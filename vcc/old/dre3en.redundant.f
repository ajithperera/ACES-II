
c This routine drives the third-order MBPT energy calculation.

      subroutine dre3en(iCore,iCoreDim,iUHF,iDoPPL)
      implicit none

      integer iCore(*), iCoreDim, iUHF, iDoPPL

      INTEGER        iFlags(100)
      COMMON /FLAGS/ iFlags

      LOGICAL         MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC
      COMMON /SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC

      LOGICAL       ROHF4,ITRFLG
      COMMON /ROHF/ ROHF4,ITRFLG

c ----------------------------------------------------------------------

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

      if (iFlags(2).gt.9.and.SING1.and.(.not.QCISD).or.
     &    (ROHF4.and.ITRFLG)
     &   ) then
        call t12int2(iCore,iCoreDim,iUHF)
      end if

      if (iUHF.eq.0) then
         call drrng(iCore,iCoreDim,3,iUHF)
      else
         call drrng(iCore,iCoreDim,1,iUHF)
         call drrng(iCore,iCoreDim,2,iUHF)
         call drrng(iCore,iCoreDim,3,iUHF)
      end if

      return
c     end subroutine dre3en
      end

