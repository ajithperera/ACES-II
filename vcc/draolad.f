
c This routine drives the AO-based particle-particle ladder contraction(s).

      subroutine draolad(iCore,iCoreDim,iUHF,bLambda,irp_x,iFlag_list,
     &                   list_mo,list_mo_inc,list_ao,list_ao_inc)
      implicit none

      integer iCore(*), iCoreDim, iUHF, irp_x, iFlag_list
      integer list_mo, list_mo_inc, list_ao, list_ao_inc
      logical bLambda

      integer iMethod
      logical bTau, bSing

      LOGICAL         MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC
      COMMON /SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC

      INTEGER        iFlags(100)
      COMMON /FLAGS/ iFlags

      LOGICAL       ROHF4,ITRFLG
      COMMON /ROHF/ ROHF4,ITRFLG

c ----------------------------------------------------------------------

      iMethod = iFlags(2)
      bSing = iMethod.ne. 2.and.
     &        iMethod.ne. 5.and.
     &        iMethod.ne. 8.and.
     &        iMethod.ne.20.and.
     &        iMethod.ne.24
      bTau = iMethod.gt.9.and.SING1.and.(.not.QCISD).or.
     &       (ROHF4.and..not.ITRFLG)

c ----------------------------------------------------------------------

      call aoladlst(iUHF,iFlag_list,irp_x)

      call t2toao(iCore,iCoreDim,iUHF,bTau,list_mo,list_ao,irp_x)

      if (iFlags(95).eq.1) then
         call aolad2(iCore,iCoreDim,iUHF,bTau,irp_x,list_ao,list_ao_inc)
      else
         call aolad3(iCore,iCoreDim,iUHF,bTau,irp_x,list_ao,list_ao_inc)
      end if

      call z2tomo(iCore,iCoreDim,iUHF,bTau,irp_x,
     &            list_mo,list_mo_inc,list_ao_inc,.false.)

      if (bSing) call ziaao(iCore,iCoreDim,iUHF,irp_x,90)

      return
c     end subroutine draolad
      end

