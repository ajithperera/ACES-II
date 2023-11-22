
c This routine drives the AO-based particle-particle ladder contraction(s).

      subroutine pccd_draolad(iCore,iCoreDim,iUHF,bLambda,irp_x,
     &                        iFlag_list,list_mo,list_mo_inc,
     &                        list_ao,list_ao_inc)
      implicit none

      integer iCore(*), iCoreDim, iUHF, irp_x, iFlag_list
      integer list_mo, list_mo_inc, list_ao, list_ao_inc
      logical bLambda
      logical bTau, bSing
      INTEGER        iFlags(100)
      COMMON /FLAGS/ iFlags

c ----------------------------------------------------------------------

      bTau = .False.
      bSing = .False.
c ----------------------------------------------------------------------

      call pccd_aoladlst(iUHF,iFlag_list,irp_x)

      call pccd_t2toao(iCore,iCoreDim,iUHF,bTau,list_mo,list_ao,irp_x)

      if (iFlags(95).eq.1) then
         call pccd_aolad2(iCore,iCoreDim,iUHF,bTau,irp_x,list_ao,
     +                    list_ao_inc)
      else
         call pccd_aolad3(iCore,iCoreDim,iUHF,bTau,irp_x,list_ao,
     +                    list_ao_inc)
      end if

      call pccd_z2tomo(iCore,iCoreDim,iUHF,bTau,irp_x,
     +                 list_mo,list_mo_inc,list_ao_inc,.false.)

      return
c     end subroutine draolad
      end

