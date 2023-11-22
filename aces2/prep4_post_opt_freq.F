#include "flags.h"
      Subroutine Prep4_post_opt_freq
      
      Implicit none
c
#include "aces2.h" 
#include "linelen.par"

      Integer Iflag_store1(100), Length, I_havegeom, iOPTARC,
     &        Iflag_store2(500), Itmp, I_lastgeom
      Character*(linelen) szOPTARC, szTmp
      Logical OPTARC_presnt 
c
      Call A2getrec(10, 'JOBARC', 'IFLAGS  ', 100, Iflag_store1)
      Call A2getrec(10, 'JOBARC', 'IFLAGS2 ', 500, Iflag_store2)
      Call A2getrec(10,  'JOBARC', 'HAVEGEOM', 1, I_havegeom)
c
      Call gfname('OPTARC',szOPTARC,iOPTARC)
      Inquire(file=szOPTARC(1:iOPTARC),exist=OPTARC_presnt)
c
#ifdef _DEBUG
      Print*, "I_havegeom,OPTARC's presence,Iflags(h_IFLAGS_vib)"
      Print*,  I_havegeom, OPTARC_presnt, Iflags(h_IFLAGS_vib)
#endif

      If (I_havegeom .GT. 0 .AND. OPTARC_presnt .AND. 
     &    Iflags(h_IFLAGS_vib) .GT. 0)
     &    Call A2putrec(10, 'JOBARC', 'VIB_POPT', 1, 1)
c
      Iflag_store2(h_IFLAGS2_geom_opt) = 0
      Call A2getrec(10, 'JOBARC', 'NFDIRREP', 1, ITmp)
      If (Itmp .NE. 0) Itmp = 0
      Call A2putrec(10, 'JOBARC', 'NFDIRREP', 1, Itmp)
      Call A2getrec(10, 'JOBARC', 'FNDFDONE', 1, iTmp)
      If (Itmp .NE. 0) Itmp = 0
      Call A2putrec(10, 'JOBARC', 'FNDFDONE', 1, iTmp)
c
#ifdef _DEBUG
      Print*, "@Pre4_Freq: vib, geom_opt, deriv_lev, grad_calc,
     &         fd_stepsize flags"
      Print*, Iflags(h_IFLAGS_vib), Iflags2(h_IFLAGS2_geom_opt),
     & Iflags(h_IFLAGS_deriv_lev), Iflags2(h_IFLAGS2_grad_calc),
     & Iflags(h_IFLAGS_fd_stepsize)
#endif 
c 
c   o Prepare the derivative level for vibrational frequency.
c
      If (Iflags(h_IFLAGS_vib) .GE. 3) 
     &    Iflag_store1(h_IFLAGS_deriv_lev) = 1
c
c   o VIB=EXACT need especial care
c
      If (Iflags(h_IFLAGS_calc) .EQ. 1 .AND. Iflags(h_IFLAGS_vib) 
     &   .EQ. 1) Iflag_store1(h_IFLAGS_deriv_lev) = 2
c
      If ((Iflags(h_IFLAGS_vib) .EQ. 1)        .AND. 
     &    (Iflags2(h_IFLAGS2_geom_opt) .NE. 0) .AND. 
     &     I_havegeom .GT. 0) 
     &     Iflag_store1(h_IFLAGS_deriv_lev) = 2
c
c   o Set the step size for the numerical vibrational frequency calculation
c
      If (Iflags(h_IFLAGS_vib) .GE. 3) Then
         If (Iflags2(h_IFLAGS2_grad_calc).eq.2) then
c          o numerical gradients
                Iflag_store1(h_IFLAGS_fd_stepsize) = 200
         Else
                Iflag_store1(h_IFLAGS_fd_stepsize) = 50
         Endif
      Endif
#ifdef _DEBUG_LVL0

       Print*, "The fd_stepsize is geting updated (new value)"
       Print*, Iflag_store1(h_IFLAGS_fd_stepsize)
#endif
c
      Call A2putrec(10, 'JOBARC', 'IFLAGS  ', 100, Iflag_store1)
      Call A2putrec(10, 'JOBARC', 'IFLAGS2 ', 500, Iflag_store2)
c
c    o Truncate the JOBARC from JODAOUT for so that lower symmetry
c      points can pproperly run
c
      Call Aces_ja_init
      Call Aces_ja_truncate('JODAOUT ',1)
      Call Aces_ja_fin

      If ((Iflags(h_IFLAGS_vib) .GE. 3) .AND. 
     &     Iflags2(h_IFLAGS2_grad_calc) .eq. 2) Then
           Call A2putrec(1, 'JOBARC', 'FNDFDONE', 1, 0)
           Call A2putrec(1, 'JOBARC', 'POSTOPFD', 1, 1)
      Endif
c
      Call Runit("xjoda") 
c
      If ((Iflags(h_IFLAGS_vib) .GE. 3) .AND. 
     &     Iflags2(h_IFLAGS2_grad_calc) .eq. 2) Then
           Call A2putrec(1, 'JOBARC', 'POSTOPFD', 1, 0)
      Endif
c
c    o Restore the records that got deleted for future use.
c
      Call A2putrec(10, 'JOBARC', 'LASTGEOM', 1, 0)
      Call A2putrec(10, 'JOBARC', 'HAVEGEOM', 1, I_havegeom)
c
      Call Runit("rm -f OPTARC")
c
#ifdef _DEBUG
      Call A2getrec(10, 'JOBARC', 'IFLAGS  ', 100, Iflag_store1)
      Print*, "Test the jobarc, derivative level:",
     &         iflag_store1(3)
#endif
      Return
      End

