#include "flags.h"
      Subroutine Tdcc_form_lists(Work,Memleft,Iuhf)

      Implicit Double Precision (A-H, O-Z)
      Integer AAAA_LENGTH_IJAB,BBBB_LENGTH_IJAB,AABB_LENGTH_IJAB
      Integer PHA_Length,PHB_Length,HHA_length,HHB_Length
      Integer PPA_Length,PPB_Length

      Dimension Work(Memleft)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "parcc.com"
#include "flags.com"

       Imode = 0

C There are three type of amplitudes of singles and doubles dimension.
C Let us call them and create pointers as follows
C 
C 390,392,394 are singles lists of mu_0, mu_t, mu_d (dot) and mu_td
C amplitudes (i,a).

      MAX_PHA_Length = 0
      MAX_PHB_Length = 0
      Do Irrep = 1, Nirrep 
         MAX_PHA_Length = Max(MAX_PHA_Length,Irpdpd(Irrep,9))
         MAX_PHB_Length = Max(MAX_PHB_Length,Irpdpd(Irrep,10))
      Enddo 
      
      CALL Updmoi(1,MAX_PHA_length,1,390,0,0)
      CALL Updmoi(1,MAX_PHA_length,1,392,0,0)
      CALL Updmoi(1,MAX_PHA_length,1,394,0,0)
      CALL Updmoi(1,MAX_PHA_length,1,396,0,0)
      Call Aces_list_memset(1,390,0)
      Call Aces_list_memset(1,392,0)
      Call Aces_list_memset(1,394,0)
      Call Aces_list_memset(1,396,0)

      If (Iuhf .NE. 0) Then
         CALL Updmoi(1,MAX_PHB_length,2,390,0,0)
         CALL Updmoi(1,MAX_PHB_length,2,392,0,0)
         CALL Updmoi(1,MAX_PHB_length,2,394,0,0)
         CALL Updmoi(1,MAX_PHB_length,2,396,0,0)
         Call Aces_list_memset(2,390,0)
         Call Aces_list_memset(2,392,0)
         Call Aces_list_memset(2,394,0)
         Call Aces_list_memset(2,396,0)
      Endif

C The doubles list for mubar,mubar_tilde,mubar_dot,mubar_tilde_dot.
C At this instance left-hand side is arbitrary. So irrep_arb can
C be any value of 1-8.

      Irrep_arb = 1
      If (Iuhf .EQ. 0) Then
         Call Inipck2(Irrep_arb,13,14,316,Imode,0,1)
         Call Inipck2(Irrep_arb,13,14,326,Imode,0,1)
         Call Inipck2(Irrep_arb,13,14,336,Imode,0,1)
         Call Inipck2(Irrep_arb,13,14,346,Imode,0,1)
      Else
C Mubar
         Call Inipck2(Irrep_arb,1,3,314,Imode,0,1)
         Call Inipck2(Irrep_arb,2,4,315,Imode,0,1)
         Call Inipck2(Irrep_arb,13,14,316,Imode,0,1)
C Mubar_tilde
         Call Inipck2(Irrep_arb,1,3,324,Imode,0,1)
         Call Inipck2(Irrep_arb,2,4,325,Imode,0,1)
         Call Inipck2(Irrep_arb,13,14,326,Imode,0,1)
C Mubar_dot
         Call Inipck2(Irrep_arb,1,3,334,Imode,0,1)
         Call Inipck2(Irrep_arb,2,4,335,Imode,0,1)
         Call Inipck2(Irrep_arb,13,14,336,Imode,0,1)

C Mubar_tilde_dot
         Call Inipck2(Irrep_arb,1,3,344,Imode,0,1)
         Call Inipck2(Irrep_arb,2,4,345,Imode,0,1)
         Call Inipck2(Irrep_arb,13,14,346,Imode,0,1)
      Endif

C Temporary storage of (ai,bi) quantities 
    
      If (Iuhf .EQ. 0) Then
         Call Inipck2(Irrep_arb,9,10,302,Imode,0,1)
      Else
         Call Inipck2(Irrep_arb,9,9,300,Imode,0,1)
         Call Inipck2(Irrep_arb,10,10,301,Imode,0,1)
         Call Inipck2(Irrep_arb,9,10,302,Imode,0,1)
      Endif

C Lists for MO dipole moment integrals (OO,OV and VV blocks).

      MAX_HHA_Length = 0
      MAX_HHB_Length = 0
      MAX_PPA_Length = 0
      MAX_PPB_Length = 0

      Do irrep = 1, Nirrep 
         MAX_HHA_Length = Max(MAX_HHA_Length,Irpdpd(Irrep,21))
         MAX_HHB_Length = Max(MAX_HHB_Length,Irpdpd(Irrep,22))
         MAX_PPA_Length = Max(MAX_PPA_Length,Irpdpd(Irrep,19))
         MAX_PPB_Length = Max(MAX_PPB_Length,Irpdpd(Irrep,20))
      Enddo 

      CALL Updmoi(1,MAX_HHA_length,1,380,0,0)
      CALL Updmoi(1,MAX_PPA_length,1,382,0,0)
      CALL Updmoi(1,MAX_PHA_length,1,384,0,0)
      Call Aces_list_memset(1,380,0)
      Call Aces_list_memset(1,382,0)
      Call Aces_list_memset(1,384,0)

      If (Iuhf .NE. 0) Then
         CALL Updmoi(1,MAX_HHB_length,2,380,0,0)
         CALL Updmoi(1,MAX_PPB_length,2,382,0,0)
         CALL Updmoi(1,MAX_PHB_length,2,384,0,0)
         Call Aces_list_memset(2,380,0)
         Call Aces_list_memset(2,382,0)
         Call Aces_list_memset(2,384,0)
      Endif

      Return
      End
