#include "flags.h"
      Subroutine Tdee_form_lists(Work,Memleft,Iuhf)

      Implicit Double Precision (A-H, O-Z)
      Integer AAAA_LENGTH_IJAB,BBBB_LENGTH_IJAB,AABB_LENGTH_IJAB
      Integer PHA_Length,PHB_Length,HHA_length,HHB_Length
      Integer PPA_Length,PPB_Length
      LOGICAL MBPT2,CC,SS,SD,DS,DD,CCD,RPA,DRPA,CCD_DUMMY,MBPT2_DUMMY
      LOGICAL LCCD,LCCSD,CC2,CCSD

      COMMON /METH/MBPT2_DUMMY,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD_DUMMY,
     +             QCISD,CCSD,UCC
      COMMON /REFTYPE/ MBPT2,CC,CCD,RPA,DRPA,LCCD,LCCSD,CC2

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
      MAX_PPA_Length = 0
      MAX_PPB_Length = 0
      MAX_HHA_Length = 0
      MAX_HHB_Length = 0

      Do Irrep = 1, Nirrep 
         MAX_PHA_Length = Max(MAX_PHA_Length,Irpdpd(Irrep,9))
         MAX_PHB_Length = Max(MAX_PHB_Length,Irpdpd(Irrep,10))
         MAX_PPA_Length = Max(MAX_PPA_Length,Irpdpd(Irrep,19))
         MAX_PPB_Length = Max(MAX_PPB_Length,Irpdpd(Irrep,20))
         MAX_HHA_Length = Max(MAX_HHA_Length,Irpdpd(Irrep,21))
         MAX_HHB_Length = Max(MAX_HHB_Length,Irpdpd(Irrep,22))
      Enddo 
      
      CALL Updmoi(1,MAX_PHA_length,1,390,0,0)
      CALL Updmoi(1,MAX_PHA_length,3,390,0,0)
      CALL Updmoi(1,MAX_PHA_length,1,392,0,0)
      CALL Updmoi(1,MAX_PHA_length,3,392,0,0)
      CALL Updmoi(1,MAX_PHA_length,1,394,0,0)
      CALL Updmoi(1,MAX_PHA_length,1,396,0,0)
      CALL Updmoi(1,MAX_HHA_length,1,381,0,0)
      CALL Updmoi(1,MAX_PPA_length,1,382,0,0)

      Call Aces_list_memset(1,390,0)
      Call Aces_list_memset(3,390,0)
      Call Aces_list_memset(1,392,0)
      Call Aces_list_memset(3,392,0)
      Call Aces_list_memset(1,394,0)
      Call Aces_list_memset(1,396,0)
      Call Aces_list_memset(1,381,0)
      Call Aces_list_memset(1,382,0)

      If (Iuhf .NE. 0) Then
         CALL Updmoi(1,MAX_PHB_length,2,390,0,0)
         CALL Updmoi(1,MAX_PHB_length,4,390,0,0)
         CALL Updmoi(1,MAX_PHB_length,2,392,0,0)
         CALL Updmoi(1,MAX_PHB_length,4,392,0,0)
         CALL Updmoi(1,MAX_PHB_length,2,394,0,0)
         CALL Updmoi(1,MAX_PHB_length,2,396,0,0)
         CALL Updmoi(1,MAX_HHB_length,2,381,0,0)
         CALL Updmoi(1,MAX_PPB_length,2,382,0,0)

         Call Aces_list_memset(2,390,0)
         Call Aces_list_memset(4,390,0)
         Call Aces_list_memset(2,392,0)
         Call Aces_list_memset(4,392,0)
         Call Aces_list_memset(2,394,0)
         Call Aces_list_memset(2,396,0)
         Call Aces_list_memset(2,381,0)
         Call Aces_list_memset(2,382,0)
      Endif

C The list 190(1 and 2) L1 lists are not created for CCD in CC code.
C Lets create them here and set them to zero as they should be.

      If (CCD_DUMMY) Then 
         CALL Updmoi(1,MAX_PHA_length,1,190,0,0)
         Call Aces_list_memset(1,190,0)
         If (Iuhf .NE. 0) Then
            CALL Updmoi(1,MAX_PHB_length,2,190,0,0)
            Call Aces_list_memset(2,190,0)
         Endif 
      Endif 

C The doubles lists for mu_t, mu_d and mu_td (ij,ab). (there are no
C mu_0 doiuble lists). 
C At this instance left-hand side is arbitrary. So irrep_arb can
C be any value of 1-8.

      Irrep_arb = 1
      If (Iuhf .EQ. 0) Then
         Call Inipck2(Irrep_arb,13,14,316,Imode,0,1)
         Call Inipck2(Irrep_arb,13,14,319,Imode,0,1)
         Call Inipck2(Irrep_arb,13,14,326,Imode,0,1)
         Call Inipck2(Irrep_arb,13,14,329,Imode,0,1)
         Call Inipck2(Irrep_arb,13,14,336,Imode,0,1)
         Call Inipck2(Irrep_arb,13,14,346,Imode,0,1)
      Else
C Mu
         Call Inipck2(Irrep_arb,1,3,314,Imode,0,1)
         Call Inipck2(Irrep_arb,2,4,315,Imode,0,1)
         Call Inipck2(Irrep_arb,13,14,316,Imode,0,1)

         Call Inipck2(Irrep_arb,1,3,317,Imode,0,1)
         Call Inipck2(Irrep_arb,2,4,318,Imode,0,1)
         Call Inipck2(Irrep_arb,13,14,319,Imode,0,1)

C Mu_tilde
         Call Inipck2(Irrep_arb,1,3,324,Imode,0,1)
         Call Inipck2(Irrep_arb,2,4,325,Imode,0,1)
         Call Inipck2(Irrep_arb,13,14,326,Imode,0,1)

         Call Inipck2(Irrep_arb,1,3,327,Imode,0,1)
         Call Inipck2(Irrep_arb,2,4,328,Imode,0,1)
         Call Inipck2(Irrep_arb,13,14,329,Imode,0,1)
C Mu_dot
         Call Inipck2(Irrep_arb,1,3,334,Imode,0,1)
         Call Inipck2(Irrep_arb,2,4,335,Imode,0,1)
         Call Inipck2(Irrep_arb,13,14,336,Imode,0,1)

C Mu_tilde_dot
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

C AO ladder lists and other auxilary lists.
    
      If (Iuhf .EQ. 0) Then
         If (Iflags(93) .Eq. 2) Then
            Call Inipck2(Irrep_arb,15,14,377,Imode,0,1)
            Call Inipck2(Irrep_arb,15,14,387,Imode,0,1)
         Endif 
      Else
         If (Iflags(93) .Eq. 2) Then
            Call Inipck2(Irrep_arb,15,3,375,Imode,0,1)
            Call Inipck2(Irrep_arb,15,4,376,Imode,0,1)
            Call Inipck2(Irrep_arb,15,14,377,Imode,0,1)

            Call Inipck2(Irrep_arb,15,3,385,Imode,0,1)
            Call Inipck2(Irrep_arb,15,4,386,Imode,0,1)
            Call Inipck2(Irrep_arb,15,14,387,Imode,0,1) 
          Endif 
      Endif 

C Other miscellaneous lists 

      Call Inipck2(Irrep_arb,9,9,304,Imode,0,1)
      Call Inipck2(Irrep_arb,9,10,307,Imode,0,1)
      Call Inipck2(Irrep_arb,11,12,309,Imode,0,1)
      Call Inipck2(Irrep_arb,9,9,350,Imode,0,1)
      Call Inipck2(Irrep_arb,10,10,351,Imode,0,1)
      Call Inipck2(Irrep_arb,9,10,352,Imode,0,1)
      Call Inipck2(Irrep_arb,11,12,353,Imode,0,1)
      Call Inipck2(Irrep_arb,9,9,360,Imode,0,1)
      Call Inipck2(Irrep_arb,10,10,361,Imode,0,1)
      Call Inipck2(Irrep_arb,9,10,362,Imode,0,1)
      Call Inipck2(Irrep_arb,11,12,363,Imode,0,1)

      If (Iuhf .NE. 0) Then
         Call Inipck2(Irrep_arb,10,10,305,Imode,0,1)
         Call Inipck2(Irrep_arb,10,9,306,Imode,0,1)
         Call Inipck2(Irrep_arb,12,11,308,Imode,0,1)
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

      CALL Updmoi(1,MAX_HHA_length,1,370,0,0)
      CALL Updmoi(1,MAX_PPA_length,1,371,0,0)
      CALL Updmoi(1,MAX_PHA_length,1,372,0,0)
      Call Aces_list_memset(1,370,0)
      Call Aces_list_memset(1,371,0)
      Call Aces_list_memset(1,372,0)

      If (Iuhf .NE. 0) Then
         CALL Updmoi(1,MAX_HHB_length,2,370,0,0)
         CALL Updmoi(1,MAX_PPB_length,2,371,0,0)
         CALL Updmoi(1,MAX_PHB_length,2,372,0,0)
         Call Aces_list_memset(2,370,0)
         Call Aces_list_memset(2,371,0)
         Call Aces_list_memset(2,372,0)
      Endif

      Return
      End
