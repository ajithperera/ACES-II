#include "flags.h"
      Subroutine Tdee_reset_lists(Work,Memleft,Iuhf,Irrepx)

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

C There are three type of amplitudes of singles and doubles dimension.
C Let us call them and resized the pointers as follows
C 
C 390,392,394 are singles lists of mu_0, mu_t, mu_d (dot) and mu_td
C amplitudes (i,a).

      PHA_Length = Irpdpd(Irrepx,9)
      PHB_Length = Irpdpd(Irrepx,10)
      PPA_Length = Irpdpd(Irrepx,19)
      PPB_Length = Irpdpd(Irrepx,20)
      HHA_Length = Irpdpd(Irrepx,21)
      HHB_Length = Irpdpd(Irrepx,22)
      
      Call Aces_list_resize(1,390,PHA_Length)
      Call Aces_list_resize(3,390,PHA_Length)
      Call Aces_list_resize(1,392,PHA_Length)
      Call Aces_list_resize(3,392,PHA_Length)
      Call Aces_list_resize(1,394,PHA_Length)
      Call Aces_list_resize(1,396,PHA_Length)
      Call Aces_list_resize(1,381,HHA_Length)
      Call Aces_list_resize(1,382,PPA_Length)
      If (CCD_DUMMY) Call Aces_list_resize(1,190,PHA_Length)

      If (Iuhf .NE. 0) Then
         Call Aces_list_resize(2,390,PHB_Length)
         Call Aces_list_resize(4,390,PHB_Length)
         Call Aces_list_resize(2,392,PHB_Length)
         Call Aces_list_resize(4,392,PHB_Length)
         Call Aces_list_resize(2,394,PHB_Length)
         Call Aces_list_resize(2,396,PHB_Length)
         Call Aces_list_resize(2,381,HHB_Length)
         Call Aces_list_resize(2,382,PPB_Length)
         If (CCD_DUMMY) Call Aces_list_resize(2,190,PHA_Length)
      Endif
 

C The doubles lists for mu_t, mu_d and mu_td (ij,ab). (there are no
C mu_0 doiuble lists). Originally these were initialized to irrep1,
C but can accomodate any irrep.  

      If (Iuhf .EQ. 0) Then
         Call Newtyp2(Irrepx,316,13,14,.True.)
         Call Newtyp2(Irrepx,319,13,14,.True.)
         Call Newtyp2(Irrepx,326,13,14,.True.)
         Call Newtyp2(Irrepx,329,13,14,.True.)
         Call Newtyp2(Irrepx,336,13,14,.True.)
         Call Newtyp2(Irrepx,346,13,14,.True.)
      Else
C Mu
         Call Newtyp2(Irrepx,314,1,3,.True.)
         Call Newtyp2(Irrepx,315,2,4,.True.)
         Call Newtyp2(Irrepx,316,13,14,.True.)

         Call Newtyp2(Irrepx,317,1,3,.True.)
         Call Newtyp2(Irrepx,318,2,4,.True.)
         Call Newtyp2(Irrepx,319,13,14,.True.)
C Mu_tilde
         Call Newtyp2(Irrepx,324,1,3,.True.)
         Call Newtyp2(Irrepx,325,2,4,.True.)
         Call Newtyp2(Irrepx,326,13,14,.True.)

         Call Newtyp2(Irrepx,327,1,3,.True.)
         Call Newtyp2(Irrepx,328,2,4,.True.)
         Call Newtyp2(Irrepx,329,13,14,.True.)
C Mu_dot
         Call Newtyp2(Irrepx,334,1,3,.True.)
         Call Newtyp2(Irrepx,335,2,4,.True.)
         Call Newtyp2(Irrepx,336,13,14,.True.)
C Mu_tilde_dot
         Call Newtyp2(Irrepx,344,1,3,.True.)
         Call Newtyp2(Irrepx,345,2,4,.True.)
         Call Newtyp2(Irrepx,346,13,14,.True.)
      Endif

C Temporary storage of (ai,bi) quantities 
    
      If (Iuhf .EQ. 0) Then
         Call Newtyp2(Irrepx,302,9,10,.True.)
      Else
         Call Newtyp2(Irrepx,300,9,9,.True.)
         Call Newtyp2(Irrepx,301,10,10,.True.)
         Call Newtyp2(Irrepx,302,9,10,.True.)
      Endif

C Ao ladder and other auxilary lists 

      If (Iuhf .EQ. 0) Then
         If (Iflags(93) .Eq. 2) then 
            Call Newtyp2(Irrepx,377,15,14,.True.)
            Call Newtyp2(Irrepx,387,15,14,.True.)
         Endif 
      Else
         If (Iflags(93) .Eq. 2) then 
            Call Newtyp2(Irrepx,375,15,3,.True.)
            Call Newtyp2(Irrepx,376,15,4,.True.)
            Call Newtyp2(Irrepx,377,15,14,.True.)
            Call Newtyp2(Irrepx,385,15,3,.True.)
            Call Newtyp2(Irrepx,386,15,4,.True.)
            Call Newtyp2(Irrepx,387,15,14,.True.)
         Endif 
      Endif

C Other miscellaneous lists

      Call Newtyp2(Irrepx,304,9,9,.True.)
      Call Newtyp2(Irrepx,307,9,10,.True.)
      Call Newtyp2(Irrepx,309,11,12,.True.)
      Call Newtyp2(Irrepx,350,9,9,.True.)
      Call Newtyp2(Irrepx,351,10,10,.True.)
      Call Newtyp2(Irrepx,352,9,10,.True.)
      Call Newtyp2(Irrepx,353,11,12,.True.)
      Call Newtyp2(Irrepx,360,9,9,.True.)
      Call Newtyp2(Irrepx,361,10,10,.True.)
      Call Newtyp2(Irrepx,362,9,10,.True.)
      Call Newtyp2(Irrepx,363,11,12,.True.)

      If (Iuhf .NE. 0) Then
         Call Newtyp2(Irrepx,305,10,10,.True.)
         Call Newtyp2(Irrepx,306,10,9,.True.)
         Call Newtyp2(Irrepx,308,12,11,.True.)
      Endif

C Lists for MO dipole moment integrals (OO,OV and VV blocks).

      HHA_Length = Irpdpd(Irrepx,21)
      HHB_Length = Irpdpd(Irrepx,22)
      PPA_Length = Irpdpd(Irrepx,19)
      PPB_Length = Irpdpd(Irrepx,20)

      Call Aces_list_resize(1,370,HHA_Length)
      Call Aces_list_resize(1,371,PPA_Length)
      Call Aces_list_resize(1,372,PHA_Length)

      If (Iuhf .NE. 0) Then
         Call Aces_list_resize(2,370,HHB_Length)
         Call Aces_list_resize(2,371,PPB_Length)
         Call Aces_list_resize(2,372,PHB_Length)
      Endif

      Return
      End
