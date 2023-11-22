#include "flags.h"
#include "disttype.h"

      Subroutine Scrnc_reset_lrspn_cc_lists(Work,Memleft,Iuhf,Irrepx)

      Implicit Double Precision (A-H, O-Z)
      Integer AAAA_LENGTH_IJAB,BBBB_LENGTH_IJAB,AABB_LENGTH_IJAB
      Integer PHA_Length,PHB_Length,HHA_length,HHB_Length
      Integer PPA_Length,PPB_Length
      Integer Ttyper,Ttypel
      Logical Ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,Drccd

      Dimension Work(Memleft)
      Common /Extrap/Maxexp,Nreduce,Ntol,Nsizec
      Common /Eominfo/Ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,
     +                Drccd

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "parcc.com"
#include "flags.com"

C 490 (Ispin=1,2) keeps Mbar(a,i) elements


      PHA_Length = Irpdpd(Irrepx,9)
      PHB_Length = Irpdpd(Irrepx,10)
      PPA_Length = Irpdpd(Irrepx,19)
      PPB_Length = Irpdpd(Irrepx,20)
      HHA_Length = Irpdpd(Irrepx,21)
      HHB_Length = Irpdpd(Irrepx,22)
      
      Call Aces_list_resize(1,480,PHA_Length)
      Call Aces_list_resize(1,482,PHA_Length)
      Call Aces_list_resize(9,448,PHA_Length)
      Call Aces_list_resize(1,490,PHA_Length)
      Call Aces_list_resize(3,490,PHA_Length)
      Call Aces_list_resize(1,493,PHA_Length)
      Call Aces_list_resize(1,491,HHA_Length)
      Call Aces_list_resize(1,492,PPA_Length)


      If (Iuhf .NE. 0) Then
         Call Aces_list_resize(2,480,PHB_Length)
         Call Aces_list_resize(2,482,PHB_Length)
         Call Aces_list_resize(9,449,PHB_Length)
         Call Aces_list_resize(2,490,PHB_Length)
         Call Aces_list_resize(4,490,PHB_Length)
         Call Aces_list_resize(2,493,PHB_Length)
         Call Aces_list_resize(1,491,HHB_Length)
         Call Aces_list_resize(2,492,PPB_Length)
      Endif

C Doubles lists for (ab,ij). 

      If (Iuhf .EQ. 0) Then
         Call Newtyp2(Irrepx,446,13,14,.True.)
         Call Newtyp2(Irrepx,456,13,14,.True.)
      Else
         Call Newtyp2(Irrepx,444,1,3,.True.)
         Call Newtyp2(Irrepx,445,2,4,.True.)
         Call Newtyp2(Irrepx,446,13,14,.True.)

         Call Newtyp2(Irrepx,454,1,3,.True.)
         Call Newtyp2(Irrepx,455,2,4,.True.)
         Call Newtyp2(Irrepx,456,13,14,.True.)
      Endif

C Temporary storage of (ai,bi) quantities 
    
      If (Iuhf .EQ. 0) Then
         Call Newtyp2(Irrepx,442,9,10,.True.)
      Else
         Call Newtyp2(Irrepx,440,9,9,.True.)
         Call Newtyp2(Irrepx,441,10,10,.True.)
         Call Newtyp2(Irrepx,442,9,10,.True.)
      Endif

C Q(ab) and Q(ij) three-body intermediates 

      PPA_Length = Irpdpd(Irrepx,19)
      PPB_Length = Irpdpd(Irrepx,20)
      HHA_Length = Irpdpd(Irrepx,21)
      HHB_Lenght = Irpdpd(Irrepx,22)
  
      Call Aces_list_resize(1,492,PPA_Length)
      Call Aces_list_resize(1,491,HHA_Length)
      If (Iuhf .Ne. 0) Then 
         Call Aces_list_resize(2,492,PPB_Length)
         Call Aces_list_resize(2,491,HHB_Length)
      Endif 

C Davidson lists 

      Len=0
      Do Ispin = 1, Iuhf+1
         Len = Len + Irpdpd(Irrepx,8+Ispin)
      Enddo
      Len = Len + Idsymsz(Irrepx,Isytyp(1,46),Isytyp(2,46))
      If (Iuhf .Ne.0) Then
          Len = Len + Idsymsz(Irrepx,Isytyp(1,44),Isytyp(2,44))
          Len = Len + Idsymsz(Irrepx,Isytyp(1,45),Isytyp(2,45))
      Endif
      Call Aces_list_resize(1,470,Len)
      Call Aces_list_resize(1,471,Len)
      Call Aces_list_resize(1,472,Len)
      Call Aces_list_resize(2,470,Len)
      Call Aces_list_resize(2,471,Len)
      Call Aces_list_resize(2,472,Len)

C Denominator and Ao-ladder lists

      Do Ispin=3,3-2*Iuhf,-1
         Ttypel=Isytyp(1,43+Ispin)
         TtypeR=Isytyp(2,43+Ispin)
         Call Newtyp2(Irrepx,443+Ispin,Ttypel,Ttyper,.True.)
         Call Newtyp2(Irrepx,447+Ispin,Ttypel,Ttyper,.True.)
         Call Newtyp2(Irrepx,460+Ispin,Ttypel,Ttyper,.True.)
         If (AOladder) Then
            Ttypel=15
            Call Newtyp2(Irrepx,413+Ispin,Ttypel,Ttyper,.True.)
            Call Newtyp2(Irrepx,463+Ispin,Ttypel,Ttyper,.True.)
         Endif
      Enddo

C Resorted R2 lisits 

      If (Iuhf .eq. 0) Then
         Do Ilist=37,43,2
            Ttypel=Isytyp(1,Ilist)
            TtypeR=Isytyp(2,Ilist)
            Call Newtyp2(Irrepx,400+Ilist,Ttypel,Ttyper,.True.)
         Enddo
         Ttypel=Isytyp(1,42)
         TtypeR=Isytyp(2,42)
         Call Newtyp2(Irrepx,442,Ttypel,Ttyper,.True.)
      Else
         Do Ilist=34,43
            Ttypel=Isytyp(1,Ilist)
            TtypeR=Isytyp(2,Ilist)
            Call Newtyp2(Irrepx,400+Ilist,Ttypel,Ttyper,.True.)
         Enddo
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

C T2 resorted and other ring like lists 

      CALL NEWTYP2(IRREPX,434,h_DistType_AI,h_DistType_AI,.TRUE.)
      CALL NEWTYP2(IRREPX,435,h_DistType_ai,h_DistType_ai,.TRUE.)
      CALL NEWTYP2(IRREPX,436,h_DistType_ai,h_DistType_AI,.TRUE.)
      CALL NEWTYP2(IRREPX,437,h_DistType_AI,h_DistType_ai,.TRUE.)
      CALL NEWTYP2(IRREPX,438,h_DistType_aI,h_DistType_Ai,.TRUE.)
      CALL NEWTYP2(IRREPX,439,h_DistType_Ai,h_DistType_aI,.TRUE.)
      CALL NEWTYP2(IRREPX,440,h_DistType_AI,h_DistType_AI,.TRUE.)
      CALL NEWTYP2(IRREPX,441,h_DistType_ai,h_DistType_ai,.TRUE.)
      CALL NEWTYP2(IRREPX,442,h_DistType_AI,h_DistType_ai,.TRUE.)
      CALL NEWTYP2(IRREPX,443,h_DistType_Ai,h_DistType_aI,.TRUE.)


      Return
      End
