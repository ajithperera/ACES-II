#include "flags.h"
      Subroutine Scrnc_form_lrspn_cc_lists(Work,Memleft,Iuhf,Nperts)

      Implicit Double Precision (A-H, O-Z)
      Integer AAAA_LENGTH_IJAB,BBBB_LENGTH_IJAB,AABB_LENGTH_IJAB
      Integer PHA_Length,PHB_Length,HHA_length,HHB_Length
      Integer PPA_Length,PPB_Length
      Integer Ttypel,Ttyper
      Integer Length(8)
      Logical Ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,Drccd

      Dimension Work(Memleft)

      Common /Extrap/Maxexp,Nreduce,Ntol,Nsizec
      Common /Eominfo/Ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,Drccd

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "parcc.com"
#include "flags.com"

       Imode = 0

C 490 (Ispin=1,2) keeps Mbar(a,i) elements 

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
      
      CALL Updmoi(1,MAX_PHA_length,1,490,0,0)
      CALL Updmoi(1,MAX_PHA_length,3,490,0,0)
      CALL Updmoi(1,MAX_HHA_length,1,491,0,0)
      CALL Updmoi(1,MAX_PPA_length,1,492,0,0)
      CALL Updmoi(1,MAX_PHA_length,1,480,0,0)
      CALL Updmoi(1,MAX_PHA_length,1,482,0,0)
      CALL Updmoi(1,MAX_PHA_length,1,493,0,0)
      CALL Updmoi(1,MAX_PHA_length,9,448,0,0)

      Call Aces_list_memset(1,480,0)
      Call Aces_list_memset(1,482,0)
      Call Aces_list_memset(9,448,0)
      Call Aces_list_memset(1,490,0)
      Call Aces_list_memset(3,490,0)
      Call Aces_list_memset(1,493,0)
      Call Aces_list_memset(1,491,0)
      Call Aces_list_memset(1,492,0)

      If (Iuhf .NE. 0) Then
         CALL Updmoi(1,MAX_PHB_length,2,490,0,0)
         CALL Updmoi(1,MAX_PHB_length,4,490,0,0)
         CALL Updmoi(1,MAX_HHB_length,2,491,0,0)
         CALL Updmoi(1,MAX_PPB_length,2,492,0,0)
         CALL Updmoi(1,MAX_PHB_length,2,480,0,0)
         CALL Updmoi(1,MAX_PHB_length,2,482,0,0)
         CALL Updmoi(1,MAX_PHB_length,9,449,0,0)
         CALL Updmoi(1,MAX_PHB_length,2,493,0,0)

         Call Aces_list_memset(2,480,0)
         Call Aces_list_memset(2,482,0)
         Call Aces_list_memset(9,449,0)
         Call Aces_list_memset(2,490,0)
         Call Aces_list_memset(4,490,0)
         Call Aces_list_memset(2,493,0)
         Call Aces_list_memset(2,491,0)
         Call Aces_list_memset(2,492,0)
      Endif

C The doubles lists for (ab,ij). At this instance left-hand side 
C is arbitrary. So irrep_arb can be any value of 1-8.

      Irrep_arb = 1
      If (Iuhf .EQ. 0) Then
         Call Inipck2(Irrep_arb,13,14,446,Imode,0,1)

         Call Inipck2(Irrep_arb,13,14,456,Imode,0,1)
      Else
         If (Rccd .OR. Drccd) Then
            Call Inipck2(Irrep_arb,19,21,444,Imode,0,1)
            Call Inipck2(Irrep_arb,20,22,445,Imode,0,1)
            Call Inipck2(Irrep_arb,19,21,454,Imode,0,1)
            Call Inipck2(Irrep_arb,20,22,455,Imode,0,1)
         Else
            Call Inipck2(Irrep_arb,1,3,444,Imode,0,1)
            Call Inipck2(Irrep_arb,2,4,445,Imode,0,1)
            Call Inipck2(Irrep_arb,1,3,454,Imode,0,1)
            Call Inipck2(Irrep_arb,2,4,455,Imode,0,1)
         Endif 

         Call Inipck2(Irrep_arb,13,14,446,Imode,0,1)

         Call Inipck2(Irrep_arb,13,14,456,Imode,0,1)
      Endif

C Temporary storage of (ai,bi) quantities 
    
      If (Iuhf .EQ. 0) Then
         Call Inipck2(Irrep_arb,9,10,442,Imode,0,1)
      Else
         Call Inipck2(Irrep_arb,9,9,440,Imode,0,1)
         Call Inipck2(Irrep_arb,10,10,441,Imode,0,1)
         Call Inipck2(Irrep_arb,9,10,442,Imode,0,1)
      Endif
      
C Create list of Q(ab) and Q(ij) three-body intermediates.

      Do Ispin = 1, Iuhf+1
         Maxab=0
         Maxij=0
         Do Irrep =1, Nirrep 
            Maxab = Max(Maxab,Irpdpd(Irrep,18+Ispin))
            Maxij = Max(Maxij,Irpdpd(Irrep,20+Ispin))
         Enddo
         Call Updmoi(1,Maxab,Ispin,492,Imode,0)
         Call Updmoi(1,Maxij,Ispin,491,Imode,0)
         Call Updmoi(1,Maxab,Ispin+2,492,Imode,0)
         Call Updmoi(1,Maxij,Ispin+2,491,Imode,0)
      Enddo 

C Davidson lists 

      MaxLen=0
      Do Irrepx = 1, Nirrep
         Len=0
         Do Ispin = 1, Iuhf+1
            Len = Len + Irpdpd(Irrepx,8+Ispin)
         Enddo
         Len = Len + Idsymsz(Irrepx,Isytyp(1,46),Isytyp(2,46))
         If (Iuhf .Ne.0) Then
             Len = Len + Idsymsz(Irrepx,Isytyp(1,44),Isytyp(2,44))
             Len = Len + Idsymsz(Irrepx,Isytyp(1,45),Isytyp(2,45))
         Endif
         MaxLen=Max(Len,Maxlen)
      Enddo 

      Numlst=2
      Do I = 1, Numlst
         If (.NOT. Nodavid) Then
            Call Updmoi(Maxexp,Maxlen,I,470,Imode,0)
            Call Updmoi(Maxexp,Maxlen,I,471,Imode,0)
         Endif 
         Call Updmoi(3+Nreduce,Maxlen,I,472,Imode,0)
      Enddo
C         
      Do Ispin=3,3-2*Iuhf,-1

C Denominator lists 

         Ttypel=Isytyp(1,43+Ispin)
         TtypeR=Isytyp(2,43+Ispin)
         Call Inipck2(1,Ttypel,Ttyper,460+Ispin,Imode,0,1)
         Call Inipck2(1,Ttypel,Ttyper,447+Ispin,Imode,0,1)
CSSS         Call Inipck2(1,Ttypel,Ttyper,443+Ispin,Imode,0,1)

C Ao ladder lists 

         If (AOladder) Then
            Ttypel=15
            Call Inipck2(1,Ttypel,Ttyper,413+Ispin,Imode,0,1)
            Call Inipck2(1,Ttypel,Ttyper,463+Ispin,Imode,0,1) 
         Endif 
      Enddo 

C Resorted R2 lists 

      If (Iuhf .eq. 0) Then
         Do Ilist=37,43,2
            Ttypel=Isytyp(1,Ilist)
            TtypeR=Isytyp(2,Ilist)
            Call Inipck2(1,Ttypel,Ttyper,400+Ilist,Imode,0,1)
            Do Irrep =1, Nirrep 
               Call aces_list_memset(Irrep,400+Ilist,0)
            Enddo
         Enddo
         Call Inipck2(1,Isytyp(1,42),Isytyp(2,42),442,Imode,0,1)
         Do Irrep =1, Nirrep 
            Call aces_list_memset(Irrep,442,0)
         Enddo
       Else
          Do Ilist=34,43
             Ttypel=Isytyp(1,Ilist)
             TtypeR=Isytyp(2,Ilist)
             Call Inipck2(1,Ttypel,Ttyper,400+Ilist,Imode,0,1)
             Do Irrep =1, Nirrep 
               Call aces_list_memset(Irrep,400+Ilist,0)
             Enddo
          Enddo
       Endif
C 
       Do Irrep = 1, Nirrep 
          Length(Irrep) = Irpdpd(Irrep,9)
          If (Iuhf .EQ. 0) Then
             Length(Irrep) = Length(Irrep) + Idsymsz(Irrep,13,14)
          Else
             Length(Irrep) = Length(Irrep) + Irpdpd(Irrep,9)
             Length(Irrep) = Length(Irrep) + Idsymsz(Irrep,1,3)
             Length(Irrep) = Length(Irrep) + Idsymsz(Irrep,2,4)
             Length(Irrep) = Length(Irrep) + Idsymsz(Irrep,13,14)
          Endif
          Call Updmoi(Nperts,Length(Irrep),Irrep,373,0,0)
          Call Updmoi(Nperts,Length(Irrep),Irrep,374,0,0)
          Call Updmoi(Nperts,Length(Irrep),Irrep,375,0,0)
          Call Updmoi(Nperts,Length(Irrep),Irrep,376,0,0)
       Enddo
       Write(6,*) "Final"
         
       Return
       End 


 
