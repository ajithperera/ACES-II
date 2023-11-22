#include "flags.h"
      Subroutine Scrnc_form_lrspn_rcc_lists(Work,Memleft,Iuhf,Nperts)

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
      
      CALL Updmoi(1,MAX_PHA_length,1,480,0,0)
      CALL Updmoi(1,MAX_PHA_length,1,482,0,0)

      Call Aces_list_memset(1,480,0)
      Call Aces_list_memset(1,482,0)

      If (Iuhf .NE. 0) Then
         CALL Updmoi(1,MAX_PHB_length,2,480,0,0)
         CALL Updmoi(1,MAX_PHB_length,2,482,0,0)

         Call Aces_list_memset(2,480,0)
         Call Aces_list_memset(2,482,0)

      Endif
C 
      Do Irrep = 1, Nirrep 
          Length(Irrep) = Irpdpd(Irrep,9)
          If (Iuhf .NE. 0) Then
             Length(Irrep) = Length(Irrep) + Irpdpd(Irrep,9)
          Endif
          Call Updmoi(Nperts,Length(Irrep),Irrep,373,0,0)
          Call Updmoi(Nperts,Length(Irrep),Irrep,374,0,0)
          Call Updmoi(Nperts,Length(Irrep),Irrep,375,0,0)
          Call Updmoi(Nperts,Length(Irrep),Irrep,376,0,0)
      Enddo
         
      Return
      End 
 
