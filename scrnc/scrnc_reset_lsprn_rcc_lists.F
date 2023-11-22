#include "flags.h"
#include "disttype.h"

      Subroutine Scrnc_reset_lrspn_rcc_lists(Work,Memleft,Iuhf,Irrepx)

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

      If (Iuhf .NE. 0) Then
         Call Aces_list_resize(2,480,PHB_Length)
         Call Aces_list_resize(2,482,PHB_Length)
      Endif

      Return
      End
