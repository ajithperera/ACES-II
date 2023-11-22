      Subroutine Scrnc_rccd_dummy_inits(Work,Maxmem,Iuhf)
      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Maxmem)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "sym.com"
#include "flags2.com"
#include "flags.com"
#include "info.com"

C Rccd and Drccd has no singles. Since all the routines that have
C been developed earier assumed CCSD (hence singles). Lets create
C singles list and set them to zero, so all of them still work
C correctly (cost consideration are minuscule). 

      Do Ispin = 1, Iuhf+1
         Call Updmoi(1,Nt(Ispin),Ispin,90,0,0)
         Call Updmoi(1,Nt(Ispin),Ispin,190,0,0)
         Call Updmoi(1,Nfmi(Ispin),Ispin,191,0,0)
         Call Updmoi(1,Nfea(Ispin),Ispin,192,0,0)

         Call Aces_list_memset(Ispin,90,0)
         Call Aces_list_memset(Ispin,190,0)
         Call Aces_list_memset(Ispin,191,0)
         Call Aces_list_memset(Ispin,192,0)
      Enddo 

      Return
      End 
