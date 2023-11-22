      Subroutine Pccd_set2null(Work,Maxcor,IUhf,List)
     
      Implicit Double Precision(A-H,O-Z)
      Integer Dissiz
      
#include "syminf.com"
#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "sympop.com"
#include "info.com"

      Dimension Work(Maxcor)

      Data Ione /1/

      Irrepx  = Ione
      List2_h = List
      List2_t = List
      
      Do Irrep = 1, Nirrep

         Irrepr = Irrep
         Irrepl = Dirprd(Irrepr,Irrepx) 
            
          Numdis = Irpdpd(Irrepr,14)
          Dissiz = Irpdpd(Irrepl,15)

          I000 = Ione
          I010 = I000 + Numdis*Dissiz
          Iend = I010 + Numdis*Dissiz
          If (Iend .Gt. Maxcor) Call Insmem("transform_4index_ab",
     +                                       Iend,Maxcor)
          Call Dzero(Work(I000),Numdis*Dissiz)
          Call Putlst(Work(I000),1,Numdis,1,Irrepr,List2_t)     

#ifdef _DEBUG_LVL0
      call checksum("T2-ot   :",Work(I000),Numdis*Dissiz)
#endif
      Enddo 

      Return
      End
            
