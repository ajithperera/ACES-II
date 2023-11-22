      Subroutine Form_new_t2pamps_ab(T2amps,Work,Iwork,Maxcor,Imaxcor,
     +                               List2,T2ln,T2ln_aa,T2ln_bb,
     +                               Ioffset,Irrepx,Iuhf)
      
      Implicit Double Precision(A-H,O-Z)
      Integer Dissiz
      Integer T2ln,T2off,T2ln_aa,T2ln_bb
      
#include "syminf.com"
#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "sympop.com"
#include "info.com"

      Dimension Work(Maxcor)
      Dimension T2amps(T2ln)
      Dimension Iwork(Imaxcor)

      Data Ione /1/

      T2off  = Ione + T2ln_aa + T2ln_bb
      
      Do Irrep = 1, Nirrep

         Irrepr = Irrep
         Irrepl = Dirprd(Irrepr,Irrepx) 
            
          Numdis = Irpdpd(Irrepr,14)
          Dissiz = Irpdpd(Irrepl,13)

          I000 = Ione
          I010 = I000 + Numdis*Dissiz
          I020 = I010 + Numdis*Dissiz
          Iend = I020 + Numdis*Dissiz
          If (Iend .Gt. Maxcor) Call Insmem("form_new_amps_ab",
     +                                       Iend,Maxcor)

C T2(Ab,Ij)

          Call Getlst(Work(I000),1,Numdis,1,Irrepr,List2)
          Call Dcopy(Numdis*Dissiz,T2amps(T2off),1,Work(I010),1)

#ifdef _DEBUG_LVLM
      call checksum("T2-in:",Work(I000),Numdis*Dissiz)
      call checksum("T2-in:",Work(I010),Numdis*Dissiz)
#endif
          Call Form_pvvoo_ab(Work(I000),Work(I010),Work(I020),
     +                       Iwork,Imaxcor,Numdis,Dissiz,Irrepr,
     +                       Irrepl,Ioffset)

#ifdef _DEBUG_LVLM
      call checksum("T2-ot:",Work(I000),Numdis*Dissiz)
#endif 
          Call Putlst(Work(I000),1,Numdis,1,Irrepr,List2)     

          T2off = T2off + Dissiz*Numdis

      Enddo 

      Return
      End
            
