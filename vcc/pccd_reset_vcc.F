      Subroutine Pccd_reset_vcc(Work,Maxcor,IUhf,List)
     
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

      Nsize_aa = Isymsz(Isytyp(1,61),Isytyp(2,61))
      Call Getall(Work,Nsize_aa,Irrepx,61)
      Call Dzero(Work,Nsize_aa)
      Call Putall(Work,Nsize_aa,Irrepx,61)

      If (Iuhf .Ne. 0) Then
         Nsize_bb = Isymsz(Isytyp(1,62),Isytyp(2,62))
         Call Getall(Work,Nsize_bb,Irrepx,62)
         Call Dzero(Work,Nsize_bb)
         Call Putall(Work,Nsize_bb,Irrepx,62)
      Endif 

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
C T2(Ab,Ij)
              Call Getlst(Work(I000),1,Numdis,1,Irrepr,List2_h)

#ifdef _DEBUG_LVL0
      call checksum("T2-in   :",Work(I000),Numdis*Dissiz)
#endif
          Do I = 1, Numdis

             Ioff = (I-1)*Dissiz 
             Call pccd_reset_vv_ab(Work(I000+Ioff),Dissiz,Irrepl)
          Enddo 

C T2(Ab,Ij)->T2(Ij,Ab)

          Call Transp(Work(I000),Work(I010),Numdis,Dissiz)

          do I = 1, Dissiz

             Joff = (I-1)*Numdis 
              Call pccd_reset_oo_ab(Work(I010+Joff),Numdis,Irrepr)
          Enddo
            
C T2(Ij,Ab)->T2(Ab,Ij)

          Call Transp(Work(I010),Work(I000),Dissiz,Numdis)

          Call Putlst(Work(I000),1,Numdis,1,Irrepr,List2_t)     

#ifdef _DEBUG_LVL0
      call checksum("T2-ot   :",Work(I000),Numdis*Dissiz)
#endif
      Enddo 

      Return
      End
            
