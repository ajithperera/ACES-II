      Subroutine Transform_4index_ab(T2amps,Coo,Cvv,Work,Maxcor,List2_h,
     +                               List2_t,T2ln,T2ln_aa,T2ln_bb,Iuhf,
     +                               Cct2,Type)
      
      Implicit Double Precision(A-H,O-Z)
      Integer Dissiz
      Integer T2ln,T2off,T2ln_aa,T2ln_bb
      Character*3 Type
      Logical CCt2
      
#include "syminf.com"
#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "sympop.com"
#include "info.com"

      Dimension Coo(Nfmi(1)+Iuhf*Nfmi(2))
      Dimension Cvv(Nfea(1)+Iuhf*Nfea(2))
      Dimension Work(Maxcor)
      Dimension T2amps(T2ln)

      Data Ione /1/

      Irrepx = Ione
      T2off  = Ione + T2ln_aa + T2ln_bb
      
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

          If (Cct2) Then
              Call Getlst(Work(I000),1,Numdis,1,Irrepr,List2_h)
          Else 
              Call Dcopy(Numdis*Dissiz,T2amps(T2off),1,Work(I000),1)
          Endif 

#ifdef _DEBUG_LVL0
      call checksum("T2-in:",Work(I000),Numdis*Dissiz)
#endif
          Do I = 1, Numdis

             Ioff = (I-1)*Dissiz 
             Call Trans_vv_ab(Work(I000+Ioff),Work(I010+Ioff),Cvv,
     +                        Dissiz,Irrepl,Iuhf,Type)
          Enddo 

C T2(Ab,Ij)->T2(Ij,Ab)

          Call Transp(Work(I000),Work(I010),Numdis,Dissiz)

          do I = 1, Dissiz

             Joff = (I-1)*Numdis 
             Call Trans_oo_ab(Work(I010+Joff),Work(I000+Joff),Coo,
     +                        Numdis,Irrepr,Iuhf,Type)
          Enddo
            
C T2(Ij,Ab)->T2(Ab,Ij)

          Call Transp(Work(I010),Work(I000),Dissiz,Numdis)

          If (CCt2) Then
             Call Putlst(Work(I000),1,Numdis,1,Irrepr,List2_t)     
          Else 
             Call Dcopy(Numdis*Dissiz,Work(I000),1,T2amps(T2off),1)
          Endif 

#ifdef _DEBUG_LVL0
      call checksum("T2-ot:",Work(I000),Numdis*Dissiz)
#endif
          T2off = T2off + Dissiz*Numdis

      Enddo 

      Return
      End
            
