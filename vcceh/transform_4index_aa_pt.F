      Subroutine Transform_4index_aa_pt(T2amps,Coo,Cvv,Work,Maxcor,
     +                                  List2_h,List2_t,T2ln,T2ln_aa, 
     +                                  Iuhf,Irrepx,Cct2,Type)
      
      Implicit Double Precision(A-H,O-Z)
      Integer Dissiz_trn,Dissiz_sqr
      Integer T2ln,T2ln_aa,T2off
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

      Do Ispin =1, 1+Iuhf

         T2off =  (Ispin-1)*T2ln_aa + Ione

         Do Irrep = 1, Nirrep

            Irrepr = Irrep
            Irrepl = Dirprd(Irrepr,Irrepx) 

            Numdis_trn = Irpdpd(Irrepr,2+Ispin)
            Dissiz_trn = Irpdpd(Irrepl,Ispin)
            Numdis_sqr = Irpdpd(Irrepr,20+Ispin)
            Dissiz_sqr = Irpdpd(Irrepl,18+Ispin)
            
            I000 = Ione
            I010 = I000 + Numdis_sqr*Dissiz_sqr
            Iend = I010 + Numdis_sqr*Dissiz_sqr
            If (Iend .Gt. Maxcor) Call Insmem("transform_4index_aa",
     +                                         Iend,Maxcor)
C T2(A<B,I<J)
            If (Cct2) Then
                Call Getlst(Work(I000),1,Numdis_trn,1,Irrepr,
     +                      List2_h+Ispin)
            Else
                Call Dcopy(Dissiz_trn*Numdis_trn,T2amps(T2off),Ione,
     +                      Work(I000),Ione)
            Endif
#ifdef _DEBUG_LVLM
      call checksum("T2-in:",Work(I000),Numdis_trn*Dissiz_trn)
#endif

C T2(A<B,I<J)->T2(A<B,IJ)

            Call Symexp(Irrepr,Pop(1,Ispin),Dissiz_trn,Work(I000))

C T2(A<B,IJ)->T2(IJ,A<B)

            Call Transp(Work(I000),Work(I010),Numdis_sqr,
     +                  Dissiz_trn)

            Do I = 1,  Dissiz_trn
               
               Ioff =  (I-1)*Numdis_sqr 

               Call Trans_oo_aa(Work(I010+Ioff),Work(I000+Ioff),Coo,
     +                          Numdis_sqr,Irrepr,Iuhf,Ispin,Type)
            Enddo 

C T2(IJ,A<B)->T2(IJ,AB)

            Call Symexp(Irrepl,Vrt(1,Ispin),Numdis_sqr,Work(I010))

C T2(IJ,AB)->T2(AB.IJ)

            Call Transp(Work(I010),Work(I000),Dissiz_sqr,
     +                  Numdis_sqr)

            Do I = 1, Numdis_sqr

               Joff = (I-1)*Dissiz_sqr 
               Call Trans_vv_aa(Work(I000+Joff),Work(I010+Joff),Cvv,
     +                          Dissiz_sqr,Irrepl,Iuhf,Ispin,Type)
            Enddo

C T2(AB,IJ)->T2(A<B,IJ)
            Call Sqsym(Irrepl,Vrt(1,Ispin),Dissiz_trn,Dissiz_sqr,
     +                 Numdis_sqr,Work(I010),Work(I000))

C T2(A<B,IJ)->T2(IJ,A<B)

            Call Transp(Work(I010),Work(I000),Numdis_sqr,
     +                   Dissiz_trn)

C T2(IJ,AB)->T2(I<J,A<B)

            Call Sqsym(Irrepr,Pop(1,Ispin),Numdis_trn,Numdis_sqr,
     +                 dissiz_trn,Work(I010),Work(I000))

C T2(I<J,A<B)->T2(A<B,I<J)

            Call Transp(Work(I010),Work(I000),Dissiz_trn,
     +                  Numdis_trn)

            If (Cct2) Then
                Call Putlst(Work(I000),1,Numdis_trn,1,Irrepr,
     +                      List2_t+Ispin)
            Else
               Call Dcopy(Dissiz_trn*Numdis_trn,Work(I000),Ione,
     +                      T2amps(T2off),Ione)
            Endif 

            T2off = T2off + Dissiz_trn*Numdis_trn
#ifdef _DEBUG_LVLM
      call checksum("T2-ot:",Work(I000),Numdis_trn*Dissiz_trn)
#endif
        Enddo
      Enddo 

      Return
      End
            
