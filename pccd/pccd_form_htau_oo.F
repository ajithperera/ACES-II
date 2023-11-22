      Subroutine Pccd_form_htau_oo(Doo,Dvv,Dvo,Hoo_pq,Hoo_qp,Work,
     +                             Maxcor,Nocc,Nvrt,Noccsq,Nvrtsq,
     +                             Nvrtocc,Iuhf,Nonhf)

      Implicit DOuble Precision(A-H,O-Z)
      Logical Non_hf_terms 
      Logical Nonhf
      Logical Symmetry

      Dimension Doo(Nocc,Nocc),Dvv(Nvrt,Nvrt),Dvo(Nvrt,Nocc)
      Dimension Hoo_pq(Nocc*Nocc),Hoo_qp(Nocc*Nocc)
      Dimension Work(Maxcor)

#include "info.com"
#include "machsp.com"
#include "sym.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
       Common /Sym/Symmetry

       Data Ione /1/

       I000 = Ione
       I010 = I000 + Nfmi(1)+Iuhf*Nfmi(2)
       I020 = I010 + Nfea(1)+Iuhf*Nfea(2)
       Iend = I020 + Nt(1)+Iuhf*Nt(2)
       If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_oo",Iend,
     +                                    Maxcor)

       Call Getrec(20,"JOBARC","DENSOO  ",Noccsq*IIntfp,Doo)
       Call Getrec(20,"JOBARC","DENSVV  ",Nvrtsq*IIntfp,Dvv)
       If (Nonhf) Call Getrec(20,"JOBARC","DENSVO  ",
     +                        Nvrtocc*IIntfp,Dvo)

       Call Analyze_Fock(Work(I000),Work(I010),Work(I020),Nfmi(1),
     +                   Nfea(1),Nt(1),Non_hf_terms)

       Call Pccd_oo(Doo,Dvv,Dvo,Hoo_pq,Work,Maxcor,Iuhf,(Non_hf_terms
     +              .or.Nonhf))
#ifdef _NOSKIP
       Irrepx = Ione
       Ispin  = Ione
       I000 = Ione
       I010 = I000 + Noccsq
       I020 = I010 + Noccsq
       Iend = I020 + Noccsq
       Call Dzero(Work(I000),Noccsq)
       Call Pccd_blockcopy(Work(I000),Hoo_pq,"OO",Work(I010),
     +                     Work(I020),Nocc,Irpdpd(Irrepx,20+Ispin),
     +                     Noccsq,Ispin)
       Call Dzero(Hoo_pq,Noccsq)
       Call Dcopy(Noccsq,Work(I000),1,Hoo_qp,1)
#endif 

       Call Putrec(20,"JOBARC","REFGRDOO",Nfmi(1)*IIntfp,Hoo_pq)

#ifdef _DEBUG_LVL0
CSSS      call pccd_check_htau("Htau_oo :",Hoo_pq,Nocc,"OO","D")
      Write(6,*)
      write(6,"(a)") "The oo contribution from ref. state"
      Write(6,"(5(1x,F15.9))") (Hoo_pq(i), i=1,Nfmi(1))
      call checksum("Hoo_pq :",Hoo_pq,Nocc*Nocc)
#endif

       Return
       End
