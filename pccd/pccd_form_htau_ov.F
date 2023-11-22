      Subroutine Pccd_form_htau_ov(Doo,Dvv,Dvo,Hov_pq,Hov_qp,Work,
     +                             Maxcor,Nocc,Nvrt,Noccsq,Nvrtsq,
     +                             Nvrtocc,Iuhf,Nonhf)

      Implicit DOuble Precision(A-H,O-Z)
      Logical Nonhf
      Logical Symmetry

      Dimension Doo(Nocc,Nocc),Dvv(Nvrt,Nvrt),Dvo(Nvrt,Nocc)
      Dimension Hov_pq(Nocc*Nvrt)
      Dimension Hov_qp(Nocc*Nvrt)
      Dimension Work(Maxcor)

#include "info.com"
#include "machsp.com"
#include "sym.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"

       Common /Sym/Symmetry

       Data Ione /1/

       Call Pccd_ov(Doo,Dvv,Dvo,Hov_pq,Work,Maxcor,Iuhf,Nonhf)
       Call Putrec(20,"JOBARC","REFGRDOV",Nt(1)*IIntfp,Hov_pq)

#ifdef _DEBUG_LVL0
      Write(6,*)
      write(6,"(a)") "The X(ov) contribution"
      Write(6,"(5(1x,F15.9))") (Hov_pq(i), i=1,Nt(1))
      call checksum("Hov_pq :",Hov_pq,Nocc*Nvrt)
#endif

       Return
       End
