      Subroutine Pccd_proc_htau(Htau_pq,Work,Maxcor,Nbas,Nocc,Nvrt)

      Implicit Double Precision(A-H,O-Z)

      Logical Sym_packed 
      Dimension Htau_pq(Nbas,Nbas)
      Dimension Work(Maxcor)

#include "machsp.com"
#include "sym.com"
#include "sympop.com"

      Data Ione,Izero,Onem /1,0,-1.0D0/

      Irrepx   = Ione
      Length   = Nbas*Nbas
      Lenoo    = Irpdpd(Irrepx,21)
      Lenvv    = Irpdpd(Irrepx,19)
      Lenvo    = Irpdpd(Irrepx,9)

      Nocc2 = Nocc*Nocc
      Nvrt2 = Nvrt*Nvrt
      Nvo2  = Nvrt*Nocc 

      I000 = Ione
      I010 = I000 + Nocc2
      I020 = I010 + Nvrt2
      I030 = I020 + Nvo2
      Iend = I030 + Nvo2
      Memleft = Maxcor - Iend

      Call Pccd_frmblocks(Htau_pq,Work(I000),Work(I010),Work(I020),
     +                    Work(I030),Work(Iend),Memleft,
     +                    Lenoo,Lenvv,Lenvo,Nocc,Nvrt,Nbas)

      I040 = Iend
      I050 = I040 + Nocc2
      Iend = I050 + Nvrt2
      Memleft = Maxcor - Iend

      Call Transp(Work(I000),Work(I040),Nocc,Nocc)
      Call Transp(Work(I010),Work(I050),Nvrt,Nvrt)
      Call Dscal(Nvrt*Nocc,Onem,Work(I030),1)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "G1(oo)=Dpo(t)*F(oo) contribution"
      call output(Work(I040),1,Nocc,1,Nocc,Nocc,Nocc,1)
      Write(6,*)
      Write(6,"(a)") "G1(vv)=Dvv*F(vv) contribution"
      call output(Work(I050),1,Nvrt,1,Nvrt,Nvrt,Nvrt,1)
#endif
      Sym_packed = .False.
      Call Pccd_putblocks(Htau_pq,Work(I040),Work(I050),Work(I020),
     +                    Work(I030),Work(Iend),Memleft,
     +                    Lenoo,Lenvv,Lenvo,Nocc,Nvrt,Nbas,
     +                    "Vo_like",Sym_packed)

      Return
      End 
