      Subroutine Pccd_form_htau_1d_2(Htau_pq,Htau_qp,Doo,Dvv,Dvo,
     +                               Fds,Work,Maxcor,Nocc,Nvrt,
     +                               Nbas,Iuhf)
      
      Implicit Double Precision(A-H,O-Z)
      Integer A,B
      Logical Sym_packed,Nonhf_ref

      Dimension Hpq(Nbas,Nbas)
      Dimension Doo(Nocc*Nocc)
      Dimension Dvv(Nvrt*Nvrt)
      Dimension Dvo(Nvrt*Nocc)
      Dimension Fds(Nbas*Nbas)
      Dimension Htau_pq(Nbas,Nbas)
      Dimension Htau_qp(Nbas,Nbas)
      Dimension Work(Maxcor)
      Dimension Ioffo(8)

      Common/Nhfref/Nonhf_ref

#include "sym.com"
#include "sympop.com"
#include "syminf.com"
#include "info.com"

      Data Ione,One,Onem,DNull,Two,Half /1,1.0D0,-1.0D0,0.0D0,2.0D0,
     +                                   0.50D0/

#ifdef _DEBUG_LVLM
      write(6,*)
      write(6,"(a)") "Fock diagonals"
      call output(Fds,1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif 
C Extract the occupied and virtual diagonal to two seperate blocks.

      I000 = Ione
      I010 = I000 + Nocc*Nocc
      I020 = I010 + Nvrt*Nvrt
      I030 = I020 + Nvrt*Nocc
      I040 = I030 + Nvrt*Nocc
      I050 = I040 + Nocc*Nocc
      I060 = I050 + Nvrt*Nvrt
      Iend = I060 + Nbas*Nbas
      If (Iend .Gt. Maxcor) Call Insmem("pccd_form_htau_1d_2",Iend,
     +                                   Maxcor)

      Call Dzero(Work(I000),Nocc*Nocc)
      Call Dzero(Work(I010),Nvrt*Nvrt)
      Call Dcopy(Nocc,Fds,Nbas+1,Work(I000),Nocc+1)
      Call Dcopy(Nvrt,Fds(Nocc*Nbas+Nocc+1),Nbas+1,Work(I010),Nvrt+1)

#ifdef _DEBUG_LVL0
      write(6,*)
      write(6,"(a)") "OO Fock diagonals"
      call output(Work(I000),1,Nocc,1,Nocc,Nocc,Nocc,1)
      Write(6,*)
      write(6,"(a)") "VV Fock diagonals"
      call output(Work(I010),1,Nvrt,1,Nvrt,Nvrt,Nvrt,1)
#endif 
      Call Dgemm("T","N",Nocc,Nvrt,Nvrt,One,Dvo,Nvrt,Work(I010),Nvrt,
     +            Dnull,Work(I020),Nocc)

      Call Dgemm("N","N",Nvrt,Nocc,Nocc,One,Dvo,Nvrt,Work(I000),Nocc,
     +            Dnull,Work(I030),Nvrt)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "G2(vo)=F(vv)*Dvo contribution"
      call output(Work(I020),1,Nocc,1,Nvrt,Nocc,Nvrt,1)
      Write(6,*)
      Write(6,"(a)") "G2(ov)=F(oo)*D(vo)(t) contribution"
      call output(Work(I030),1,Nvrt,1,Nocc,Nvrt,Nocc,1)
#endif 

      Call Pccd_proc_ovvo(Work(I020),Work(I030),Work(Iend),Maxcor,
     +                    Nocc,Nvrt)
      
      Call Dgemm("N","N",Nocc,Nocc,Nocc,One,Work(I000),Nocc,Doo,Nocc,
     +            Dnull,Work(I040),Nocc)
      Call Daxpy(Nocc,One,Work(I000),Nocc+1,Work(I040),Nocc+1)

      Call Dgemm("N","N",Nvrt,Nvrt,Nvrt,One,Work(I010),Nvrt,Dvv,Nvrt,
     +            Dnull,Work(I050),Nvrt)
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "G2(oo)=Dpo(t)*F(oo) contribution"
      call output(Work(I040),1,Nocc,1,Nocc,Nocc,Nocc,1)
      Write(6,*)
      Write(6,"(a)") "G2(vv)=Dvv*F(vv) contribution"
      call output(Work(I050),1,Nvrt,1,Nvrt,Nvrt,Nvrt,1)
#endif 
      Sym_packed = .False.
      Call Pccd_frmful(Work(I060),Work(I040),Work(I050),Work(I030),
     +                 Work(I020),Work(Iend),Maxcor,Nocc,Nvrt,Nbas,
     +                 "Ov_like",Sym_paked)
#ifdef _DEBUG_LVL0
      write(6,*)
      write(6,"(2a)") "F(p,q) =  F(p,r)D(r,q) contributions from the", 
     +                " diagonal of the Fock matrix"
      call output(Work(I060),1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif 

      Call Dscal(Nbas*Nbas,Onem,Htau_pq,1)
      Call Daxpy(Nbas*Nbas,One,Work(I060),1,Htau_pq,1)

      Do A = Nocc+1, Nbas
         Do I = 1, Nocc
            Index = (A-1)*Nbas + I
            Htau_pq(I,A) = -Htau_pq(I,A)
         Enddo
      Enddo

      Do I = 1, Nocc
         Do A = Nocc+1, Nbas
            Index = (I-1)*Nbas + A
            Htau_pq(A,I) = -Htau_pq(A,I)
         Enddo
      Enddo


#ifdef _DEBUG_LVL0
      write(6,*)
      write(6,"(2a)")"F(p,q) =  F(p,r)D(r,q) contributions from the", 
     +               " digonal of the Fock matrix added"
      write(6,"(a)") " to the contributions from the non diagonal Fock."
      call output(Htau_pq,1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif 

      Return
      End
   
