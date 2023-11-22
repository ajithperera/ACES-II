      Subroutine Pccd_form_htau_1d(Htau_pq,Htau_qp,Work,Maxcor,
     +                             Hpq,Dpq,Dhf,Doo,Dvv,Dov,Nocc,
     +                             Nvrt,Nbas,Iuhf,E)
      
      Implicit Double Precision(A-H,O-Z)
      Integer A,B
      Logical Sym_packed 

      Dimension Hpq(Nbas,Nbas)
      Dimension Dpq(Nbas,Nbas)
      Dimension Dhf(Nbas,Nbas)
      Dimension Doo(Nocc,Nocc)
      Dimension Dvv(Nvrt,Nvrt)
      Dimension Dov(Nocc,Nvrt)
      Dimension Htau_pq(Nbas,Nbas)
      Dimension Htau_qp(Nbas,Nbas)
      Dimension Work(Maxcor)
      Dimension Ioffo(8)

#include "info.com"
#include "machsp.com"
#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Data Ione,One,Onem,DNull,Two,Half /1,1.0D0,-1.0D0,0.0D0,2.0D0,
     +                                   0.50D0/

C The expression coded here is H(p,r)D(r,q) - H(r,q)D(p,r). This is
C antisymmetric (from p+q)

#ifdef _NO_NRML_ORD

C Reorder the MO one particle Hamiltonian to correlated-order.
      I000 = Ione 
      I010 = I000 + Nbas
      Iend = I010 + Nbas*Nbas 
      Call Getrec(20,"JOBARC","REORDERA",Nbas,Work(I000))

#ifdef _DEBUG_LVLM
      Write(6,"(a)") "The reording array"
      Write(6,"(6(1x,i3))") (Work(I000+i-1),i=1,Nbas)
#endif 

      Call Dcopy(Nbas*Nbas,Hpq,1,Work(I010),1)
      Call Pccd_reordr(Work(I000),Work(I010),Hpq,Nbas,Nbas)

#ifdef _DEBUG_LVLM
      Write(6,"(a)") "The reordered one electron Hamiltonian"
      Call output(Hpq,1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif
#endif 

      I000 = Ione
      I010 = I000 + Nocc*Nocc
      I020 = I010 + Nvrt*Nvrt
      I030 = I020 + Nvrt*Nocc
      Iend = I030 + Nvrt*Nocc
      Maxcor = Maxcor - Iend 

      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d_1",Iend,
     +                                   Maxcor)

      Call Dgemm("T","T", Nbas,Nbas,Nbas,One,Hpq,Nbas,Dpq,Nbas,
     +            Dnull,Htau_qp,Nbas)

C Compute the One-Particle contribution to the energy. When normal-order
C this is one particle contribution to the correlation energy.

#ifdef _DEBUG_LVLM
      write(6,*)
      write(6,"(a)") "F(q,p) =  F(r,q)D(p,r) contribution" 
      call output(Htau_qp,1,Nbas,1,Nbas,Nbas,Nbas,1)
      write(6,*)

      E=dnull
      do i = 1, Nbas 
        E = E + Htau_qp(i,i)
      enddo 
      write(6,"(2a,(2x,F15.10))") "One particle contribution to the",
     +                            " energy :", E 
#endif 

      Call Analyze_fock(Work(I000),Work(I010),Work(I030),Nfmi(1),
     +                  Nfea(1),Nt(1),Nonhf)

      Call Pccd_form_fd(Doo,Dvv,Dov,Work(I000),Work(I010),Work(I020),
     +                  Work(I030),Work(Iend),Maxcor,Iuhf,Nonhf)

       Sym_packed = .False. 
       Call Pccd_frmful(Htau_pq,Work(I000),Work(I010),Work(I020),
     +                  Work(I030),Work(Iend),Maxcor,Nocc,Nvrt,
     +                  Nbas,"Ov_like",Sym_packed)

#ifdef _DEBUG_LVL0
      write(6,*)
      write(6,"(2a)") "F(p,q) =  F(p,r)D(r,q) for HF like Fock matrices"
     +                ," contributions added"
      call output(Htau_pq,1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif 

      Return
      End
   
