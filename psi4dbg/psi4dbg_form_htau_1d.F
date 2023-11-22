      Subroutine Psi4dbg_form_htau_1d(Htau_pq,Htau_qp,Work,Maxcor,
     +                                Hpq,Dpq,Nocc,Nvrt,Nbas,E)
      
      Implicit Double Precision(A-H,O-Z)

      Dimension Hpq(Nbas,Nbas)
      Dimension Dpq(Nbas,Nbas)
      Dimension Htau_pq(Nbas,Nbas)
      Dimension Work(Maxcor)
      Dimension Ioffo(8)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"
#include "info.com"

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

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The reordered one electron Hamiltonian"
      Call output(Hpq,1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif

#ifdef _NOSKIP 
C This block was originally used to built HF density  for One electron 
C Hamiltonian in SCF ordering). Instead, Now the one electron 
C Hamiltonian is reordered. 

      Ioffo(1) = 0
      Do Irrep =2, Nirrep
         Ioffo(Irrep)=Ioffo(Irrep-1)+(Pop(irrep-1,1)+Vrt(Irrep-1,1))
      Enddo
      Call Dzero(Dpq,Nbas*Nbas)
      Do Irrep = 1, Nirrep
         N = Pop(Irrep,1)
         Do I = 1, N
            Dpq(Ioffo(Irrep)+I,Ioffo(Irrep)+I) =  Two
         Enddo 
      Enddo
#endif 
#endif 

      I000 = Ione
      I010 = I000 + Nbas*Nbas
      Iend = I010 + Nbas*Nbas

      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d",Iend,Maxcor)

      Call Dgemm("N","N", Nbas,Nbas,Nbas,One,Hpq,Nbas,Dpq,Nbas,
     +            Dnull,Htau_pq,Nbas)

      Call Dgemm("T","T", Nbas,Nbas,Nbas,One,Hpq,Nbas,Dpq,Nbas,
     +            Dnull,Htau_qp,Nbas)

C Compute the One-Particle contribution to the energy. When normal-order
C this is one particle contribution to the correlation energy.

#ifdef _DEBUG_LVL0
      write(6,*)
      write(6,"(a)") "F(p,q) =  F(p,r)D(r,q) contribution" 
      call output(Htau_pq,1,Nbas,1,Nbas,Nbas,Nbas,1)
      write(6,*)
      write(6,"(a)") "F(q,p) =  F(r,q)D(p,r) contribution" 
      call output(Htau_qp,1,Nbas,1,Nbas,Nbas,Nbas,1)
      write(6,*)
#endif 

      E=dnull
      do i = 1, Nbas 
        E = E + Htau_pq(i,i)
      enddo 
      write(6,"(2a,(2x,F15.10))") "One particle contribution to the",
     +                            " MBPT(2) correlation energy :", 
     +                              E

#ifdef _NOSKIP
      I000 = Ione
      I010 = I000 + Nbas*Nbas
      Iend = I010 + Nbas*Nbas

      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_1d",Iend,Maxcor)

      Call Dgemm("N","N", Nbas,Nbas,Nbas,One,Hpq,Nbas,Dpq,Nbas,
     +            Dnull,Htau_pq,Nbas)

      Call Dgemm("T","T", Nbas,Nbas,Nbas,Onem,Hpq,Nbas,Dpq,Nbas,
     +            Dnull,Htau_qp,Nbas)

#ifdef _DEBUG_LVL0
      write(6,*)
      write(6,"(a)") "F(p,q) =  F(p,r)D(r,q) contribution" 
      call output(Htau_pq,1,Nbas,1,Nbas,Nbas,Nbas,1)
      write(6,*)
      write(6,"(a)") "F(q,p) =  F(r,q)D(p,r) contribution" 
      call output(Htau_qp,1,Nbas,1,Nbas,Nbas,Nbas,1)
      write(6,*)
#endif 

      Call Dscal(Nbas*Nbas,Two,Htau_pq,1)
      Call Transp(Work(I000),Work(I010),Nbas,Nbas)
      Call Daxpy(Nbas*Nbas,One,Work(I010),1,Htau_pq,1)

C This is -H(q,r)D(r,p) + H(r,p)D(q,r), which is from q+p, which
C is simply identical to symmetrizing (1/2(pq+qp)). 

      Call Dgemm("N","N", Nbas,Nbas,Nbas,Onem,Hpq,Nbas,Dpq,Nbas,
     +            Dnull,Work(I000),Nbas)
      Call Transp(Work(I000),Work(I010),Nbas,Nbas)

      Call Dgemm("T","T", Nbas,Nbas,Nbas,One,Hpq,Nbas,Dpq,Nbas,
     +            Dnull,Htau_pq,Nbas)

      Call Daxpy(Nbas*Nbas,One,Work(I010),1,Htau_pq,1)

      write(6,"(a)") "Htau_1d(p,q)-2"
      Call output(Htau_pq,1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif 

      Return
      End
   
