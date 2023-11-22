      Subroutine Psi4dbg_form_dpq(Dcpq,Dspq,Doo,Dvv,Dvo,Work,Maxcor,
     +                            Nocc,Nvrt,Nbas,Iuhf,Non_hf)

      Implicit Double Precision (A-H,O-Z)
      Logical pCCD,CCD,LCCD
      Logical Non_hf

      Dimension Dcpq(Nbas,Nbas)
      Dimension Dspq(Nbas,Nbas)
      Dimension Dvv(Nvrt,Nvrt) 
      Dimension Doo(Nocc,Nocc) 
      Dimension Dvo(Nvrt,Nocc) 
      Dimension Work(Maxcor)

#include "machsp.com" 
      Common /Meth/pCCD,CCD,LCCD

      Data Ione,Inull,One,Two /1,0,1.0D0,2.0D0/

      Call Psi4dbg_frmful(Doo,Dvv,Dvo,Dcpq,Work,Maxcor,Nbas,0,
     +                    Non_hf)

      Scfocc = One
CSSS      Call Dscal(Nbas*Nbas,Two,Dcpq,1)

C This density matrix is post-HF ordered (i.e. all occupied for all irreps
C fist then followed by virtuals in all irreps).

#ifdef _NO_NRML_ORD
#ifndef _HF_ONLY
      Do I = 1, Nocc
         Dcpq(I,I) = Dcpq(I,I) + Scfocc
      Enddo 
#ifdef _DEBUG_LVL0
      If (pCCD) Write(6,"(a)") "The pCCD + SCF density matrix"
      If (CCD) Write(6,"(a)")  "The CCD + SCF density matrix"
      If (LCCD) Write(6,"(a)") "The LCCD + SCF density matrix"
      call output(Dcpq,1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif 
#endif 
#else
      If (pCCD) Write(6,"(a)") "The pCCD density matrix"
      If (CCD)  Write(6,"(a)") "The CCD density matrix"
      If (LCCD) Write(6,"(a)") "The LCCD density matrix"
      call output(Dcpq,1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif 

#ifdef _HF_ONLY
      Scfocc = One
      Call Dzero(Dcpq,Nbas*Nbas)
      Do I = 1, Nocc
         Dcpq(I,I) = Scfocc
      Enddo 
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The HF density matrix"
      call output(Dcpq,1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif 
#endif 

      Return
      End
      
