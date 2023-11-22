      Subroutine Pccd_form_dpq(Dcpq,Dspq,Doo,Dvv,Dvo,Work,Maxcor,
     +                         Nocc,Nvrt,Nbas,IUhf,Non_hf)

      Implicit Double Precision (A-H,O-Z)
      Logical pCCD,CCD,LCCD
      LOgical Non_hf

      Dimension Dcpq(Nbas,Nbas)
      Dimension Dspq(Nbas,Nbas)
      Dimension Dvv(Nvrt,Nvrt) 
      Dimension Doo(Nocc,Nocc) 
      Dimension Dvo(Nvrt,Nocc) 
      Dimension Work(Maxcor)

#include "machsp.com" 

      Common /Calc/pCCD,CCD,LCCD

      Data Ione,Inull,One,Two /1,0,1.0D0,2.0D0/

      Call Pccd_frmful_ext(Doo,Dvv,Dvo,Dcpq,Work,Maxcor,Nbas,0,
     +                     Non_hf)

      Scfocc = One 
      If (Iuhf .Eq. 1) Scfocc = One
      Call Dzero(Dspq,Nbas*Nbas)
      Do I = 1, Nocc
         Dspq(I,I) = Scfocc
      Enddo 

C Let's work with the alpha block of the Fock and the density matrices.
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
#ifdef _DEBUG_LVL0
      If (pCCD) Write(6,"(a)") "The pCCD density matrix"
      If (CCD)  Write(6,"(a)") "The CCD density matrix"
      If (LCCD) Write(6,"(a)") "The LCCD density matrix"
      call output(Dcpq,1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif
#endif 

#ifdef _HF_ONLY
      Scfocc = Two
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
      
