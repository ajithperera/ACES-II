      Subroutine Pccd_form_htau(Htau_qp,Htau_pq,Hoo_qp,Hoo_pq,Hvv_qp,
     +                          Hvv_pq,Hvo,Hov,Work,Maxcor,Nocc,Nvrt,
     +                          Nbas)

      Implicit Double Precision(A-H,O-Z)
      Logical Sym_packed 
      Logical Symmetry
      Logical Opt_orbs,OO_constr,OV_constr

      Dimension Htau_pq(Nbas,Nbas)
      Dimension Htau_qp(Nbas,Nbas)
      Dimension Hoo_pq(Nocc*Nocc)
      Dimension Hoo_qp(Nocc*Nocc)
      Dimension Hvv_pq(Nvrt*Nvrt)
      Dimension Hvv_qp(Nvrt*Nvrt)
      Dimension Hvo(Nvrt*Nocc)
      Dimension Hov(Nocc*Nvrt)
      Dimension Work(Maxcor)

      Data Ione,Onem,One,Dnull,Half,Two/1,-1.0D0,1.0D0,0.0D0,0.50D0,
     +                                  2.0D0/

#include "syminf.com"
#include "sym.com"
#include "machsp.com"
#include "sympop.com"

      Common /Symm/Symmetry 
      Common /OO_info/Opt_orbs,OO_constr,OV_constr 

      Nmax   = Max(Nvrt,Nocc)
      Irrepx = Ione 

      If (Symmetry) Then
         Ndim_oo = Irpdpd(Irrepx,21)
         Ndim_vv = Irpdpd(Irrepx,19)
         Ndim_vo = Irpdpd(Irrepx,9)
      Else
         Ndim_oo   = Nocc*Nocc
         Ndim_vv   = Nvrt*Nvrt
         Ndim_vo   = Nocc*Nvrt
      Endif

      I000 = Ione
      I010 = I000 + Nmax*Nmax
      I020 = I010 + Nbas*Nbas
      I030 = I020 + Nmax*Nmax
      Iend = I030 + Nbas*Nbas
      Maxcor = Maxcor - Iend 
      If (Iend.Gt.Maxcor) Call Insmem("Pccd_form_htau",Iend,Maxcor)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " The Htau_pq (from one-electron Hamiltonian)"
      Call output(Htau_pq,1,Nbas,1,Nbas,Nbas,Nbas,1)
      Write(6,"(a)") " The Htau_qp(from one-electron Hamiltonian)"
      Call output(Htau_qp,1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif 
C-----------------Energy as trace of Gfock for consistency check-------------

      Sym_packed = .True.
      Call Pccd_frmful(Work(I010),Hoo_qp,Hvv_qp,Hov,Hvo,Work(Iend),
     +                 Maxcor,Nocc,Nvrt,Nbas,"Ov_like",Sym_packed)

#ifndef _HF_ONLY
C Scale the diagonals by half (correlation contribution only).

      Call Dscal(Nbas,Half,Work(I010),Nbas+1)
#endif
      Call Daxpy(Nbas*Nbas,One,Work(I010),1,Htau_qp,1)

      Call Getrec(20, 'JOBARC', 'NUCREP', Iintfp, Zrepl)

      e=dnull
#ifndef _HF_ONLY
      do i = 1, Nbas
         e = e + Htau_qp(i,i)
      enddo
#else 
      do i = 1, Nbas
         e = e + Htau_qp(i,i)
      enddo
#endif 

#ifdef _DEBUG_LVL0
      write(6,*)
      write(6,"(a,2(2x,F15.10))") " The electronic energy :", e
#endif

      write(6,*)

#ifdef _HF_ONLY
      e = e + Zrepl
      write(6,"(a,2(2x,F15.10))") " The reference energy  :", e
      Call DScal(Nbas,Dnull,Htau_qp,Nbas+1)
      Call Pccd_antisymmetrize(Htau_qp,Work,Maxcor,Nbas,Nocc,Nvrt)
#ifdef _DEBUG_LVL0
      Write(6,"(a)") " The HF orbital rotation gradient"
      Call output(Htau_qp,1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif 
      Call Dcopy(Nbas*Nbas,Htau_qp,1,Htau_pq,1)
      Return
#else
#ifdef _DEBUG_LVL0
      e = e
      write(6,"(a,(2x,F15.10))") " The correlation energy:", e
      write(6,*)
      write(6,"(2a)") " Warning! The correlation energy computed",
     +                " as the trace of the extended Fock matrix"
      Write(6,"(a)")  " is printed for verification purposes only"
#endif 
#endif
#ifdef _NO_NRML_ORD
      e = e + Zrepl
      write(6,"(a,2(2x,F15.10))") " The total energy:", e
#endif
      write(6,*)
C--------Energy as trace of Gfock for consistency check is complete-------------

      Call Getrec(20,"JOBARC","REFGRDOO",Nfmi(1)*Iintfp,Work(I000))
      Call Dscal(Ndim_oo,Onem,Hoo_qp,1)
      Call Dscal(Ndim_vv,Onem,Hvv_qp,1)
      Call Dscal(Ndim_vo,Onem,Hov,1)
      Call Dscal(Mdim_vo,Onem,Hvo,1)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(2a)") " The OO=(1/2) D(pq)[<pi||qj>",
     +               "+<pj||qi>]) term"
      Write(6,"(5(1x,F15.9))") (Work(I000+i-1), i=1,Nfmi(1))
      call checksum("Hoo_pq :",Work(I000),Nfmi(1))
#endif 

C This contribution is zero for orbital rotations. It is safer to add them
C and let the antisymmetrization set them to zero. This way we can have 
C exact identities to ACES II analytic gradient code.

      Call Daxpy(Ndim_oo,One,Work(I000),1,Hoo_qp,1)

      Call Getrec(20,"JOBARC","REFGRDOV",Nt(1)*Iintfp,Work(I000))
      Call Daxpy(Ndim_vo,One,Work(I000),1,Hvo,1)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(2a)") " The VO=[<ip||aq> + <iq//ap>]D(p,q)",
     +                " term"
      Write(6,"(5(1x,F15.9))") (Work(I000+i-1), i=1,Nt(1))
      call checksum("Hvo_pq :",Work(I000),Nt(1))
#endif 

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") " The Hoo_qp"
      Write(6,"(5(1x,F15.9))") (Hoo_qp(i),i=1,Ndim_oo)
      Write(6,"(a)") " The Hvv_qp"
      Write(6,"(5(1x,F15.9))") (Hvv_qp(i),i=1,Ndim_vv)
      Write(6,"(a)") " The VO and OV blocks"
      Write(6,"(5(1x,F15.9))") (Hvo(i),i=1,Ndim_vo)
      Write(6,*)
      Write(6,"(5(1x,F15.9))") (Hov(i),i=1,Ndim_vo)
#endif 
      Call Pccd_frmful(Work(I010),Hoo_qp,Hvv_qp,Hvo,Hov,Work(Iend),
     +                 Maxcor,Nocc,Nvrt,Nbas,"Ov_like",Sym_packed) 
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "-------Htau_pq construction------"
      Write(6,"(a)") "The correlation + OO contribution"
      Call output(Work(I010),1,Nbas,1,Nbas,Nbas,Nbas,1)
CSSS      c = Dnull
CSSS      do i = 1, Nbas
CSSS         Ioff = (i-1)*nbas + (i-1)
CSSS         c = c + Work(I010+Ioff)
CSSS      enddo 
CSSS      write(6,"(a,2(2x,F15.10))") " The test energy :", c
      write(6,*)
      Write(6,"(a)") "The Htau_pq=D(p,r)*F(r,q) contribution"
      Call output(Htau_pq,1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif 

C *Warning!!--- Change sign here break first step
      Call Daxpy(Nbas*Nbas,One,Work(I010),1,Htau_pq,1)

#ifdef _DEBUG_LVL0
      Write(6,"(a,a)") " The Htau_qp (after OV/VO blocks added",
     +                 " and not anti-symmetrized)"
      Call output(Htau_pq,1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif 

C Antisymmetrize to generate the orbital rotation gradients matrix.The scalling
C by minus two accomplish two things. The minus is needed to get the correct
C sign for the VO and OV blocks (VO block must be positive). The two gives the
C orbital rotaion gradient of the Alpha block (and match with half of the value
C obtained numerically). Note that proceeding steps need to be tailored to the
C fact that we are working with alpha block only.

      Call DScal(Nbas,Dnull,Htau_pq,Nbas+1)

      Call Pccd_antisymmetrize(Htau_pq,Work,Maxcor,Nbas,Nocc,Nvrt)

C This scalling by two is not necessary. For RHF we can work with Alpha
C or beta blocks. See also pccd_scale_grads for the factor 1/2. 

CSSS      Call Dscal(Nbas*Nbas,Onem,Htau_pq,1)

      Call Pccd_htau_trnspose(Htau_pq,Work(Iend),Nbas,Nocc,Nvrt)

      If (OO_constr) Call pccd_insist_no_oovv(Htau_pq,Nocc,Nvrt,Nbas)
      If (OV_constr) Call pccd_insist_no_ovvo(Htau_pq,Nocc,Nvrt,Nbas)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " The orbital rotation gradient"
      Call output(Htau_pq,1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif 

      Return
      End


