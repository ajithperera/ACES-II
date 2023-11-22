      Subroutine Pccd_nodiis(Grd,Work,Maxcor,Ispin,Nbas,Ncycle)

      Implicit Double Precision(A-H,O-Z)
    
      Dimension Work(Maxcor)
      Dimension Grd(Nbas,Nbas)
#include "machsp.com"
#include "sym.com"
#include "syminf.com"
#include "sympop.com"
      Common /Symm/Symmetry

      Data Ione,Onem,half,Dnull,One/1,-1.0D0,0.50D0,0.0D0,1.0D0/

C Form  K-K(^t). This appears double dipping but help convergence
C tremendoulsy.
      Nbas2 = Nbas*Nbas
      I000 = Ione
      I010 = I000 + Nbas2
      I020 = I010 + Nbas2
      I030 = I020 + Nbas2
      I040 = I030 + Nbas2
      I050 = I040 + Nbas2
      
      If (Iend .Gt. Maxcor) Call Insmem("pccd_nodiis",Iend,Maxcor)
      Maxcor = Maxcor-Iend
    
CSSS      Call Transp(Grd,Work(I000),Nbas,Nbas)
CSSS      Call Daxpy(Nbas*Nbas,Onem,Work(I000),1,Grd,1)

#ifdef _DEBUG_LVL0
      write(6,*)
      Write(6,"(a)") "Antisymmetrized Scalled gradients (K-K(^t))"
      Call output(Grd,1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif

C 1/2 K*K(^t)
      Call Dgemm("N","T",Nbas,Nbas,Nbas,Half,Grd,Nbas,Grd,Nbas,
     +            Dnull,Work(I000),Nbas)

      Do I = 1, Nbas
         Grd(I,I) = Grd(I,I) + One
      Enddo

      Call Daxpy(Nbas*Nbas,One,Work(I000),1,Grd,1)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "U=(1+G+1/2G*K)"
      Call output(Grd,1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif

      Call Pccd_gramschmidt(Grd,Nbas,Nbas)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Unitary check of Kappa"
      Call Dgemm("N","T",Nbas,Nbas,Nbas,One,Grd,Nbas,Grd,Nbas,
     +           Dnull,Work(I000),Nbas)
      Call output(Work(i000),1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif

      If (Ispin .Eq. 1) Call Getrec(20,"JOBARC","SCFEVECA",Nbas2*Iintfp,
     +                              Work(I010))
      If (Ispin .Eq. 2) Call Getrec(20,"JOBARC","SCFEVECB",Nbas2*Iintfp,
     +                              Work(I010))

#ifdef _BRUECK_LIKE
C Form new C(new), where C(new) =  G (C(old))

      Call Dgemm("N","N",Nbas,Nbas,Nbas,One,Grd,Nbas,Work(I010),Nbas,
     +            Dnull,Work(I000),Nbas)

C Orthogonalize C(New)

      Call Getrec(20,"JOBARC","AOOVRLAP",Nbas2*Iintfp,Work(I010))
      Call Ao2mo2(Work(I010),Work(I010),Work(I020),Work(I030),Nbas,
     +           Nbas,Ispin)
      Call Eig(Work(I010),Work(I020),Nbas,Nbas,0)

      Call Invsqt(Work(I010),Nbas+1,Nbas)
      Call Dgemm("N","N",Nbas,Nbas,Nbas,One,Work(I020),Nbas,
     +            Work(I010),Nbas,Dnull,Work(I030),Nbas)
      Call Dgemm("N","N",Nbas,Nbas,Nbas,One,Work(I030),Nbas,
     +            Work(I020),Nbas,Dnull,Work(I010),Nbas)

      Call Dgemm("N","N",Nbas,Nbas,Nbas,One,Work(I000),Nbas,
     +            Work(I010),Nbas,Dnull,Work(I020),Nbas)

      If (Ispin .Eq. 1) Call Putrec(20,"JOBARC","SCFEVCA0",Nbas2*Iintfp,
     +                              Work(I020))
      If (Ispin .Eq. 2) Call Putrec(20,"JOBARC","SCFEVCB0",Nbas2*Iintfp,
     +                              Work(I020))
#else

      If (Ispin .Eq. 1) Call Getrec(20,"JOBARC","EVECOAOA",Nbas2*Iintfp,
     +                              Work(I010))
      If (Ispin .Eq. 2) Call Getrec(20,"JOBARC","EVECOAOB",Nbas2*Iintfp,
     +                              Work(I010))

      If (Ispin .Eq. 1) Call Putrec(20,"JOBARC","EVCOAOXA",Nbas2*Iintfp,
     +                              Work(I010))
      If (Ispin .Eq. 2) Call Putrec(20,"JOBARC","EVCOAOXB",Nbas2*Iintfp,
     +                              Work(I010))

      Call Mo2oao(Grd,Work(I000),Work(I010),Work(I020),Nbas,Ispin)
      call output(Work(i000),1,Nbas,1,Nbas,Nbas,Nbas,1)

      Call Dgemm("N","N",Nbas,Nbas,Nbas,One,Work(I000),Nbas,
     +            Work(I010),Nbas,Dnull,Work(I020),Nbas)

      If (Ispin .Eq. 1) Call Putrec(20,"JOBARC","SCFEVCA0",Nbas2*Iintfp,
     +                              Work(I020))
      If (Ispin .Eq. 2) Call Putrec(20,"JOBARC","SCFEVCB0",Nbas2*Iintfp,
     +                              Work(I020))
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Orthogonalized new rotated Vectos in OAO basis"
      call output(Work(i020),1,Nbas,1,Nbas,Nbas,Nbas,1)
      call Dgemm("T","N",Nbas,Nbas,Nbas,One,Work(I020),Nbas,Work(I020),
     +            Nbas,Dnull,Work(I030),Nbas)
      Ones = Ddot(Nbas*Nbas,Work(I030),1,Work(I030),1)
      If (Int(Ones) .Eq. Nbas) Write(6,"(2a)") "Normalization test",
     +                                         " passes" 
#endif 
#endif
CSSS      Call Pccd_rotate(Grd,Work(Iend),Maxcor,Ispin,Nbas,Ncycle)

      Return
      End
