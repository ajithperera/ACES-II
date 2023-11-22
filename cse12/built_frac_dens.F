      Subroutine Built_frac_dens(Dens_aa,Dens_bb,Dens_tt,Evec,
     +                          Docc,Scr,Nbasis)

      Implicit Double Precision(A-H,O-Z)

      Dimension Dens_aa(Nbasis,Nbasis) 
      Dimension Dens_bb(Nbasis,Nbasis) 
      Dimension Dens_tt(Nbasis,Nbasis)

      Dimension Docc(Nbasis,Nbasis) 
      Dimension Evec(Nbasis*Nbasis)
      Dimension Scr(Nbasis*Nbasis)

      Integer P,Q,R,S

#include "cse12_vars.com"
#include "machsp.com"

      Nbasis2 = Nbasis*Nbasis 
      Call Dzero(Docc,Nbasis2)

C built AA AO density matrix

      Do P = 1, Nbasis
         Docc(P,P) = Ocn_oa(P)
      Enddo

      Call Getrec(20,"JOBARC","SCFEVCA0",Nbasis*Nbasis*Iintfp,
     +            Evec)
      Call Form_dens(Dens_aa,Evec,Scr,Docc,Nbasis)

C built BB AO density matrix 

      Do P = 1, Nbasis 
         Docc(P,P) = Ocn_ob(P)
      Enddo

      Call Getrec(20,"JOBARC","SCFEVCB0",Nbasis*Nbasis*Iintfp,
     +               Evec)
      Call Form_dens(Dens_bb,Evec,Scr,Docc,Nbasis)
     
      Call Dcopy(Nbasis2,Dens_aa,1,Dens_tt,1)
      Call Daxpy(Nbasis2,1.0D0,Dens_bb,1,Dens_tt,1)
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Alpha/beta and total density"
      call output(Dens_aa,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
      call output(Dens_bb,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
      call output(Dens_tt,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
#endif

      Return
      End 


