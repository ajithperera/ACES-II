      Subroutine Get_oneh(Oneh1,Oneh2,Oneh_aa,Oneh_bb,Scr1,
     +                    Scr2,E1_aa,E1_bb,Nbasis)

      Implicit Double Precision(A-H,O-Z)

      Dimension Oneh1(Nbasis*(Nbasis+1)/2)
      Dimension Oneh2(Nbasis,Nbasis) 
      Dimension Oneh_aa(Nbasis,Nbasis)
      Dimension Oneh_bb(Nbasis,Nbasis)
 
      Dimension Scr1(Nbasis*Nbasis),Scr2(Nbasis*Nbasis)

#include "cse12_vars.com"

      Ldim = Nbasis*(Nbasis+1)/2
      Call Getrec(20,"JOBARC","ONEHAO  ",Ldim,Oneh1)
      Call Expnd2(Oneh1,Oneh2,Nbasis)

#ifdef _DEBUG_LVLM
      Write(6,"(a)") "The AO basis one electron operator"
      call output(Oneh2,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
#endif

      Nbasis2 = Nbasis*Nbasis 
      If (Iend .GT. Maxcor) Call Insmem("prep_fock",Iend,Maxcor)

      Call Form_mo_fock(Oneh2,Oneh_aa,Scr1,Scr2,Nbasis,Nbasis,1)
      Call Form_mo_fock(Oneh2,Oneh_bb,Scr1,Scr2,Nbasis,Nbasis,2)

      E1_aa = 0.0D0
      DO I = 1, Nbasis 
         C = Ocn_oa(I)
         E1_aa = E1_aa + Oneh_aa(I,I)*C
      Enddo

      E1_bb = 0.0D0
      DO I = 1, Nbasis 
         C = Ocn_ob(i)
         E1_bb = E1_bb + oneh_bb(i,i)*C
      Enddo

#ifdef _DEBUG_LVLM
      Write(6,"(a,2(1x,F15.10))") "One electron energies: ", 
     +                             E1_aa,E1_bb
      Write(6,"(a)") "Alpha/beta one electron Hamiltonian"
      call output(oneh_aa,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
      call output(oneh_bb,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
#endif 

      Return
      End
