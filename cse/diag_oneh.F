      Subroutine Diag_oneh(Oneh_aa,Oneh_bb,Dens_aa,Dens_bb,Repuls,
     +                     Nbasis,Naobfns,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,
     +                     Work,Memleft)

      Implicit Double Precision (A-H, O-Z)

      Dimension Work(Memleft)
      Dimension Oneh_aa(Nbasis,Nbasis),Oneh_bb(Nbasis,Nbasis)
      Dimension Dens_aa(Nbasis,Nbasis),Dens_bb(Nbasis,Nbasis)
      Dimension Nocca(8),Noccb(8)

      Data Ione /1/

#include "cse_vars.com"

      Nbasis2 = Nbasis*Nbasis

      Ievalr = Ione 
      Ievall = Ievalr + Nbasis
      Ievecr = Ievall + Nbasis
      Ievecl = Ievecr + Nbasis2
      Iscr1  = Ievecl + Nbasis2
      Iscr2  = Iscr1  + 4*Nbasis 
      Ifoccn = Iscr2  + Nbasis2
      Iend   = Ifoccn + Nbasis2

      If (Iend .Ge. Memleft) Call Insmem("Diag_oneh",Iend,
     +                                   Memleft)
 
C Diagonalize the AA block and form the AA right density.

      Call Dgeev("V","V",Nbasis,Oneh_aa,Nbasis,Work(Ievalr),
     +            Work(Ievall),Work(Ievecl),Nbasis,Work(Ievecr),
     +            Nbasis,Work(Iscr1),4*Nbasis,Ierr)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Left Alpha eigenvectors"
      call output(Work(Ievecl),1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
      Write(6,"(a)") "Right Alpha eigenvectors"
      call output(Work(Ievecr),1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
#endif 

      Call Dzero(Work(Ifoccn),Nbasis2) 
      Call Getrec(20,"JOBARC",'NIRREP  ',Ione  ,Nirrep)
      Call Getrec(20,"JOBARC",'OCCUPYA0',Nirrep,Nocca)
      Call Getrec(20,"JOBARC",'OCCUPYB0',Nirrep,Noccb)

      Nalpha = 0
      Nbeta  = 0
      Do I = 1, Nirrep
         Nalpha = Nalpha + Nocca(I)
         Nbeta  = Nbeta  + Noccb(I)
      Enddo 
      
      Do I = 1, Nalpha 
         Ioff =  (I-1)*(Nbasis + 1)
         Work(Ifoccn + Ioff) = 1.0D0
      Enddo 

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Canonical alpha  occupation matrix"
      call output(Work(Ifoccn),1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
#endif 

      Call Reord_eigs(Work(Ievalr),Nbasis)

      Write(6,*)
      Write(6,"(a)") "The alpha eigenvalues of one hamiltonian"
      Write(6,"(5(1x,F15.10))") (Work(Ievalr+I-1),I=1,Nbasis)

      E = 0.0D0
      Do I = 1, Nbasis
         E = E + Work(Ievalr+I-1)*Work(Ifoccn+I-1)
         Oneh_aa(I,I) = Work(Ievalr+I-1) 
      Enddo 

      Call Form_dens(Dens_aa,Work(Ievecr),Work(Iscr2),Work(Ifoccn),
     +               Nbasis) 

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Alpha Density  matrix"
      call output(Dens_aa,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
#endif 
 
C Diagonalize the BB block and form the BB right density.

      Call Dgeev("V","V",Nbasis,Oneh_bb,Nbasis,Work(Ievalr),
     +            Work(Ievall),Work(Ievecl),Nbasis,Work(Ievecr),
     +            Nbasis,Work(Iscr1),4*Nbasis,Ierr)

      Do I = 1, Nbeta 
         Ioff =  (I-1)*(Nbasis + 1)
         Work(Ifoccn + Ioff) = 1.0D0
      Enddo 

      Call Reord_eigs(Work(Ievalr),Nbasis)
      Write(6,*)
      Write(6,"(a)") "The beta eigenvalues of one hamiltonian"
      Write(6,"(5(1x,F15.10))") (Work(Ievalr+I-1),I=1,Nbasis)
      Write(6,*)

      Do I = 1, Nbasis
         E = E + Work(Ievalr+I-1)*Work(Ifoccn+I-1)
         Oneh_bb(I,I) = Work(Ievalr+I-1)
      Enddo

      Write(6,"(2a)") "The reference energy printed here is",
     +                " for canonical occupations (i.e. not read"
      Write(6,"(a)")  "from occupation number arrays)."
      Write(6,"(a,(1X,F15.10))") "The reference energy :",E+Repuls
      Write(6,*)

      Call Form_dens(Dens_bb,Work(Ievecr),Work(Iscr2),Work(Ifoccn),
     +               Nbasis) 

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Beta Density  matrix"
      call output(Dens_aa,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
#endif 

      Return
      End
     
