      Subroutine Form_cumldens(T1aa,T1bb,Fockpq_aa,Fockpq_bb,Work,
     +                         Maxcor,Nbasis,Ref_energy,Nocc_a,
     +                         Nocc_b,Nvrt_a,Nvrt_b,Frac_occ)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
     
      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)
      Dimension Fockpq_aa(Nocc_a,Nocc_a)
      Dimension Fockpq_bb(Nocc_b,Nocc_b)

      Logical  Frac_occ

      Data Ione /1/
   
      Nbasis2   = Nbasis*Nbasis 
      Idensp_aa = Ione
      Idensp_bb = Idensp_aa + Nbasis2
      Idensm_aa = Idensp_bb + Nbasis2
      Idensm_bb = Idensm_aa + Nbasis2
      Iend      = Idensm_bb + Nbasis2

      If (Iend .GT. Maxcor) Call Insmem("form_cumldens",Iend,Maxcor)

      Call Dens_frmt1(Work(Idensp_aa),Work(Idensp_bb),Work(Idensm_aa),
     +                Work(Idensm_bb),T1aa,T1bb,Fockpq_aa,Fockpq_bb,
     +                Ref_energy,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,Nbasis,
     +                Frac_occ) 

C Lets get the natural orbitals and the eigenvalues. 

      Inatorb_r = Iend 
      Inatorb_l = Inatorb_r + Nbasis2
      Ieval_r   = Inatorb_l + Nbasis2
      Ieval_l   = Ieval_r   + Nbasis
      Iscr1     = Ieval_l   + Nbasis 
      Iscr2     = Iscr1     + Nbasis2
      Iend      = Iscr2     + 4*Nbasis 
     
      Call Dcopy(Nbasis2,Work(Idensp_aa),1,Work(Iscr1),1)
      If (Iend .GT. Maxcor) Call Insmem("form_cumldens",Iend,Maxcor)

      Call Dgeev("V","V",Nbasis,Work(Iscr1),Nbasis,Work(Ieval_r),
     +            Work(Ieval_l),Work(Inatorb_r),Nbasis,
     +            Work(Inatorb_l),Nbasis,Work(Iscr2),4*Nbasis,Ierr)


#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The right eigenvalues of the cumulant density"
      Write(6,"(6(1x,F15.10))") (Work(Ieval_r+I-1),i=1,Nbasis)
      Write(6,*)
#endif

      Return 
      End 
