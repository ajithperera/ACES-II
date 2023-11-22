      Subroutine form_3body_intms1(Qij_aa,Qij_bb,Qef_aa,Qef_bb,
     +                             Hbar_aibc1,Hbar_aibc2,Hbar_aibc3,
     +                             Hbar_aibc4,Hbar_ijka1,Hbar_ijka2,
     +                             Hbar_ijka3,Hbar_ijka4,T1aa,T1bb,
     +                             Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension Hbar_aibc1(Nvrt_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension Hbar_aibc2(Nvrt_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension Hbar_aibc3(Nvrt_a,Nocc_b,Nvrt_a,Nvrt_b)
      Dimension Hbar_aibc4(Nvrt_b,Nocc_a,Nvrt_b,Nvrt_a)

      Dimension Hbar_ijka1(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension Hbar_ijka2(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension Hbar_ijka3(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension Hbar_ijka4(Nocc_b,Nocc_a,Nocc_b,Nvrt_a)

      Dimension Qij_aa(Nocc_a,Nocc_a)
      Dimension Qij_bb(Nocc_b,Nocc_b)
      Dimension Qef_aa(Nvrt_a,Nvrt_a)
      Dimension Qef_bb(Nvrt_a,Nvrt_a)
    
      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Call form_qef1(Qef_aa,Qef_bb,Hbar_aibc1,Hbar_aibc2,Hbar_aibc3,
     +               Hbar_aibc4,T1aa,T1bb,Nocc_a,Nocc_b,Nvrt_a,
     +               Nvrt_b)

      Call form_qij1(Qij_aa,Qij_bb,Hbar_ijka1,Hbar_ijka2,Hbar_ijka3,
     +               Hbar_ijka4,T1aa,T1bb,Nocc_a,Nocc_b,Nvrt_a,
     +               Nvrt_b)

      Return
      End 
  
     

