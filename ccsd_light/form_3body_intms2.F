      Subroutine form_3body_intms2(Qij_aa,Qij_bb,Qef_aa,Qef_bb,
     +                             Hbar_ijab1,Hbar_ijab2,Hbar_ijab3,
     +                             T2aa,T2bb,T2ab,Nocc_a,Nocc_b,
     +                             Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension Hbar_ijab1(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension Hbar_ijab2(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension Hbar_ijab3(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension Qij_aa(Nocc_a,Nocc_a)
      Dimension Qij_bb(Nocc_b,Nocc_b)
      Dimension Qef_aa(Nvrt_a,Nvrt_a)
      Dimension Qef_bb(Nvrt_a,Nvrt_a)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Call form_qef2(Qef_aa,Qef_bb,Hbar_ijab1,Hbar_ijab2,Hbar_ijab3,
     +               T2aa,T2bb,T2ab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Call form_qij2(Qij_aa,Qij_bb,Hbar_ijab1,Hbar_ijab2,Hbar_ijab3,
     +               T2aa,T2bb,T2ab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Return
      End 
  
     

