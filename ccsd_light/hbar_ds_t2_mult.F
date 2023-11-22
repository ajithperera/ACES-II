      Subroutine Hbar_ds_t2_mult(Hbar_iajk1,Hbar_iajk2,Hbar_iajk3,
     +                           Hbar_iajk4,Hbar_abci1,Hbar_abci2,
     +                           Hbar_abci3,Hbar_abci4,T1aa,T1bb,
     +                           T2aa,T2bb,T2ab,T2old_aa,T2old_bb,
     +                           T2old_ab,Ht2_abij1,Ht2_abij2,
     +                           Ht2_abij3,Qij_aa,Qij_bb,Qef_aa,
     +                           Qef_bb,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension Hbar_iajk1(Nocc_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Hbar_iajk2(Nocc_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Hbar_iajk3(Nocc_a,Nvrt_b,Nocc_a,Nocc_b)
      Dimension Hbar_iajk4(Nocc_b,Nvrt_a,Nocc_b,Nocc_a)

      Dimension Hbar_abci1(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Hbar_abci2(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Hbar_abci3(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Hbar_abci4(Nvrt_b,Nvrt_a,Nvrt_b,Nocc_a)
    
      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension T2old_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2old_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2old_ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension Qij_aa(Nocc_a,Nocc_a)
      Dimension Qij_bb(Nocc_b,Nocc_b)
      Dimension Qef_aa(Nvrt_a,Nvrt_a)
      Dimension Qef_bb(Nvrt_b,Nvrt_b)

      Dimension Ht2_abij1(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Ht2_abij2(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Ht2_abij3(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

C The T1_int2 need to be modified to accomodate the new storage pattern
C Hbar(ijka) (W ints were stored as W4(ij,ka)).

      Call T1_int2_modf(T1aa,T1bb,Ht2_abij1,Ht2_abij2,
     +                  Ht2_abij3,Hbar_abci1,Hbar_abci2,
     +                  Hbar_abci3,Hbar_abci4,Hbar_iajk1,
     +                  Hbar_iajk2,Hbar_iajk3,Hbar_iajk4,
     +                  Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Call Fae_int2(T2old_aa,T2old_bb,T2old_ab,Ht2_abij1,
     +              Ht2_abij2,Ht2_abij3,Qef_aa,Qef_bb,
     +              Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Call Fmi_int2(T2old_aa,T2old_bb,T2old_ab,Ht2_abij1,
     +              Ht2_abij2,Ht2_abij3,Qij_aa,Qij_bb,
     +              Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Return
      End 
    
   
