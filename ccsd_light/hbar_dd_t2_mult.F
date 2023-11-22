      Subroutine Hbar_dd_t2_mult(Hbar_aea,Hbar_aeb,Hbar_mia,Hbar_mib,
     +                           Hbar_abef1,Hbar_abef2,Hbar_abef3,
     +                           Hbar_ijkl1,Hbar_ijkl2,Hbar_ijkl3,
     +                           Hbar_mbej1,Hbar_mbej2,Hbar_mbej3,
     +                           Hbar_mbej4,Hbar_mbej5,Hbar_mbej6,
     +                           Ht2_abij1,Ht2_abij2,Ht2_abij3,
     +                           Qij_aa,Qij_bb,Qef_aa,Qef_bb,T2aa,
     +                           T2bb,T2ab,T2old_aa,T2old_bb,T2old_ab,
     +                           Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension Hbar_abef1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension Hbar_abef2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension Hbar_abef3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Dimension Hbar_ijkl1(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension Hbar_ijkl2(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension Hbar_ijkl3(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Dimension Hbar_mbej1(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Hbar_mbej2(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Hbar_mbej3(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Hbar_mbej4(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension Hbar_mbej5(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension Hbar_mbej6(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Dimension Ht2_abij1(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Ht2_abij2(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Ht2_abij3(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension Hbar_mia(Nocc_a,Nocc_a)
      Dimension Hbar_mib(Nocc_b,Nocc_b)
      Dimension Hbar_aea(Nvrt_a,Nvrt_a)
      Dimension Hbar_aeb(Nvrt_b,Nvrt_b)

      Dimension Qij_aa(Nocc_a,Nocc_a)
      Dimension Qij_bb(Nocc_b,Nocc_b)
      Dimension Qef_aa(Nvrt_a,Nvrt_a)
      Dimension Qef_bb(Nvrt_b,Nvrt_b)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension T2old_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2old_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2old_ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Data Onem/-1.0D0/

      Call Fae_int2(T2aa,T2bb,T2ab,Ht2_abij1,Ht2_abij2,
     +              Ht2_abij3,Hbar_aea,Hbar_aeb,Nocc_a,
     +              Nocc_b,Nvrt_a,Nvrt_b)

      Call Fmi_int2(T2aa,T2bb,T2ab,Ht2_abij1,Ht2_abij2,
     +              Ht2_abij3,Hbar_mia,Hbar_mib,Nocc_a,
     +              Nocc_b,Nvrt_a,Nvrt_b)

      Call Wmnij_int2(T2aa,T2bb,T2ab,Ht2_abij1,Ht2_abij2,
     +                Ht2_abij3,Hbar_ijkl1,Hbar_ijkl2,
     +                Hbar_ijkl3,Nocc_a,Nocc_b,Nvrt_a,
     +                Nvrt_b)

      Call Wabef_int2(T2aa,T2bb,T2ab,Ht2_abij1,Ht2_abij2,
     +                Ht2_abij3,Hbar_abef1,Hbar_abef2,
     +                Hbar_abef3,Nocc_a,Nocc_b,Nvrt_a,
     +                Nvrt_b)

      Call Wmbej_int2(T2aa,T2bb,T2ab,Ht2_abij1,Ht2_abij2,
     +                Ht2_abij3,Hbar_mbej1,Hbar_mbej2,
     +                Hbar_mbej3,Hbar_mbej4,Hbar_mbej5,
     +                Hbar_mbej6,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Call Fae_int2(T2old_aa,T2old_bb,T2old_ab,Ht2_abij1,
     +              Ht2_abij2,Ht2_abij3,Qef_aa,Qef_bb,Nocc_a,
     +              Nocc_b,Nvrt_a,Nvrt_b)

      Call Fmi_int2(T2old_aa,T2old_bb,T2old_ab,Ht2_abij1,
     +              Ht2_abij2,Ht2_abij3, Qij_aa,Qij_bb,Nocc_a,
     +              Nocc_b,Nvrt_a,Nvrt_b)

      Return
      End 
    
   
