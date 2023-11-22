      Subroutine Hbar_ss_t1_mult(Hbar_aea,Hbar_aeb,Hbar_mia,Hbar_mib,
     +                           Hbar_mbej1,Hbar_mbej2,Hbar_mbej3,
     +                           Hbar_mbej4,Hbar_mbej5,Hbar_mbej6,
     +                           T1aa,T1bb,Ht1_aiaa,Ht1_aibb,
     +                           Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension Hbar_mbej1(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Hbar_mbej2(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Hbar_mbej3(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Hbar_mbej4(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension Hbar_mbej5(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension Hbar_mbej6(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Dimension Hbar_mia(Nocc_a,Nocc_a)
      Dimension Hbar_mib(Nocc_b,Nocc_b)
      Dimension Hbar_aea(Nvrt_a,Nvrt_a)
      Dimension Hbar_aeb(Nvrt_b,Nvrt_b)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension Ht1_aiaa(Nvrt_a,Nocc_a)
      Dimension Ht1_aibb(Nvrt_b,Nocc_b)

      Call T1_int1_modf(T1aa,T1bb,Ht1_aiaa,Ht1_aibb,Hbar_mbej1,
     +                  Hbar_mbej2,Hbar_mbej3,Hbar_mbej4,
     +                  Hbar_mbej5,Hbar_mbej6,Nocc_a,Nocc_b,
     +                  Nvrt_a,Nvrt_b)

      Call Fae_int1(T1aa,T1bb,Ht1_aiaa,Ht1_aibb,Hbar_aea,Hbar_aeb,
     +              Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Call Fmi_int1(T1aa,T1bb,Ht1_aiaa,Ht1_aibb,Hbar_mia,Hbar_mib,
     +              Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

 
      Return
      End
