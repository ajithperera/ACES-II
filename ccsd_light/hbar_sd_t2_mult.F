      Subroutine Hbar_sd_t2_mult(Hbar_mea,Hbar_meb,Hbar_aibc1,
     +                           Hbar_aibc2,Hbar_aibc3,Hbar_aibc4,
     +                           Hbar_ijka1,Hbar_ijka2,Hbar_ijka3,
     +                           Hbar_ijka4,T2aa,T2bb,T2ab,Ht1_aiaa,
     +                           Ht1_aibb,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension Hbar_mea(Nocc_a,Nvrt_a)
      Dimension Hbar_meb(Nocc_b,Nvrt_b)

      Dimension Hbar_aibc1(Nvrt_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension Hbar_aibc2(Nvrt_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension Hbar_aibc3(Nvrt_a,Nocc_b,Nvrt_a,Nvrt_b)
      Dimension Hbar_aibc4(Nvrt_b,Nocc_a,Nvrt_b,Nvrt_a)

      Dimension Hbar_ijka1(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension Hbar_ijka2(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension Hbar_ijka3(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension Hbar_ijka4(Nocc_b,Nocc_a,Nocc_b,Nvrt_a)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)
      
      Dimension Ht1_aiaa(Nvrt_a,Nocc_a)
      Dimension Ht1_aibb(Nvrt_b,Nocc_b)

      Call Fme_int1(T2aa,T2bb,T2ab,Ht1_aiaa,Ht1_aibb,Hbar_mea,
     +              Hbar_meb,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C The T2_int1 need to be modified to accomodate the new storage pattern
C Hbar(ai,bc) (W4 ints were stored as W4(ab,ci)). 

      Call T2_int1_modf(T2aa,T2bb,T2ab,Hbar_aibc1,Hbar_aibc2,
     +                  Hbar_aibc3,Hbar_aibc4,Hbar_ijka1,
     +                  Hbar_ijka2,Hbar_ijka3,Hbar_ijka4,
     +                  Ht1_aiaa,Ht1_aibb,Nocc_a,Nocc_b,Nvrt_a,
     +                  Nvrt_b)

      Return
      End 
