      Subroutine Form_t2intms(T2old_aa,T2old_bb,T2old_ab,W_aa,
     +                        W_bb,W_ab,Fmi_a,Fmi_b,Fae_a,Fae_b,
     +                        Wmnij_1,Wmnij_2,Wmnij_3,Wabef_1,
     +                        Wabef_2,Wabef_3,Wmbej_1,Wmbej_2,
     +                        Wmbej_3,Wmbej_4,Wmbej_5,Wmbej_6,
     +                        Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T2old_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2old_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2old_ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_b)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension Fmi_a(Nocc_a,Nocc_a)
      Dimension Fmi_b(Nocc_b,Nocc_b)
      Dimension Fae_a(Nvrt_a,Nvrt_a)
      Dimension Fae_b(Nvrt_b,Nvrt_b)

      Dimension Wmnij_1(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension Wmnij_2(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension Wmnij_3(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Dimension Wabef_1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension Wabef_2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension Wabef_3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)
   
      Dimension Wmbej_1(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Wmbej_2(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Wmbej_3(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Wmbej_4(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension Wmbej_5(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension Wmbej_6(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Call T2_infmi(T2old_aa,T2old_bb,T2old_ab,Fmi_a,Fmi_b,
     +              W_aa,W_bb,W_ab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Call T2_infae(T2old_aa,T2old_bb,T2old_ab,Fae_a,Fae_b,
     +              W_aa,W_bb,W_ab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Call T2_inmnij(T2old_aa,T2old_bb,T2old_ab,W_aa,W_bb,
     +               W_ab,Wmnij_1,Wmnij_2,Wmnij_3,Nocc_a,
     +               Nocc_b,Nvrt_a,Nvrt_b)

      Call T2_inabef(T2old_aa,T2old_bb,T2old_ab,W_aa,W_bb,
     +               W_ab,Wabef_1,Wabef_2,Wabef_3,Nocc_a,
     +               Nocc_b,Nvrt_a,Nvrt_b)
  
      Call T2_inmbej(T2old_aa,T2old_bb,T2old_ab,W_aa,W_bb,
     +               W_ab,Wmbej_1,Wmbej_2,Wmbej_3,Wmbej_4,
     +               Wmbej_5,Wmbej_6,Nocc_a,Nocc_b,Nvrt_a,
     +               Nvrt_b)

      Return
      End



