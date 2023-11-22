      Subroutine Form_htabij(T1_aa,T1_bb,T1tau_aa,T1tau_bb,T1tau_ab,
     +                       Htmnij_1,Htmnij_2,Htmnij_3,Htabef_1,
     +                       Htabef_2,Htabef_3,Htmbej_1,Htmbej_2,
     +                       Htmbej_3,Htmbej_4,Htmbej_5,Htmbej_6,
     +                       W0_aa,W0_bb,W0_ab,W3_aaaa,W3_bbbb,
     +                       W3_abab,W3_baba,W3_abba,W3_baab,W4_aa,
     +                       W4_bb,W4_ab,W4_ba,W5_aa,W5_bb,W5_ab,
     +                       W5_ba,Htabij_1,Htabij_2,Htabij_3,
     +                       Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T1_aa(Nvrt_a,Nocc_a)
      Dimension T1_bb(Nvrt_b,Nocc_b)

      Dimension T1tau_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T1tau_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T1tau_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension W0_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension W0_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension W0_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension W3_aaaa(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W3_bbbb(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W3_abab(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W3_baba(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension W3_abba(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension W3_baab(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Dimension Htmnij_1(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension Htmnij_2(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension Htmnij_3(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Dimension Htabef_1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension Htabef_2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension Htabef_3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Dimension W4_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W4_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W4_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W4_ba(Nvrt_a,Nvrt_b,Nocc_a,Nvrt_b)
    
      Dimension W5_aa(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension W5_bb(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension W5_ab(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension W5_ba(Nocc_a,Nocc_b,Nvrt_a,Nocc_b)

      Dimension Htabij_1(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Htabij_2(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Htabij_3(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Call Wabij_inhtabij(W0_aa,W0_bb,W0_ab,Htabij_1,Htabij_2,
     +                   Htabij_3,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Call T1t1_inhtabij(T1_aa,T1_bb,Htabij_1,Htabij_2,Htabij_3,
     +                   W3_aaaa,W3_bbbb,W3_abab,W3_baba,W3_abba,
     +                   W3_baab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)
 
      Call T1_inhtabij(T1_aa,T1_bb,Htabij_1,Htabij_2,Htabij_3,
     +                 W4_aa,W4_bb,W4_ab,W4_ba,W5_aa,W5_bb,
     +                 W5_ab,W5_ba,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Call Htmnij_inhtabij(T1tau_aa,T1tau_bb,T1tau_ab,Htabij_1,
     +                     Htabij_2,Htabij_3,Htmnij_1,Htmnij_2,
     +                     Htmnij_3,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Call Htabef_inhtabij(T1tau_aa,T1tau_bb,T1tau_ab,Htabij_1,
     +                     Htabij_2,Htabij_3,Htabef_1,Htabef_2,
     +                     Htabef_3,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Return
      End  
