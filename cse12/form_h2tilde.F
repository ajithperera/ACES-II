      SUbroutine Form_h2tilde(T1_aa,T1_bb,T1tau_aa,T1tau_bb,T1tau_ab,
     +                        W0_aa,W0_bb,W0_ab,W1_aa,W1_bb,W1_ab,
     +                        W_aa,W_bb,W_ab,W2_aa,W2_bb,W2_ab,
     +                        W3_aaaa,W3_bbbb,W3_abab,W3_baba,W3_abba,
     +                        W3_baab,W4_aa,W4_bb,W4_ab,W4_ba,W5_aa,
     +                        W5_bb,W5_ab,W5_ba,Htmnij_1,Htmnij_2,
     +                        Htmnij_3,Htabef_1,Htabef_2,Htabef_3,
     +                        Htmbej_1,Htmbej_2,Htmbej_3,Htmbej_4,
     +                        Htmbej_5,Htmbej_6,Htabij_1,Htabij_2,
     +                        Htabij_3,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T1_aa(Nvrt_a,Nocc_a)
      Dimension T1_bb(Nvrt_b,Nocc_b)

      Dimension T1tau_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T1tau_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T1tau_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension W0_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension W0_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension W0_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension W1_aa(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension W1_bb(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension W1_ab(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Dimension W2_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension W2_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension W2_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Dimension W4_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W4_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W4_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W4_ba(Nvrt_a,Nvrt_b,Nocc_a,Nvrt_b)

      Dimension W5_aa(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension W5_bb(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension W5_ab(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension W5_ba(Nocc_a,Nocc_b,Nvrt_a,Nocc_b)

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

      Dimension Htabij_1(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Htabij_2(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Htabij_3(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Htmbej_1(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Htmbej_2(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Htmbej_3(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Htmbej_4(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension Htmbej_5(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension Htmbej_6(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)
 
      Data Dtwo /2.0D0/

C Ht(mn,ij) elements 

      Call form_htmnij(T1_aa,T1_bb,T1tau_aa,T1tau_bb,T1tau_ab,W1_aa,
     +                 W1_bb,W1_ab,W_aa,W_bb,W_ab,W5_aa,W5_bb,W5_ab,
     +                 W5_ba,Htmnij_1,Htmnij_2,Htmnij_3,Nocc_a,Nocc_b,
     +                 Nvrt_a,Nvrt_b)

C Ht(ab,ef) elements 

      Call form_htabef(T1_aa,T1_bb,T1tau_aa,T1tau_bb,T1tau_ab,W2_aa,
     +                 W2_bb,W2_ab,W_aa,W_bb,W_ab,W4_aa,W4_bb,W4_ab,
     +                 W4_ba,Htabef_1,Htabef_2,Htabef_3,Nocc_a,Nocc_b,
     +                 Nvrt_a,Nvrt_b)

C Ht(mb,ej) elements 

      Call form_htmbej(T1_aa,T1_bb,W_aa,W_bb,W_ab,W4_aa,W4_bb,
     +                 W4_ab,W4_ba,W5_aa,W5_bb,W5_ab,W5_ba,
     +                 W3_aaaa,W3_bbbb,W3_abab,W3_baba,W3_abba,
     +                 W3_baab,Htmbej_1,Htmbej_2,Htmbej_3,
     +                 Htmbej_4,Htmbej_5,Htmbej_6,Nocc_a,
     +                 Nocc_b,Nvrt_a,Nvrt_b)

C Ht(ab,ij) elements 

      Call form_htabij(T1_aa,T1_bb,T1tau_aa,T1tau_bb,T1tau_ab,Htmnij_1,
     +                 Htmnij_2,Htmnij_3,Htabef_1,Htabef_2,Htabef_3,
     +                 Htmbej_1,Htmbej_2,Htmbej_3,Htmbej_4,Htmbej_5,
     +                 Htmbej_6,W0_aa,W0_bb,W0_ab,W3_aaaa,W3_bbbb,
     +                 W3_abab,W3_baba,W3_abba,W3_baab,W4_aa,W4_bb,
     +                 W4_ab,W4_ba,W5_aa,W5_bb,W5_ab,W5_ba,Htabij_1,
     +                 Htabij_2,Htabij_3,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C Modify the Htmnij and Htabef so that the last term (T1*T1) has the
C factor of 1/2 instead of 1/4. It is acomplished by multiplling the
C entire term by 2 and then by fixing the others that should not 
C have been multiplied by 2. 

      Loooo_aa = Nocc_a*Nocc_a*Nocc_a*Nocc_a
      Loooo_bb = Nocc_b*Nocc_b*Nocc_b*Nocc_b
      Loooo_ab = Nocc_a*Nocc_b*Nocc_a*Nocc_b

      Call Dscal(Loooo_aa,Dtwo,Htmnij_1,1)
      Call Dscal(Loooo_bb,Dtwo,Htmnij_2,1)
      Call Dscal(Loooo_ab,Dtwo,Htmnij_3,1)

      Lvvvv_aa = Nvrt_a*Nvrt_a*Nvrt_a*Nvrt_a
      Lvvvv_bb = Nvrt_b*Nvrt_b*Nvrt_b*Nvrt_b
      Lvvvv_ab = Nvrt_a*Nvrt_b*Nvrt_a*Nvrt_b

      Call Dscal(Lvvvv_aa,Dtwo,Htabef_1,1)
      Call Dscal(Lvvvv_bb,Dtwo,Htabef_2,1)
      Call Dscal(Lvvvv_ab,Dtwo,Htabef_3,1)

      Call Modf_htmnij(T1_aa,T1_bb,W1_aa,W1_bb,W1_ab,W5_aa,W5_bb,W5_ab,
     +                 W5_ba,Htmnij_1,Htmnij_2,Htmnij_3,Nocc_a,Nocc_b,
     +                 Nvrt_a,Nvrt_b)

      Call Modf_htabef(T1_aa,T1_bb,W2_aa,W2_bb,W2_ab,W4_aa,W4_bb,W4_ab,
     +                 W4_ba,Htabef_1,Htabef_2,Htabef_3,Nocc_a,Nocc_b,
     +                 Nvrt_a,Nvrt_b)

      Return
      End 


