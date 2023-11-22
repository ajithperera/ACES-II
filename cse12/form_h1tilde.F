      Subroutine Form_h1tilde(T1_aa,T1_bb,T1taup_aa,T1taup_bb,
     +                        T1taup_ab,Htmi_aa,Htmi_bb,Htae_aa,
     +                        Htae_bb,Htme_aa,Htme_bb,Htem_aa,
     +                        Htem_bb,W_aa,W_bb,W_ab,W3_aaaa,
     +                        W3_bbbb,W3_abab,W3_baba,W3_abba,
     +                        W3_baab,W4_aa,W4_bb,W4_ab,W4_ba,W5_aa,
     +                        W5_bb,W5_ab,W5_ba,Fockoo_a,Fockoo_b,
     +                        Fockvv_a,Fockvv_b,Fockoo_od_a,
     +                        Fockoo_od_b,Fockvv_od_a,Fockvv_od_b,
     +                        Fockov_a,Fockov_b,Nocc_a,Nocc_b,Nvrt_a,
     +                        Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T1_aa(Nvrt_a,Nocc_a)
      Dimension T1_bb(Nvrt_b,Nocc_b)

      Dimension T1taup_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T1taup_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T1taup_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)
    
      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

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

      Dimension Htmi_aa(Nocc_a,Nocc_a)
      Dimension Htmi_bb(Nocc_b,Nocc_b)
      Dimension Htae_aa(Nvrt_a,Nvrt_a)
      Dimension Htae_bb(Nvrt_b,Nvrt_b)
      Dimension Htme_aa(Nocc_a,Nvrt_a)
      Dimension Htme_bb(Nocc_b,Nvrt_b)
      Dimension Htem_aa(Nvrt_a,Nocc_a)
      Dimension Htem_bb(Nvrt_b,Nocc_b)

      Dimension Fockov_a(Nocc_a,Nvrt_a)
      Dimension Fockov_b(Nocc_b,Nvrt_b)
      Dimension Fockvv_od_a(Nvrt_a,Nvrt_a)
      Dimension Fockvv_od_b(Nvrt_b,Nvrt_b)
      Dimension Fockoo_od_a(Nocc_a,Nocc_a)
      Dimension Fockoo_od_b(Nocc_b,Nocc_b)
      Dimension Fockoo_a(Nocc_a,Nocc_a)
      Dimension Fockoo_b(Nocc_b,Nocc_b) 
      Dimension Fockvv_a(Nvrt_a,Nvrt_a)
      Dimension Fockvv_b(Nvrt_b,Nvrt_b)
 
#include "cse12_vars.com"
C For legacy reasons, inside these routine htmi,htae,htme are 
C referred to Fmi,Fae,Fme (as in the standard CC one-paticle 
C nomenclature).  

C Htmi(a,e) elements

      Call Form_htmi(T1taup_aa,T1taup_bb,T1taup_ab,T1_aa,T1_bb,
     +               Htmi_aa,Htmi_bb,W_aa,W_bb,W_ab,W5_aa,W5_bb,
     +               W5_ab,W5_ba,Fockoo_od_a,Fockoo_od_b,Fockov_a,
     +               Fockov_b,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C Htae(a,e) elements

      Call Form_htae(T1taup_aa,T1taup_bb,T1taup_ab,T1_aa,T1_bb,
     +               Htae_aa,Htae_bb,W_aa,W_bb,W_ab,W4_aa,W4_bb,
     +               W4_ab,W4_ba,Fockvv_od_a,Fockvv_od_b,Fockov_a,
     +               Fockov_b,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C Htme(m,e) elements

      Call Form_htme(T1_aa,T1_bb,Htme_aa,Htme_bb,W_aa,W_bb,W_ab,
     +               Fockov_a,Fockov_b,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

c Htem(e,m) elements

      Call Form_htem(T1_aa,T1_bb,Htmi_aa,Htmi_bb,Htae_aa,Htae_bb,
     +               Htme_aa,Htme_bb,Htem_aa,Htem_bb,W3_aaaa,W3_bbbb,
     +               W3_abab,W3_baba,W3_abba,W3_baab,Fockov_a,
     +               Fockov_b,Fockoo_a,Fockoo_b,Fockvv_a,Fockvv_b,
     +               Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Return
      End
