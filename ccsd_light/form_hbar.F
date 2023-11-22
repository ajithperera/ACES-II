      Subroutine form_hbar(Hbar_mia,Hbar_mib,Hbar_aea,Hbar_aeb,
     +                     Hbar_mea,Hbar_meb,Hbar_ijkl1,
     +                     Hbar_ijkl2,Hbar_ijkl3,Hbar_abef1,
     +                     Hbar_abef2,Hbar_abef3,Hbar_mbej1, 
     +                     Hbar_mbej2,Hbar_mbej3,Hbar_mbej4,
     +                     Hbar_mbej5,Hbar_mbej6,Hbar_abci1,
     +                     Hbar_abci2,Hbar_abci3,Hbar_abci4,
     +                     Hbar_aibc1,Hbar_aibc2,Hbar_aibc3,
     +                     Hbar_aibc4,Hbar_ijka1,Hbar_ijka2,
     +                     Hbar_ijka3,Hbar_ijka4,Hbar_iajk1,
     +                     Hbar_iajk2,Hbar_iajk3,Hbar_iajk4,
     +                     Hbar_ijab1,Hbar_ijab2,Hbar_ijab3,
     +                     T1aa,T1bb,T2aa,T2bb,T2ab,Tau_aa,
     +                     Tau_bb,Tau_ab,w_aa,W_bb,W_ab,
     +                     W2_aa,W2_bb,W2_ab,W3_aaaa,W3_bbbb,
     +                     W3_abab,W3_baba,W3_abba,W3_baab,
     +                     W4_aa,W4_bb,W4_ab,W4_ba,
     +                     W5_aa,W5_bb,W5_ab,W5_ba,
     +                     Fockoo_a,Fockoo_b,Fockvv_a,
     +                     Fockvv_b,Fockov_a,Fockov_b,
     +                     Work,Maxcor,Nocc_a,Nocc_b,Nvrt_a,
     +                     Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

      Dimension Hbar_abef1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension Hbar_abef2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension Hbar_abef3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Dimension Hbar_ijkl1(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension Hbar_ijkl2(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension Hbar_ijkl3(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Dimension Hbar_ijab1(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension Hbar_ijab2(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension Hbar_ijab3(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension Hbar_ijka1(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension Hbar_ijka2(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension Hbar_ijka3(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension Hbar_ijka4(Nocc_b,Nocc_a,Nocc_b,Nvrt_a)

      Dimension Hbar_iajk1(Nocc_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Hbar_iajk2(Nocc_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Hbar_iajk3(Nocc_a,Nvrt_b,Nocc_a,Nocc_b)
      Dimension Hbar_iajk4(Nocc_b,Nvrt_a,Nocc_b,Nocc_a)

      Dimension Hbar_abci1(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Hbar_abci2(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Hbar_abci3(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Hbar_abci4(Nvrt_b,Nvrt_a,Nvrt_b,Nocc_a)

      Dimension Hbar_aibc1(Nvrt_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension Hbar_aibc2(Nvrt_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension Hbar_aibc3(Nvrt_a,Nocc_b,Nvrt_a,Nvrt_b)
      Dimension Hbar_aibc4(Nvrt_b,Nocc_a,Nvrt_b,Nvrt_a)

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
      Dimension Hbar_mea(Nocc_a,Nvrt_a)
      Dimension Hbar_meb(Nocc_b,Nvrt_b)

      Dimension Fockoo_a(Nocc_a,Nocc_a)
      Dimension Fockoo_b(Nocc_b,Nocc_b)
      Dimension Fockvv_a(Nvrt_a,Nvrt_a)
      Dimension Fockvv_b(Nvrt_b,Nvrt_b)
      Dimension Fockov_a(Nocc_a,Nvrt_a)
      Dimension Fockov_b(Nocc_b,Nvrt_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_b)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension W2_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension W2_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension W2_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Dimension W3_aaaa(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W3_bbbb(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W3_abab(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W3_baba(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension W3_abba(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension W3_baab(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Dimension W4_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W4_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W4_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W4_ba(Nvrt_a,Nvrt_b,Nocc_a,Nvrt_b)

      Dimension W5_aa(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension W5_bb(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension W5_ab(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension W5_ba(Nocc_a,Nocc_b,Nvrt_a,Nocc_b)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension Tau_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Tau_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Tau_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Data  Ione /1/

      Lenijkl_1 = Nocc_a*Nocc_a*Nocc_a*Nocc_a
      Lenijkl_2 = Nocc_b*Nocc_b*Nocc_b*Nocc_b
      Lenijkl_3 = Nocc_a*Nocc_b*Nocc_a*Nocc_b

      Lenabef_1 = Nvrt_a*Nvrt_a*Nvrt_a*Nvrt_a
      Lenabef_2 = Nvrt_b*Nvrt_b*Nvrt_b*Nvrt_b
      Lenabef_3 = Nvrt_a*Nvrt_b*Nvrt_a*Nvrt_b

      Lenmbej_1 = Nocc_a*Nvrt_a*Nvrt_a*Nocc_a
      Lenmbej_2 = Nocc_b*Nvrt_b*Nvrt_b*Nocc_b
      Lenmbej_3 = Nocc_a*Nvrt_b*Nvrt_a*Nocc_b
      Lenmbej_4 = Nocc_b*Nvrt_a*Nvrt_b*Nocc_a
      Lenmbej_5 = Nocc_a*Nvrt_b*Nocc_a*Nvrt_b
      Lenmbej_6 = Nocc_b*Nvrt_a*Nocc_b*Nvrt_a

C First fix Hbar_mi and Hbar_ae by adding orbital energies. Also
C form the twidle-twidle intermediate (<mb||ej> + t(fb,jn)<mn||ef>)

      Call Fix_hbarpq(Hbar_mia,Hbar_mib,Hbar_aea,Hbar_aeb,
     +                Hbar_mea,Hbar_meb,Fockoo_a,Fockoo_b,
     +                Fockvv_a,Fockvv_b,Nocc_a,Nocc_b,
     +                Nvrt_a,Nvrt_b)

      It2ijkl_1 = Ione
      It2ijkl_2 = It2ijkl_1 + Lenijkl_1
      It2ijkl_3 = It2ijkl_2 + Lenijkl_2
      Iend      = It2ijkl_3 + Lenijkl_3 
      If (Iend .Ge. Maxcor) Call Insmem("fix_hbar_ijkl",Iend,Maxcor)

C Fix the W_mnij intermediate to get Hbar_ijkl with 1/2Tau(ef,kl)*<kl||ef>

      Call Fix_hbar_ijkl(Hbar_ijkl1,Hbar_ijkl2,Hbar_ijkl3,
     +                   Work(It2ijkl_1),Work(It2ijkl_2),
     +                   Work(It2ijkl_3),W_aa,W_bb,W_ab,
     +                   Tau_aa,Tau_bb,Tau_ab,Nocc_a,Nocc_b,
     +                   Nvrt_a,Nvrt_b)

      It2abef_1 = Ione
      It2abef_2 = It2abef_1 + Lenabef_1
      It2abef_3 = It2abef_2 + Lenabef_2
      Iend      = It2abef_3 + Lenabef_3
      If (Iend .Ge. Maxcor) Call Insmem("fix_hbar_abef",Iend,Maxcor)

C Fix the W_abef intermediate to get Hbar_abef with 1/2Tau(ab,mn)*<mn||ef>

      Call Fix_hbar_abef(Hbar_abef1,Hbar_abef2,Hbar_abef3,
     +                   Work(It2abef_1),Work(It2abef_2),
     +                   Work(It2abef_3),W_aa,W_bb,W_ab,
     +                   Tau_aa,Tau_bb,Tau_ab,Nocc_a,Nocc_b,
     +                   Nvrt_a,Nvrt_b)

      It2mbej_1 = Ione
      It2mbej_2 = It2mbej_1 + Lenmbej_1
      It2mbej_3 = It2mbej_2 + Lenmbej_2
      It2mbej_4 = It2mbej_3 + Lenmbej_3
      It2mbej_5 = It2mbej_4 + Lenmbej_4
      It2mbej_6 = It2mbej_5 + Lenmbej_5
      Iend      = It2mbej_6 + Lenmbej_6
      If (Iend .Ge. Maxcor) Call Insmem("fix_hbar_mbej",Iend,Maxcor)

      Call Fix_hbar_mbej(Hbar_mbej1,Hbar_mbej2,Hbar_mbej3,Hbar_mbej4,
     +                   Hbar_mbej5,Hbar_mbej6,Work(It2mbej_1),
     +                   Work(It2mbej_2),Work(It2mbej_3),
     +                   Work(It2mbej_4),Work(It2mbej_5),
     +                   Work(It2mbej_6),W3_aaaa,W3_bbbb,
     +                   W3_abab,W3_baba,W3_abba,W3_baab,
     +                   W_aa,W_bb,W_ab,T2aa,T2bb,T2ab,
     +                   Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C For -P(ab)T1(b,m)<am|cd> term. This is part of W(abcd) that we needs
C when we compute Hbar(ab,ci).

      It2abef_1 = Iend 
      It2abef_2 = It2abef_1 + Lenabef_1
      It2abef_3 = It2abef_2 + Lenabef_2
      Iend      = It2abef_3 + Lenabef_3
      If (Iend .Ge. Maxcor) Call Insmem("Form_Wabef_intms",Iend,
     +                                   Maxcor)

      Call Form_Wabef_intms(T1aa,T1bb,W4_aa,W4_bb,W4_ab,W4_ba,
     +                      Work(It2abef_1),
     +                      Work(It2abef_2),Work(It2abef_3),
     +                      Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Call Form_hbar_ijka(Hbar_ijka1,Hbar_ijka2,Hbar_ijka3,Hbar_ijka4,
     +                    T1aa,T1bb,W_aa,W_bb,W_ab,W5_aa,W5_bb,W5_ab,
     +                    W5_ba,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C This is simply to get the abci stored integrals sorted into aibc order
C to the fist contraction to Hbar(ab,ci)

      Scale = 0.0D0
      Call Form_hbar_aibc(Hbar_aibc1,Hbar_aibc2,Hbar_aibc3,Hbar_aibc4,
     +                    T1aa,T1bb,W_aa,W_bb,W_ab,W4_aa,W4_bb,W4_ab,
     +                    W4_ba,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,Scale)

      Call Form_hbar_abci(Hbar_abci1,Hbar_abci2,Hbar_abci3,Hbar_abci4,
     +                    W4_aa,W4_bb,W4_ab,W4_ba,W5_aa,W5_bb,W5_ab,
     +                    W5_ba,Hbar_ijka1,Hbar_ijka2,Hbar_ijka3,
     +                    Hbar_ijka4,Hbar_mea,Hbar_meb,
     +                    Work(It2abef_1),Work(It2abef_2),
     +                    Work(It2abef_3),
     +                    Work(It2mbej_1),Work(It2mbej_2),
     +                    Work(It2mbej_3),Work(It2mbej_4),
     +                    Work(It2mbej_5),Work(It2mbej_6),
     +                    Hbar_aibc1,Hbar_aibc2,Hbar_aibc3,
     +                    Hbar_aibc4,T1aa,T1bb,T2aa,T2bb,T2ab,
     +                    Tau_aa,Tau_bb,Tau_ab,Nocc_a,Nocc_b,
     +                    Nvrt_a,Nvrt_b)

      Call Form_t1_2hbar_abci(Hbar_abci1,Hbar_abci2,Hbar_abci3,
     +                        Hbar_abci4,T1aa,T1bb,W2_aa,W2_bb,
     +                        W2_ab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C This actuall construct the Hbar(ai,bc)

      Scale = 1.0D0
      Call Form_hbar_aibc(Hbar_aibc1,Hbar_aibc2,Hbar_aibc3,Hbar_aibc4,
     +                    T1aa,T1bb,W_aa,W_bb,W_ab,W4_aa,W4_bb,W4_ab,
     +                    W4_ba,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,Scale)

      Call Form_hbar_iajk(Hbar_iajk1,Hbar_iajk2,Hbar_iajk3,Hbar_iajk4,
     +                    W4_aa,W4_bb,W4_ab,W4_ba,W5_aa,W5_bb,W5_ab,
     +                    W5_ba,Hbar_aibc1,Hbar_aibc2,Hbar_aibc3,
     +                    Hbar_aibc4,Hbar_mea,Hbar_meb,Hbar_ijkl1,
     +                    Hbar_ijkl2,Hbar_ijkl3,Work(It2mbej_1),
     +                    Work(It2mbej_2),Work(It2mbej_3),
     +                    Work(It2mbej_4),Work(It2mbej_5),
     +                    Work(It2mbej_6),T1aa,T1bb,T2aa,T2bb,
     +                    T2ab,Tau_aa,Tau_bb,Tau_ab,Nocc_a,
     +                    Nocc_b,Nvrt_a,Nvrt_b)

      Call Form_hbar_ijab(Hbar_ijab1,Hbar_ijab2,Hbar_ijab3,
     +                    W_aa,W_bb,W_ab,Nocc_a,Nocc_b,Nvrt_a,
     +                    Nvrt_b)

      Return
      End 
