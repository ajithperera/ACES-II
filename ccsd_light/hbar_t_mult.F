      Subroutine Hbar_t_mult(Hbar_mia,Hbar_mib,Hbar_aea,Hbar_aeb,
     +                       Hbar_mea,Hbar_meb,Hbar_ijkl1,
     +                       Hbar_ijkl2,Hbar_ijkl3,Hbar_abef1,
     +                       Hbar_abef2,Hbar_abef3,Hbar_mbej1, 
     +                       Hbar_mbej2,Hbar_mbej3,Hbar_mbej4,
     +                       Hbar_mbej5,Hbar_mbej6,Hbar_abci1,
     +                       Hbar_abci2,Hbar_abci3,Hbar_abci4,
     +                       Hbar_aibc1,Hbar_aibc2,Hbar_aibc3,
     +                       Hbar_aibc4,Hbar_ijka1,Hbar_ijka2,
     +                       Hbar_ijka3,Hbar_ijka4,Hbar_iajk1,
     +                       Hbar_iajk2,Hbar_iajk3,Hbar_iajk4,
     +                       Hbar_ijab1,Hbar_ijab2,Hbar_ijab3,
     +                       Ht1_aiaa,Ht1_aibb,Ht2_abij1,Ht2_abij2,
     +                       Ht2_abij3,T1aa,T1bb,T2aa,T2bb,T2ab,
     +                       T1old_aa,T1old_bb,T2old_aa,T2old_bb,
     +                       T2old_ab,Work,Maxcor,Nocc_a,Nocc_b,
     +                       Nvrt_a,Nvrt_b)

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

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)
      Dimension T1old_aa(Nvrt_a,Nocc_a)
      Dimension T1old_bb(Nvrt_b,Nocc_b)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)
      Dimension T2old_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2old_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2old_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Ht1_aiaa(Nvrt_a,Nocc_a)
      Dimension Ht1_aibb(Nvrt_b,Nocc_b)

      Dimension Ht2_abij1(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Ht2_abij2(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Ht2_abij3(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Data  Ione /1/

      L_aa = Nocc_a*Nvrt_a
      L_bb = Nocc_b*Nvrt_b
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b

      Call Dzero(Ht1_aiaa,L_aa)
      Call Dzero(Ht1_aibb,L_bb)
      Call Dzero(Ht2_abij1,L_aaaa)
      Call Dzero(Ht2_abij2,L_bbbb)
      Call Dzero(Ht2_abij3,L_abab)

#ifdef _DEBUG_LVL0
      Write(6,*)
      call checksum("T1aa      :",T1aa,L_aa)
      call checksum("T1bb      :",T1bb,L_bb)
      call checksum("T2aa      :",T2aa,L_aaaa)
      call checksum("T2bb      :",T2bb,L_bbbb)
      call checksum("T2ab      :",T2ab,L_abab)
#endif 

      Call Hbar_ss_t1_mult(Hbar_aea,Hbar_aeb,Hbar_mia,Hbar_mib,
     +                     Hbar_mbej1,Hbar_mbej2,Hbar_mbej3,
     +                     Hbar_mbej4,Hbar_mbej5,Hbar_mbej6,
     +                     T1aa,T1bb,Ht1_aiaa,Ht1_aibb,
     +                     Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Call Hbar_sd_t2_mult(Hbar_mea,Hbar_meb,Hbar_aibc1,Hbar_aibc2,
     +                     Hbar_aibc3,Hbar_aibc4,Hbar_ijka1,Hbar_ijka2,
     +                     Hbar_ijka3,Hbar_ijka4,T2aa,T2bb,T2ab,
     +                     Ht1_aiaa,Ht1_aibb,Nocc_a,Nocc_b,Nvrt_a,
     +                     Nvrt_b)
      Lenoo_aa = Nocc_a*Nocc_a
      Lenoo_bb = Nocc_b*Nocc_b
      Lenvv_aa = Nvrt_a*Nvrt_a
      Lenvv_bb = Nvrt_b*Nvrt_b

      Iq_ija = Ione
      Iq_ijb = Iq_ija + Lenoo_aa
      Iq_efa = Iq_ijb + Lenoo_bb
      Iq_efb = Iq_efa + Lenvv_aa
      Iend   = Iq_efb + Lenvv_bb
      If (Iend .Ge. Maxcor) Call Insmem("hbar_t_mult",Iend,Maxcor)

      Call Form_3body_intms1(Work(Iq_ija),Work(Iq_ijb),Work(Iq_efa),
     +                       Work(Iq_efb),Hbar_aibc1,Hbar_aibc2,
     +                       Hbar_aibc3,Hbar_aibc4,Hbar_ijka1,
     +                       Hbar_ijka2,Hbar_ijka3,Hbar_ijka4,T1aa,
     +                       T1bb,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Call Hbar_ds_t2_mult(Hbar_iajk1,Hbar_iajk2,Hbar_iajk3,Hbar_iajk4,
     +                     Hbar_abci1,Hbar_abci2,Hbar_abci3,Hbar_abci4,
     +                     T1aa,T1bb,T2aa,T2bb,T2ab,T2old_aa,T2old_bb,
     +                     T2old_ab,Ht2_abij1,Ht2_abij2,Ht2_abij3,
     +                     Work(Iq_ija),Work(Iq_ijb),Work(Iq_efa),
     +                     Work(Iq_efb),Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Call Form_3body_intms2(Work(Iq_ija),Work(Iq_ijb),Work(Iq_efa),
     +                       Work(Iq_efb),Hbar_ijab1,Hbar_ijab2,
     +                       Hbar_ijab3,T2aa,T2bb,T2ab,Nocc_a,Nocc_b,
     +                       Nvrt_a,Nvrt_b)

      Call Hbar_dd_t2_mult(Hbar_aea,Hbar_aeb,Hbar_mia,Hbar_mib,
     +                     Hbar_abef1,Hbar_abef2,Hbar_abef3,
     +                     Hbar_ijkl1,Hbar_ijkl2,Hbar_ijkl3,
     +                     Hbar_mbej1,Hbar_mbej2,Hbar_mbej3,
     +                     Hbar_mbej4,Hbar_mbej5,Hbar_mbej6,Ht2_abij1,
     +                     Ht2_abij2,Ht2_abij3,Work(Iq_ija),
     +                     Work(Iq_ijb),Work(Iq_efa),Work(Iq_efb),
     +                     T2aa,T2bb,T2ab,T2old_aa,T2old_bb,T2old_ab,
     +                     Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Return
      end 
