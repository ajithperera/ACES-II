      Subroutine Fix_hbar_mbej(Hbar_mbej1,Hbar_mbej2,Hbar_mbej3,
     +                         Hbar_mbej4,Hbar_mbej5,Hbar_mbej6,
     +                         Wtt_mbej1,Wtt_mbej2,Wtt_mbej3,
     +                         Wtt_mbej4,Wtt_mbej5,Wtt_mbej6,
     +                         W3_aaaa,W3_bbbb,W3_abab,W3_baba,
     +                         W3_abba,W3_baab,W_aa,W_bb,W_ab,
     +                         T2aa,T2bb,T2ab,Nocc_a,Nocc_b,Nvrt_a,
     +                         Nvrt_b)

      Implicit Double Precision(A-H,O-Z) 

      Data Ione,One,Onem /1,1.0D0,-1.00D0/

      Dimension Hbar_mbej1(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Hbar_mbej2(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Hbar_mbej3(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Hbar_mbej4(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension Hbar_mbej5(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension Hbar_mbej6(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Dimension W3_aaaa(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W3_bbbb(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W3_abab(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W3_baba(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension W3_abba(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension W3_baab(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Dimension Wtt_mbej1(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Wtt_mbej2(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Wtt_mbej3(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Wtt_mbej4(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension Wtt_mbej5(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension Wtt_mbej6(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C First form the -1/2t(fb,jn)<mn||ef> contribution.

      Call Form_t2mbej(T2aa,T2bb,T2ab,W_aa,W_bb,W_ab,Wtt_mbej1,
     +                 Wtt_mbej2,Wtt_mbej3,Wtt_mbej4,Wtt_mbej5,
     +                 Wtt_mbej6,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C Now add that to the existing Hbar to get -t(fb,jn)<mn||ef> 

      Lenmbej_1 = Nocc_a*Nvrt_a*Nvrt_a*Nocc_a
      Lenmbej_2 = Nocc_b*Nvrt_b*Nvrt_b*Nocc_b
      Lenmbej_3 = Nocc_a*Nvrt_b*Nvrt_a*Nocc_b
      Lenmbej_4 = Nocc_b*Nvrt_a*Nvrt_b*Nocc_a
      Lenmbej_5 = Nocc_a*Nvrt_b*Nocc_a*Nvrt_b
      Lenmbej_6 = Nocc_b*Nvrt_a*Nocc_b*Nvrt_a

      Call Daxpy(Lenmbej_1,One,Wtt_mbej1,1,Hbar_mbej1,1)
      Call Daxpy(Lenmbej_2,One,Wtt_mbej2,1,Hbar_mbej2,1)
      Call Daxpy(Lenmbej_3,One,Wtt_mbej3,1,Hbar_mbej3,1)
      Call Daxpy(Lenmbej_4,One,Wtt_mbej4,1,Hbar_mbej4,1)
      Call Daxpy(Lenmbej_5,One,Wtt_mbej5,1,Hbar_mbej5,1)
      Call Daxpy(Lenmbej_6,One,Wtt_mbej6,1,Hbar_mbej6,1)

#ifdef _DEBUG_LVL0
      call checksum("Form_wMBEJ:",Hbar_mbej1,Lenmbej_1)
      call checksum("Form_wmbej:",Hbar_mbej2,Lenmbej_2)
      call checksum("Form_wMbEj:",Hbar_mbej3,Lenmbej_3)
      call checksum("Form_wmBeJ:",Hbar_mbej4,Lenmbej_4)
      call checksum("Form_wmBEj:",Hbar_mbej5,Lenmbej_5)
      call checksum("Form_wMbeJ:",Hbar_mbej6,Lenmbej_6)
#endif

C Now form the  <mb||ej> + t(fb,jn)<mn||ef> 

      Call Form_wtt_mbej(Wtt_mbej1,Wtt_mbej2,Wtt_mbej3,
     +                   Wtt_mbej4,Wtt_mbej5,Wtt_mbej6,
     +                   W3_aaaa,W3_bbbb,W3_abab,W3_baba,
     +                   W3_abba,W3_baab,Nocc_a,Nocc_b,
     +                   Nvrt_a,Nvrt_b)
    
      Return
      End
      
  
