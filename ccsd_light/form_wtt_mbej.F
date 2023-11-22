      Subroutine Form_wtt_mbej(Wtt_mbej1,Wtt_mbej2,Wtt_mbej3,
     +                         Wtt_mbej4,Wtt_mbej5,Wtt_mbej6,
     +                         W3_aaaa,W3_bbbb,W3_abab,W3_baba,
     +                         W3_abba,W3_baab,Nocc_a,Nocc_b,
     +                         Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

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

      Data Twom /2.0D0/

      Integer M,N,A,F,E,B,J

#include "ccsdlight_vars.com"

      Lenmbej_1 = Nocc_a*Nvrt_a*Nvrt_a*Nocc_a
      Lenmbej_2 = Nocc_b*Nvrt_b*Nvrt_b*Nocc_b
      Lenmbej_3 = Nocc_a*Nvrt_b*Nvrt_a*Nocc_b
      Lenmbej_4 = Nocc_b*Nvrt_a*Nvrt_b*Nocc_a
      Lenmbej_5 = Nocc_a*Nvrt_b*Nocc_a*Nvrt_b
      Lenmbej_6 = Nocc_b*Nvrt_a*Nocc_b*Nvrt_a

      Call Dscal(Lenmbej_1,Twom,Wtt_mbej1,1)
      Call Dscal(Lenmbej_2,Twom,Wtt_mbej2,1)
      Call Dscal(Lenmbej_3,Twom,Wtt_mbej3,1)
      Call Dscal(Lenmbej_4,Twom,Wtt_mbej4,1)
      Call Dscal(Lenmbej_5,Twom,Wtt_mbej5,1)
      Call Dscal(Lenmbej_6,Twom,Wtt_mbej6,1)

      Do J = 1, Nocc_a
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do M = 1, Nocc_a
         C = (1.0D0-Ocn_va(B))*Ocn_oa(J)
         Wtt_mbej1(M,B,E,J) = Wtt_mbej1(M,B,E,J) + W3_aaaa(M,B,E,J)*C
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_b
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do M = 1, Nocc_b
         C = (1.0D0-Ocn_vb(b))*Ocn_ob(j)
         Wtt_mbej2(m,b,e,j) = Wtt_mbej2(m,b,e,j)+ W3_bbbb(m,b,e,j)*C
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do M = 1, Nocc_a
         C = (1.0D0-Ocn_vb(b))*Ocn_ob(j)
         Wtt_mbej3(M,b,E,j) = Wtt_mbej3(M,b,E,j) + W3_abab(M,b,E,j)*C
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_a
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_a
      Do M = 1, Nocc_b
         C = (1.0D0-Ocn_va(B))*Ocn_oa(J)
         Wtt_mbej4(m,B,e,J) = Wtt_mbej4(m,B,e,J) + W3_baba(m,B,e,J)*C
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_a
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do M = 1, Nocc_a
         C = (1.0D0-Ocn_vb(b))*Ocn_oa(J)
         Wtt_mbej5(M,b,J,e) = Wtt_mbej5(M,b,J,e) + W3_abba(M,b,J,e)*C
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do M = 1, Nocc_b
         C = (1.0D0-Ocn_va(B))*Ocn_ob(j)
         Wtt_mbej6(m,B,j,E) = Wtt_mbej6(m,B,j,E) + W3_baab(m,B,j,E)*C
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      call checksum("Form_wtti1:",Wtt_mbej1,Lenmbej_1)
      call checksum("Form_wtti2:",Wtt_mbej2,Lenmbej_2)
      call checksum("Form_wtti3:",Wtt_mbej3,Lenmbej_3)
      call checksum("Form_wtti4:",Wtt_mbej4,Lenmbej_4)
      call checksum("Form_Wtti5:",Wtt_mbej5,Lenmbej_5)
      call checksum("Form_Wtti6:",Wtt_mbej6,Lenmbej_6)
#endif
      Return
      EnD








