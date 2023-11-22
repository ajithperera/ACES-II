      Subroutine Form_hbar_abci(Hbar_abci1,Hbar_abci2,Hbar_abci3,
     +                          Hbar_abci4,W4_aa,W4_bb,W4_ab,W4_ba,
     +                          W5_aa,W5_bb,W5_ab,W5_ba,Hbar_ijka1,
     +                          Hbar_ijka2,Hbar_ijka3,Hbar_ijka4,
     +                          Hbar_mea,Hbar_meb,Hbar_abef1,
     +                          Hbar_abef2,Hbar_abef3,
     +                          Wtt_mbej1,Wtt_mbej2,Wtt_mbej3,
     +                          Wtt_mbej4,Wtt_mbej5,Wtt_mbej6,
     +                          Hbar_aibc1,Hbar_aibc2,Hbar_aibc3,
     +                          Hbar_aibc4,T1aa,T1bb,T2aa,T2bb,T2ab,
     +                          Tau_aa,Tau_bb,Tau_ab,Nocc_a,Nocc_b,
     +                          Nvrt_a,Nvrt_b)
     
      Implicit Double Precision(A-H,O-Z)

      Dimension Wtt_mbej1(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Wtt_mbej2(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Wtt_mbej3(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Wtt_mbej4(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension Wtt_mbej5(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension Wtt_mbej6(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

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

      Dimension Hbar_abef1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension Hbar_abef2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension Hbar_abef3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Dimension Hbar_aibc1(Nvrt_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension Hbar_aibc2(Nvrt_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension Hbar_aibc3(Nvrt_a,Nocc_b,Nvrt_a,Nvrt_b)
      Dimension Hbar_aibc4(Nvrt_b,Nocc_a,Nvrt_b,Nvrt_a)

      Dimension Hbar_abci1(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Hbar_abci2(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Hbar_abci3(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Hbar_abci4(Nvrt_b,Nvrt_a,Nvrt_b,Nocc_a)

      dimension Hbar_ijka1(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension Hbar_ijka2(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension Hbar_ijka3(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension Hbar_ijka4(Nocc_b,Nocc_a,Nocc_b,Nvrt_a)

      Dimension Hbar_mea(Nocc_a,Nvrt_a)
      Dimension Hbar_meb(Nocc_b,Nvrt_b)

      Integer A,B,C,I

#include "ccsdlight_vars.com"

CSSS      Length_abci1 = Nvrt_a*Nvrt_a*Nvrt_a*Nocc_a
CSSS      Length_abci2 = Nvrt_b*Nvrt_b*Nvrt_b*Nocc_b
CSSS      Length_abci3 = Nvrt_a*Nvrt_b*Nvrt_a*Nocc_b
CSSS      Length_abci4 = Nvrt_b*Nvrt_a*Nvrt_b*Nocc_a
CSSS      Call Dcopy(Length_abci1,W4_aa,1,Hbar_abci1,1)
CSSS      Call Dcopy(Length_abci2,W4_bb,1,Hbar_abci2,1)
CSSS      Call Dcopy(Length_abci3,W4_ab,1,Hbar_abci3,1)

      Do I = 1, Nocc_a
      Do C = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         D = (1.0D0-Ocn_va(A))*(1.0D0-Ocn_va(B))*Ocn_oa(I)
         Hbar_abci1(A,B,C,I) = W4_aa(A,B,C,I)*D
      Enddo
      Enddo
      Enddo
      Enddo

      Do I = 1, Nocc_b
      Do C = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         D = (1.0D0-Ocn_vb(a))*(1.0D0-Ocn_vb(b))*Ocn_ob(i)
         Hbar_abci2(a,b,c,i) = W4_bb(a,b,c,i)*D
      Enddo
      Enddo
      Enddo
      Enddo

      Do I = 1, Nocc_b
      Do C = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         D = (1.0D0-Ocn_va(A))*(1.0D0-Ocn_vb(b))*Ocn_ob(i)
         Hbar_abci3(A,b,C,i) = W4_ab(A,b,C,i)*D
      Enddo
      Enddo
      Enddo
      Enddo

      Do I = 1, Nocc_a
      Do C = 1, Nvrt_b
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_b
         D = (1.0D0-Ocn_vb(a))*(1.0D0-Ocn_va(B))*Ocn_oa(I)
         Hbar_abci4(a,B,c,I) = W4_ba(B,a,I,c)*D
      Enddo
      Enddo 
      Enddo 
      Enddo 

C Form P(ab)T(be,im)*<am||ce>  

      Call Form_Hbabci_term1(Hbar_abci1,Hbar_abci2,Hbar_abci3,
     +                       Hbar_abci4,W4_aa,W4_bb,W4_ab,W4_ba,
     +                       Hbar_aibc1,Hbar_aibc2,Hbar_aibc3,
     +                       Hbar_aibc4,T2aa,T2bb,T2ab,Nocc_a,
     +                       Nocc_b,Nvrt_a,Nvrt_b)

C Form 1/2 Tau(ab,mn)*Hbar(ci,mn)
 
      Call Form_Hbabci_term2(Hbar_abci1,Hbar_abci2,Hbar_abci3,
     +                       Hbar_abci4,Hbar_ijka1,Hbar_ijka2,
     +                       Hbar_ijka3,Hbar_ijka4,Tau_aa,Tau_bb,
     +                       Tau_ab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C Form T(ab,im)*Hbar(m,e)

      Call Form_Hbabci_term3(Hbar_abci1,Hbar_abci2,Hbar_abci3,
     +                       Hbar_abci4,Hbar_mea,Hbar_meb,
     +                       T2aa,T2bb,T2ab,Nocc_a,Nocc_b,Nvrt_a,
     +                       Nvrt_b)

C Form T(e,i)*Hbar(ab,ce)

      Call Form_Hbabci_term4(Hbar_abci1,Hbar_abci2,Hbar_abci3,
     +                       Hbar_abci4,Hbar_abef1,Hbar_abef2,
     +                       Hbar_abef3,T1aa,T1bb,Nocc_a,Nocc_b,
     +                       Nvrt_a,Nvrt_b)

C Form -P(ab)T(a,m)[<mb||ci>-T(be,ni)*<mn||ce>]-> -P(ab)T(a,m)*Wtt(mb,ci)

      Call Form_hbabci_term5(Hbar_abci1,Hbar_abci2,Hbar_abci3,
     +                       Hbar_abci4,Wtt_mbej1,Wtt_mbej2,
     +                       Wtt_mbej3,Wtt_mbej4,Wtt_mbej5,
     +                       Wtt_mbej6,T1aa,T1bb,Nocc_a,Nocc_b,
     +                       Nvrt_a,Nvrt_b)

#ifdef _DEBUG_LVL0
      L_aaaa =Nvrt_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb =Nvrt_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab =Nvrt_a*Nocc_b*Nvrt_a*Nvrt_b
      L_baba =Nvrt_b*Nocc_a*Nvrt_b*Nvrt_a
      call checksum("Hbar_abci1:",Hbar_abci1,L_aaaa)
      call checksum("Hbar_abci2:",Hbar_abci2,L_bbbb)
      call checksum("Hbar_abci3:",Hbar_abci3,L_abab)
      call checksum("Hbar_abci4:",Hbar_abci4,L_baba)
#endif

      Return
      End
     


