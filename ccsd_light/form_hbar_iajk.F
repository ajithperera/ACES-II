      Subroutine Form_hbar_iajk(Hbar_iajk1,Hbar_iajk2,Hbar_iajk3,
     +                          Hbar_iajk4,W4_aa,W4_bb,W4_ab,W4_ba,
     +                          W5_aa,W5_bb,W5_ab,W5_ba,
     +                          Hbar_aibc1,Hbar_aibc2,Hbar_aibc3,
     +                          Hbar_aibc4,Hbar_mea,Hbar_meb,
     +                          Hbar_ijkl1,Hbar_ijkl2,Hbar_ijkl3,
     +                          Wtt_mbej1,Wtt_mbej2,Wtt_mbej3,
     +                          Wtt_mbej4,Wtt_mbej5,Wtt_mbej6,T1aa,
     +                          T1bb,T2aa,T2bb,T2ab,Tau_aa,Tau_bb,
     +                          Tau_ab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)
     
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

      Dimension Hbar_iajk1(Nocc_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Hbar_iajk2(Nocc_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Hbar_iajk3(Nocc_a,Nvrt_b,Nocc_a,Nocc_b)
      Dimension Hbar_iajk4(Nocc_b,Nvrt_a,Nocc_b,Nocc_a)

      Dimension Hbar_aibc1(Nvrt_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension Hbar_aibc2(Nvrt_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension Hbar_aibc3(Nvrt_a,Nocc_b,Nvrt_a,Nvrt_b)
      Dimension Hbar_aibc4(Nvrt_b,Nocc_a,Nvrt_b,Nvrt_a)

      Dimension Hbar_mea(Nocc_a,Nvrt_a)
      Dimension Hbar_meb(Nocc_b,Nvrt_b)

      Integer A,I,J,K

#include "ccsdlight_vars.com"

      Do K = 1, Nocc_a
      Do J = 1, Nocc_a
      Do A = 1, Nvrt_a
      Do I = 1, Nocc_a
         C = Ocn_oa(J)*Ocn_oa(K)*(1.0D0-Ocn_va(A))
         Hbar_iajk1(I,A,J,K) = W5_aa(J,K,I,A)*C
      Enddo
      Enddo 
      Enddo 
      Enddo 

      Do K = 1, Nocc_b
      Do J = 1, Nocc_b
      Do A = 1, Nvrt_b
      Do I = 1, Nocc_b
         C = Ocn_ob(j)*Ocn_ob(k)*(1.0D0-Ocn_vb(a))
         Hbar_iajk2(i,a,j,k) = W5_bb(j,k,i,a)*C
      Enddo
      Enddo 
      Enddo 
      Enddo 

      Do K = 1, Nocc_b
      Do J = 1, Nocc_a
      Do A = 1, Nvrt_b
      Do I = 1, Nocc_a
         C =  Ocn_oa(J)*Ocn_ob(k)*(1.0D0-Ocn_vb(a))
         Hbar_iajk3(I,a,J,k) = W5_ab(I,k,J,a)*C 
      Enddo
      Enddo
      Enddo
      Enddo

      Do K = 1, Nocc_a
      Do J = 1, Nocc_b
      Do A = 1, Nvrt_a
      Do I = 1, Nocc_b
         C =  Ocn_ob(j)*Ocn_oa(K)*(1.0D0-Ocn_va(A))
         Hbar_iajk4(i,A,j,K) = W5_ba(K,i,A,j)*C
      Enddo
      Enddo
      Enddo
      Enddo

C Form P(ij)T(ae,km)*<im||je>  

      Call Form_Hbiajk_term1(Hbar_iajk1,Hbar_iajk2,Hbar_iajk3,
     +                       Hbar_iajk4,W5_aa,W5_bb,W5_ab,W5_ba,
     +                       T2aa,T2bb,T2ab,Nocc_a,Nocc_b,Nvrt_a,
     +                       Nvrt_b)

C Form 1/2 Tau(ef,jk)*W(ia,ef)
 
      Call Form_Hbiajk_term2(Hbar_iajk1,Hbar_iajk2,Hbar_iajk3,
     +                       Hbar_iajk4,W4_aa,W4_bb,W4_ab,
     +                       W4_ba,Tau_aa,Tau_bb,Tau_ab,
     +                       Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C Form T(ae,jk)*Hbar(i,e)

      Call Form_Hbiajk_term3(Hbar_iajk1,Hbar_iajk2,Hbar_iajk3,
     +                       Hbar_iajk4,Hbar_mea,Hbar_meb,
     +                       T2aa,T2bb,T2ab,Nocc_a,Nocc_b,Nvrt_a,
     +                       Nvrt_b)

C Form T(a,m)*Hbar(im,jk)

      Call Form_Hbiajk_term4(Hbar_iajk1,Hbar_iajk2,Hbar_iajk3,
     +                       Hbar_iajk4,Hbar_ijkl1,Hbar_ijkl2,
     +                       Hbar_ijkl3,Hbar_ijkl4,T1aa,T1bb,
     +                       Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C Form -P(ij)T(e,j)[<ia||ek>-T(af,mk)*<im||ef>]-> -P(ij)T(e,j)*Wtt(ia,ek)

      Call Form_hbiajk_term5(Hbar_iajk1,Hbar_iajk2,Hbar_iajk3,
     +                       Hbar_iajk4,Wtt_mbej1,Wtt_mbej2,
     +                       Wtt_mbej3,Wtt_mbej4,Wtt_mbej5,
     +                       Wtt_mbej6,T1aa,T1bb,Nocc_a,Nocc_b,
     +                       Nvrt_a,Nvrt_b)

#ifdef _DEBUG_LVL0
      L_aaaa =Nocc_a*Nvrt_a*Nocc_a*Nocc_a
      L_bbbb =Nocc_b*Nvrt_b*Nocc_b*Nocc_b
      L_abab =Nocc_a*Nvrt_b*Nocc_a*Nocc_b
      L_baba =Nocc_b*Nvrt_a*Nocc_b*Nocc_a 
      call checksum("Hbar_iajk1:",Hbar_iajk1,L_aaaa)
      call checksum("Hbar_iajk2:",Hbar_iajk2,L_bbbb)
      call checksum("Hbar_iajk3:",Hbar_iajk3,L_abab)
      call checksum("Hbar_iajk4:",Hbar_iajk4,L_baba)
#endif

      Return
      End
     


