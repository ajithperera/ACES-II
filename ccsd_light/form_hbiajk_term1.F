      Subroutine Form_hbiajk_term1(Hbar_iajk1,Hbar_iajk2,Hbar_iajk3,
     +                             Hbar_iajk4,W5_aa,W5_bb,W5_ab,W5_ba,
     +                             T2aa,T2bb,T2ab,Nocc_a,Nocc_b,Nvrt_a,
     +                             Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension W5_aa(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension W5_bb(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension W5_ab(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension W5_ba(Nocc_a,Nocc_b,Nvrt_a,Nocc_b)

      Dimension Hbar_iajk1(Nocc_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Hbar_iajk2(Nocc_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Hbar_iajk3(Nocc_a,Nvrt_b,Nocc_a,Nocc_b)
      Dimension Hbar_iajk4(Nocc_b,Nvrt_a,Nocc_b,Nocc_a)

      Integer I,J,K,M,E,A

#include "ccsdlight_vars.com"

C AAAA contribution
      
      Do K = 1, Nocc_a
      Do J = 1, Nocc_a
      Do A = 1, Nvrt_a
      Do I = 1, Nocc_a
         T = 0.0D0
         C = Ocn_oa(J)
      Do E = 1, Nvrt_a
      Do M = 1, Nocc_a
         T = T + T2aa(A,E,K,M)*W5_aa(I,M,J,E)*C
      Enddo
      Enddo
         Hbar_iajk1(I,A,J,K) = Hbar_iajk1(I,A,J,K) + T
         Hbar_iajk1(I,A,K,J) = Hbar_iajk1(I,A,K,J) - T
      Enddo
      Enddo
      Enddo
      Enddo

      Do K = 1, Nocc_a
      Do J = 1, Nocc_a
      Do A = 1, Nvrt_a
      Do I = 1, Nocc_a
         T = 0.0D0
         C = Ocn_oa(J)
      Do E = 1, Nvrt_b
      Do M = 1, Nocc_b
         T = T + T2ab(A,e,K,m)*W5_ab(I,m,J,e)*C
      Enddo
      Enddo
         Hbar_iajk1(I,A,J,K) = Hbar_iajk1(I,A,J,K) + T
         Hbar_iajk1(I,A,K,J) = Hbar_iajk1(I,A,K,J) - T
      Enddo
      Enddo
      Enddo
      Enddo

C BBBB contribution

      Do K = 1, Nocc_b
      Do J = 1, Nocc_b
      Do A = 1, Nvrt_b
      Do I = 1, Nocc_b
         T = 0.0D0
         C = Ocn_ob(j)
      Do E = 1, Nvrt_b
      Do M = 1, Nocc_b
         T = T + T2bb(a,e,k,m)*W5_bb(i,m,j,e)*C
      Enddo
      Enddo
         Hbar_iajk2(i,a,j,k) = Hbar_iajk2(i,a,j,k) + T
         Hbar_iajk2(i,a,k,j) = Hbar_iajk2(i,a,k,j) - T
      Enddo
      Enddo
      Enddo
      Enddo

      Do K = 1, Nocc_b
      Do J = 1, Nocc_b
      Do A = 1, Nvrt_b
      Do I = 1, Nocc_b
         T = 0.0D0
         C = Ocn_ob(j)
      Do E = 1, Nvrt_a
      Do M = 1, Nocc_a
         T = T + T2ab(E,a,M,k)*W5_ba(M,i,E,j)*C
      Enddo
      Enddo
         Hbar_iajk2(i,a,j,k) = Hbar_iajk2(i,a,j,k) + T
         Hbar_iajk2(i,a,k,j) = Hbar_iajk2(i,a,k,j) - T
      Enddo
      Enddo
      Enddo
      Enddo

C ABAB contribution

      Do K = 1, Nocc_b
      Do J = 1, Nocc_a
      Do A = 1, Nvrt_b
      Do I = 1, Nocc_a
         T = 0.0D0
         C = Ocn_oa(J)
      Do E = 1, Nvrt_b
      Do M = 1, Nocc_b
         T = T + T2bb(a,e,k,m)*W5_ab(I,m,J,e)*C
      Enddo
      Enddo
         Hbar_iajk3(I,a,J,k) = Hbar_iajk3(I,a,J,k) + T
      Enddo
      Enddo
      Enddo
      Enddo

      Do K = 1, Nocc_b
      Do J = 1, Nocc_a
      Do A = 1, Nvrt_b
      Do I = 1, Nocc_a
         T = 0.0D0
         C = Ocn_oa(J)
      Do E = 1, Nvrt_a
      Do M = 1, Nocc_a
         T = T + T2ab(E,a,M,k)*W5_aa(I,M,J,E)*C
      Enddo
      Enddo
         Hbar_iajk3(I,a,J,k) = Hbar_iajk3(I,a,J,k) + T
      Enddo
      Enddo
      Enddo
      Enddo

      Do K = 1, Nocc_b
      Do J = 1, Nocc_a
      Do A = 1, Nvrt_b
      Do I = 1, Nocc_a
         T = 0.0D0
         C = Ocn_ob(k)
      Do E = 1, Nvrt_a
      Do M = 1, Nocc_b
         T = T + T2ab(E,a,J,m)*W5_ba(I,m,E,k)*C
      Enddo
      Enddo
         Hbar_iajk3(I,a,J,k) = Hbar_iajk3(I,a,J,k) - T
      Enddo
      Enddo
      Enddo
      Enddo

C BABA contribution

      Do K = 1, Nocc_a
      Do J = 1, Nocc_b
      Do A = 1, Nvrt_a
      Do I = 1, Nocc_b
         T = 0.0D0
         C = Ocn_ob(j)
      Do E = 1, Nvrt_a
      Do M = 1, Nocc_a
         T = T + T2aa(A,E,K,M)*W5_ba(M,i,E,j)*C
      Enddo
      Enddo
         Hbar_iajk4(i,A,j,K) = Hbar_iajk4(i,A,j,K) + T
      Enddo
      Enddo
      Enddo
      Enddo

      Do K = 1, Nocc_a
      Do J = 1, Nocc_b
      Do A = 1, Nvrt_a
      Do I = 1, Nocc_b
         T = 0.0D0
         C = Ocn_ob(j)
      Do E = 1, Nvrt_b
      Do M = 1, Nocc_b
         T = T + T2ab(A,e,K,m)*W5_bb(i,m,j,e)*C
      Enddo
      Enddo
         Hbar_iajk4(i,A,j,K) = Hbar_iajk4(i,A,j,K) + T
      Enddo
      Enddo
      Enddo
      Enddo

      Do K = 1, Nocc_a
      Do J = 1, Nocc_b
      Do A = 1, Nvrt_a
      Do I = 1, Nocc_b
         T = 0.0D0
         C = Ocn_oa(K)
      Do E = 1, Nvrt_b
      Do M = 1, Nocc_a
         T = T + T2ab(A,e,M,j)*W5_ab(M,i,K,e)*C
      Enddo
      Enddo
         Hbar_iajk4(i,A,j,K) = Hbar_iajk4(i,A,j,K) - T
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVLM
      Write(6,*) "@-term1"
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
