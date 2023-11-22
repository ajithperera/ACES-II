      Subroutine Form_hbiajk_term2(Hbar_iajk1,Hbar_iajk2,Hbar_iajk3,
     +                             Hbar_iajk4,W4_aa,W4_bb,W4_ab,W4_ba,
     +                             Tau_aa,Tau_bb,Tau_ab,Nocc_a,Nocc_b,
     +                             Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension Tau_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Tau_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Tau_ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension W4_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W4_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W4_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W4_ba(Nvrt_a,Nvrt_b,Nocc_a,Nvrt_b)

      Dimension Hbar_iajk1(Nocc_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Hbar_iajk2(Nocc_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Hbar_iajk3(Nocc_a,Nvrt_b,Nocc_a,Nocc_b)
      Dimension Hbar_iajk4(Nocc_b,Nvrt_a,Nocc_b,Nocc_a)

      Integer K,J,A,I,F,E

#include "ccsdlight_vars.com"

C AAAA contribution
      
      Do K = 1, Nocc_a
      Do J = 1, Nocc_a
      Do A = 1, Nvrt_a
      Do I = 1, Nocc_a
         T = 0.0D0
         C = (1.0D0 - Ocn_va(A))
      Do F = 1, Nvrt_a
      Do E = 1, Nvrt_a
         T = T + Tau_aa(E,F,J,K)*W4_aa(F,E,A,I)*C
      Enddo
      Enddo
         Hbar_iajk1(I,A,J,K) = Hbar_iajk1(I,A,J,K) + T*0.50D0
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
         C = (1.0D0 - Ocn_vb(a))
      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_b
         T = T + Tau_bb(e,f,j,k)*W4_bb(f,e,a,i)*C
      Enddo
      Enddo
         Hbar_iajk2(i,a,j,k) = Hbar_iajk2(i,a,j,k) + T*0.50D0
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
         C = (1.0D0 - Ocn_vb(a))
      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_a
         T = T + Tau_ab(E,f,J,k)*W4_ba(E,a,I,f)*C
      Enddo
      Enddo
         Hbar_iajk3(I,a,J,k) = Hbar_iajk3(I,a,J,k) + T
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
         C = (1.0D0 - Ocn_va(A))
      Do F = 1, Nvrt_a
      Do E = 1, Nvrt_b
         T = T + Tau_ab(F,e,K,j)*W4_ab(A,e,F,i)*C
      Enddo
      Enddo
         Hbar_iajk4(i,A,j,K) = Hbar_iajk4(i,A,j,K) + T
      Enddo
      Enddo
      Enddo
      Enddo

      Return
      End
