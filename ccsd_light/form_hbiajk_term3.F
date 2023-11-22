      Subroutine Form_hbiajk_term3(Hbar_iajk1,Hbar_iajk2,Hbar_iajk3,
     +                             Hbar_iajk4,Hbar_mea,Hbar_meb,
     +                             T2aa,T2bb,T2ab,Nocc_a,Nocc_b,
     +                             Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension Hbar_iajk1(Nocc_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Hbar_iajk2(Nocc_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Hbar_iajk3(Nocc_a,Nvrt_b,Nocc_a,Nocc_b)
      Dimension Hbar_iajk4(Nocc_b,Nvrt_a,Nocc_b,Nocc_a)

      Dimension Hbar_mea(Nocc_a,Nvrt_a)
      Dimension Hbar_meb(Nocc_b,Nvrt_b)

      Integer I,A,J,K,E

#include "ccsdlight_vars.com"

C AAAA contribution
      
      Do I = 1, Nocc_a
      Do A = 1, Nvrt_a
      Do J = 1, Nocc_a
      Do K = 1, Nocc_a
         T = 0.0D0
      Do E = 1, Nvrt_a
         T = T + T2aa(A,E,J,K)*Hbar_mea(I,E)
      Enddo
         Hbar_iajk1(I,A,J,K) = Hbar_iajk1(I,A,J,K) - T
      Enddo
      Enddo
      Enddo
      Enddo

C BBBB contribution

      Do I = 1, Nocc_b
      Do A = 1, Nvrt_b
      Do J = 1, Nocc_b
      Do K = 1, Nocc_b
         T = 0.0D0
      Do E = 1, Nvrt_b
         T = T + T2bb(a,e,j,k)*Hbar_meb(i,e)
      Enddo
         Hbar_iajk2(i,a,j,k) = Hbar_iajk2(i,a,j,k) - T
      Enddo
      Enddo
      Enddo
      Enddo

C ABAB contribution

      Do I = 1, Nocc_a
      Do A = 1, Nvrt_b
      Do J = 1, Nocc_a
      Do K = 1, Nocc_b
         T = 0.0D0
      Do E = 1, Nvrt_a
         T = T + T2ab(E,a,J,k)*Hbar_mea(I,E)
      Enddo
         Hbar_iajk3(I,a,J,k) = Hbar_iajk3(I,a,J,k) + T
      Enddo
      Enddo
      Enddo
      Enddo

C BABA contribution

      Do I = 1, Nocc_b
      Do A = 1, Nvrt_a
      Do J = 1, Nocc_b
      Do K = 1, Nocc_a
         T = 0.0D0
      Do E = 1, Nvrt_b
         T = T + T2ab(A,e,K,j)*Hbar_meb(i,e)
      Enddo
         Hbar_iajk4(i,A,j,K) = Hbar_iajk4(i,A,j,K) + T
      Enddo
      Enddo
      Enddo
      Enddo

      Return
      End
