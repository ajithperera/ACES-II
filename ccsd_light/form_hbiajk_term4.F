      Subroutine Form_hbiajk_term4(Hbar_iajk1,Hbar_iajk2,Hbar_iajk3,
     +                             Hbar_iajk4,Hbar_ijkl1,Hbar_ijkl2,
     +                             Hbar_ijkl3,Hbar_ijkl4,Taa,Tbb,
     +                             Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension Taa(Nvrt_a,Nocc_a)
      Dimension Tbb(Nvrt_b,Nocc_b)

      Dimension Hbar_ijkl1(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension Hbar_ijkl2(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension Hbar_ijkl3(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Dimension Hbar_iajk1(Nocc_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Hbar_iajk2(Nocc_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Hbar_iajk3(Nocc_a,Nvrt_b,Nocc_a,Nocc_b)
      Dimension Hbar_iajk4(Nocc_b,Nvrt_a,Nocc_b,Nocc_a)

      Integer A,B,C,I,M,N

#include "ccsdlight_vars.com"

C AAAA contribution
      
      Do I = 1, Nocc_a
      Do A = 1, Nvrt_a
      Do J = 1, Nocc_a
      Do K = 1, Nocc_a
         T = 0.0D0
      Do M = 1, Nocc_a
         T = T + Taa(A,M)*Hbar_ijkl1(M,I,K,J)
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
      Do M = 1, Nocc_b
         T = T + Tbb(a,m)*Hbar_ijkl2(m,i,k,j)
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
      Do M = 1, Nocc_b
         T = T + Tbb(a,m)*Hbar_ijkl3(I,m,J,k)
      Enddo
         Hbar_iajk3(I,a,J,k) = Hbar_iajk3(I,a,J,k) - T
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
      Do M = 1, Nocc_a
         T = T + Taa(A,M)*Hbar_ijkl3(M,i,K,j)
      Enddo
         Hbar_iajk4(i,A,j,K) = Hbar_iajk4(i,A,j,K) - T
      Enddo
      Enddo
      Enddo
      Enddo

      Return
      End
