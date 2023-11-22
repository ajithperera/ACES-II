      Subroutine Form_hbiajk_term5(Hbar_iajk1,Hbar_iajk2,Hbar_iajk3,
     +                             Hbar_iajk4,Wtt_mbej1,Wtt_mbej2,
     +                             Wtt_mbej3,Wtt_mbej4,Wtt_mbej5,
     +                             Wtt_mbej6,Taa,Tbb,Nocc_a,Nocc_b,
     +                             Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension Taa(Nvrt_a,Nocc_a)
      Dimension Tbb(Nvrt_b,Nocc_b)

      Dimension Wtt_mbej1(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Wtt_mbej2(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Wtt_mbej3(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Wtt_mbej4(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension Wtt_mbej5(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension Wtt_mbej6(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

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
      Do E = 1, Nvrt_a
         T = T + Taa(E,J)*Wtt_mbej1(I,A,E,K)
      Enddo
         Hbar_iajk1(I,A,J,K) = Hbar_iajk1(I,A,J,K) + T
         Hbar_iajk1(I,A,K,J) = Hbar_iajk1(I,A,K,J) - T
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
         T = T + Tbb(e,j)*Wtt_mbej2(i,a,e,k)
      Enddo
         Hbar_iajk2(i,a,j,k) = Hbar_iajk2(i,a,j,k) + T
         Hbar_iajk2(i,a,k,j) = Hbar_iajk2(i,a,k,j) - T
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
         T = T + Taa(E,J)*Wtt_mbej3(I,a,E,k)
      Enddo
         Hbar_iajk3(I,a,J,k) = Hbar_iajk3(I,a,J,k) + T
      Enddo
      Enddo
      Enddo
      Enddo

      Do I = 1, Nocc_a
      Do A = 1, Nvrt_b
      Do J = 1, Nocc_a
      Do K = 1, Nocc_b
         T = 0.0D0
      Do E = 1, Nvrt_b
         T = T + Tbb(e,k)*Wtt_mbej5(I,a,J,e)
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
         T = T + Tbb(e,j)*Wtt_mbej4(i,A,e,K)
      Enddo
         Hbar_iajk4(i,A,j,K) = Hbar_iajk4(i,A,j,K) + T
      Enddo 
      Enddo 
      Enddo
      Enddo 

      Do I = 1, Nocc_b
      Do A = 1, Nvrt_a
      Do J = 1, Nocc_b
      Do K = 1, Nocc_a
         T = 0.0D0
      Do E = 1, Nvrt_a
         T = T + Taa(E,K)*Wtt_mbej6(i,A,j,E)
      Enddo
         Hbar_iajk4(i,A,j,K) = Hbar_iajk4(i,A,j,K) + T
      Enddo
      Enddo
      Enddo
      Enddo


      Return
      End
