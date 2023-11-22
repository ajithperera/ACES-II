      Subroutine Form_hbabci_term4(Hbar_abci1,Hbar_abci2,Hbar_abci3,
     +                             Hbar_abci4,Hbar_abef1,Hbar_abef2,
     +                             Hbar_abef3,Taa,Tbb,Nocc_a,Nocc_b,
     +                             Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension Taa(Nvrt_a,Nocc_a)
      Dimension Tbb(Nvrt_b,Nocc_b)

      Dimension Hbar_abef1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension Hbar_abef2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension Hbar_abef3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Dimension Hbar_abci1(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Hbar_abci2(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Hbar_abci3(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Hbar_abci4(Nvrt_b,Nvrt_a,Nvrt_b,Nocc_a)

      Integer A,B,C,I,E

#include "ccsdlight_vars.com"

C AAAA contribution
      
      Do I = 1, Nocc_a
      Do C = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do E = 1, Nvrt_a
         T = T + Taa(E,I)*Hbar_abef1(A,B,C,E)
      Enddo
         Hbar_abci1(A,B,C,I) = Hbar_abci1(A,B,C,I) + T
      Enddo
      Enddo
      Enddo
      Enddo

C BBBB contribution

      Do I = 1, Nocc_b
      Do C = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         T = 0.0D0
      Do E = 1, Nvrt_b
         T = T + Tbb(e,i)*Hbar_abef2(a,b,c,e)
      Enddo
         Hbar_abci2(a,b,c,i) = Hbar_abci2(a,b,c,i) + T
      Enddo
      Enddo
      Enddo
      Enddo

C ABAB contribution

      Do I = 1, Nocc_b
      Do C = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do E = 1, Nvrt_b
         T = T + Tbb(e,i)*Hbar_abef3(A,b,C,e)
      Enddo
         Hbar_abci3(A,b,C,i) = Hbar_abci3(A,b,C,i) + T
      Enddo
      Enddo
      Enddo
      Enddo

C BABA contribution

      Do I = 1, Nocc_a
      Do C = 1, Nvrt_b
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_b
         T = 0.0D0
      Do E = 1, Nvrt_a
         T = T + Taa(E,I)*Hbar_abef3(B,a,E,c)
      Enddo
         Hbar_abci4(a,B,c,I) = Hbar_abci4(a,B,c,I) + T
      Enddo
      Enddo
      Enddo
      Enddo

      Return
      End
