      Subroutine Form_hbabci_term3(Hbar_abci1,Hbar_abci2,Hbar_abci3,
     +                             Hbar_abci4,Hbar_mea,Hbar_meb,
     +                             T2aa,T2bb,T2ab,Nocc_a,Nocc_b,
     +                             Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension Hbar_abci1(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Hbar_abci2(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Hbar_abci3(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Hbar_abci4(Nvrt_b,Nvrt_a,Nvrt_b,Nocc_a)

      Dimension Hbar_mea(Nocc_a,Nvrt_a)
      Dimension Hbar_meb(Nocc_b,Nvrt_b)

      Integer A,B,C,I,M,E

#include "ccsdlight_vars.com"

C AAAA contribution
      
      Do I = 1, Nocc_a
      Do C = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
         T = T + T2aa(A,B,I,M)*Hbar_mea(M,C)
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
      Do M = 1, Nocc_b
         T = T + T2bb(a,b,i,m)*Hbar_meb(m,c)
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
      Do M = 1, Nocc_a
         T = T + T2ab(A,b,M,i)*Hbar_mea(M,C)
      Enddo
         Hbar_abci3(A,b,C,i) = Hbar_abci3(A,b,C,i) - T
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
      Do M = 1, Nocc_b
         T = T + T2ab(B,a,I,m)*Hbar_meb(m,c)
      Enddo
         Hbar_abci4(a,B,c,I) = Hbar_abci4(a,B,c,I) - T
      Enddo
      Enddo
      Enddo
      Enddo

      Return
      End
