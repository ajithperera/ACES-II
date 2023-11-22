      Subroutine Form_hbabci_term2(Hbar_abci1,Hbar_abci2,Hbar_abci3,
     +                             Hbar_abci4,Hbar_ijka1,Hbar_ijka2,
     +                             Hbar_ijka3,Hbar_ijka4,Tau_aa,
     +                             Tau_bb,Tau_ab,Nocc_a,Nocc_b,Nvrt_a,
     +                             Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension Tau_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Tau_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Tau_ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension Hbar_abci1(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Hbar_abci2(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Hbar_abci3(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Hbar_abci4(Nvrt_b,Nvrt_a,Nvrt_b,Nocc_a)

      Dimension Hbar_ijka1(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension Hbar_ijka2(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension Hbar_ijka3(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension Hbar_ijka4(Nocc_b,Nocc_a,Nocc_b,Nvrt_a)

      Integer A,B,C,I,M,N

#include "ccsdlight_vars.com"

C AAAA contribution
      
      Do I = 1, Nocc_a
      Do C = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do N = 1, Nocc_a
      Do M = 1, Nocc_a
         T = T + Tau_aa(A,B,M,N)*Hbar_ijka1(N,M,I,C)
      Enddo
      Enddo
         Hbar_abci1(A,B,C,I) = Hbar_abci1(A,B,C,I) + T*0.50D0
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
      Do N = 1, Nocc_b
      Do M = 1, Nocc_b
         T = T + Tau_bb(a,b,m,n)*Hbar_ijka2(n,m,i,c)
      Enddo
      Enddo
         Hbar_abci2(a,b,c,i) = Hbar_abci2(a,b,c,i) + T*0.50D0
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
      Do N = 1, Nocc_b
      Do M = 1, Nocc_a
         T = T + Tau_ab(A,b,M,n)*Hbar_ijka4(n,M,i,C)
      Enddo
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
      Do N = 1, Nocc_a
      Do M = 1, Nocc_b
         T = T + Tau_ab(B,a,N,m)*Hbar_ijka3(N,m,I,c)
      Enddo
      Enddo
         Hbar_abci4(a,B,c,I) = Hbar_abci4(a,B,c,I) + T
      Enddo
      Enddo
      Enddo
      Enddo

      Return
      End
