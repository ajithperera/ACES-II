      Subroutine Form_t1_2hbar_abci(Hbar_abci1,Hbar_abci2,
     +                              Hbar_abci3,Hbar_abci4,
     +                              T1aa,T1bb,W2_aa,W2_bb,
     +                              W2_ab,Nocc_a,Nocc_b,
     +                              Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)
     
      Dimension W2_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension W2_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension W2_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Dimension Hbar_abci1(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Hbar_abci2(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Hbar_abci3(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Hbar_abci4(Nvrt_b,Nvrt_a,Nvrt_b,Nocc_a)

      Integer A,B,C,I,E

#include "ccsdlight_vars.com"

C AAAA block 

      Do I = 1, Nocc_a
      Do C = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
        T = 0.0D0
        D = (1.0D0-Ocn_va(A))*(1.0D0-Ocn_va(B))
      Do E = 1, Nvrt_a 
        T = T + W2_aa(A,B,C,E)*T1aa(E,I)
      Enddo 
        Hbar_abci1(A,B,C,I) = Hbar_abci1(A,B,C,I) + T*D
      Enddo 
      Enddo 
      Enddo 
      Enddo 

C BBBB block 

      Do I = 1, Nocc_b
      Do C = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         T = 0.0D0
         D = (1.0D0-Ocn_vb(a))*(1.0D0-Ocn_vb(b))
      Do E = 1, Nvrt_b
         T = T + W2_bb(a,b,c,e)*T1bb(e,i)
      Enddo
         Hbar_abci2(a,b,c,i) = Hbar_abci2(a,b,c,i) + T*D
      Enddo
      Enddo
      Enddo
      Enddo

C ABAB block 

      Do I = 1, Nocc_b
      Do C = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do E = 1, Nvrt_b
         D = (1.0D0-Ocn_va(A))*(1.0D0-Ocn_vb(b))
         T = T + W2_ab(A,b,C,e)*T1bb(e,i)
      Enddo
      Hbar_abci3(A,b,C,i) = Hbar_abci3(A,b,C,i) + T*D
      Enddo
      Enddo
      Enddo
      Enddo

      Do I = 1, Nocc_a
      Do C = 1, Nvrt_b
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_b
         T = 0.0D0
         D = (1.0D0-Ocn_vb(a))*(1.0D0-Ocn_va(B))
      Do E = 1, Nvrt_a
         T = T + W2_ab(B,a,E,c)*T1aa(E,I)
      Enddo
      Hbar_abci4(a,B,c,I) = Hbar_abci4(a,B,c,I) + T*D
      Enddo
      Enddo
      Enddo
      Enddo

      Return
      End 
