      Subroutine  Form_hbar_ijab(Hbar_ijab1,Hbar_ijab2,Hbar_ijab3,
     +                           W_aa,W_bb,W_ab,Nocc_a,Nocc_b,
     +                           Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension Hbar_ijab1(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension Hbar_ijab2(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension Hbar_ijab3(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_b)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Integer I,J,A,B

#include "ccsdlight_vars.com"

C AAAA block 

      DO J = 1, Nocc_a
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_a
      DO A = 1, Nvrt_a
CSSS         C = Ocn_oa(I)*Ocn_oa(j)*(1.0D0-Ocn_va(A))*(1.0D0-Ocn_va(B))
         Hbar_ijab1(I,J,A,B) = W_aa(I,J,A,B)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C BBBB block 

      DO J = 1, Nocc_b
      DO I = 1, Nocc_b
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_b
CSSS         C = Ocn_ob(I)*Ocn_ob(j)*(1.0D0-Ocn_vb(a))*(1.0D0-Ocn_vb(b))
         Hbar_ijab2(i,j,a,b) = W_bb(i,j,a,b)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C ABAB block 

      DO J = 1, Nocc_b
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_a
CSSS         C = Ocn_oa(I)*Ocn_ob(j)*(1.0D0-Ocn_va(A))*(1.0D0-Ocn_vb(b))
         Hbar_ijab3(I,j,A,b) = W_ab(I,j,A,b)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      Return
      End
