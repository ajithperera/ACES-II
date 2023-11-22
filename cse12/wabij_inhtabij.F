      Subroutine Wabij_inhtabij(W0_aa,W0_bb,W0_ab,Htabij_aa,Htabij_bb,
     +                         Htabij_ab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension Htabij_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Htabij_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Htabij_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension W0_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension W0_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension W0_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Integer I,J,A,B

#include "cse12_vars.com"

C AAAA block 

      DO J = 1, Nocc_a
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_a
      DO A = 1, Nvrt_a
         C = Ocn_oa(I)*Ocn_oa(j)*(1.0D0-Ocn_va(A))*(1.0D0-Ocn_va(B))
         Htabij_aa(A,B,I,J) = W0_aa(A,B,I,J)*C
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C BBBB block 

      DO J = 1, Nocc_b
      DO I = 1, Nocc_b
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_b
       C = Ocn_ob(I)*Ocn_ob(j)*(1.0D0-Ocn_vb(a))*(1.0D0-Ocn_vb(b))
       Htabij_bb(a,b,i,j) = W0_bb(a,b,i,j)*C 
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C ABAB block 


      DO J = 1, Nocc_b
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_a
         C = Ocn_oa(I)*Ocn_ob(j)*(1.0D0-Ocn_va(A))*(1.0D0-Ocn_vb(b))
         Htabij_ab(A,b,I,j) = W0_ab(A,b,I,j)*C
      ENDDO
      ENDDO
      ENDDO
      ENDDO

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " @-Wabij_inhtabij"
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("Htabij_aa :",Htabij_aa,L_aaaa)
      call checksum("Htabij_bb :",Htabij_bb,L_bbbb)
      call checksum("Htabij_ab :",Htabij_ab,L_abab)
#endif

      Return
      End
