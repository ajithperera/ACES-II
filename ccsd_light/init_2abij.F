      Subroutine Init_2abij(T2resid_aa,T2resid_bb,T2resid_ab,
     +                      T1resid_aa,T1resid_bb,W0_aa,W0_bb,
     +                      W0_ab,Fockov_aa,Fockov_bb,Nocc_a,
     +                      Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)
      
      Dimension T2resid_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2resid_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2resid_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension T1resid_aa(Nvrt_a,Nocc_a)
      Dimension T1resid_bb(Nvrt_b,Nocc_b)

      Dimension W0_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension W0_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension W0_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Fockov_aa(Nocc_a,Nvrt_a) 
      Dimension Fockov_bb(Nocc_b,Nvrt_b) 

      Integer I,J,A,B

#include "ccsdlight_vars.com"

C AAAA block 

      DO J = 1, Nocc_a
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_a
      DO A = 1, Nvrt_a
         C = Ocn_oa(I)*Ocn_oa(j)*(1.0D0-Ocn_va(A))*(1.0D0-Ocn_va(B))
         T2resid_aa(A,B,I,J) = W0_aa(A,B,I,J)*C
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
         T2resid_bb(a,b,i,j) = W0_bb(a,b,i,j)*C
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
         T2resid_ab(A,b,I,j) = W0_ab(A,b,I,j)*C
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C AA block 

      DO I = 1, Nocc_a
      DO A = 1, Nvrt_a
         C = Ocn_oa(I)*(1.0D0-Ocn_va(A))
         T1resid_aa(A,I) = fockov_aa(I,A)*C
#ifdef _DEBUG_LVLM
         Write(6,"(a,2I2,1x,3(1x,F15.10))")"A,I,T1resid,C1,C2: ",
     +                                      A,I,fockov_aa(A,I),
     +                                      Ocn_oa(I),(1.0D0-Ocn_va(A))
#endif 
      ENDDO
      ENDDO

C BB block

      DO I = 1, Nocc_b
      DO A = 1, Nvrt_b
         C = Ocn_ob(i)*(1.0D0-Ocn_vb(a))
         T1resid_bb(a,i) = fockov_bb(i,a)*C 
      ENDDO
      ENDDO

#ifdef _DEBUG_LVL0
      Write(6,*) 
      L_aa   = Nocc_a*Nvrt_a
      L_bb   = Nocc_b*Nvrt_b
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("T1Resid_aa:",T1Resid_aa,L_aa)
      call checksum("T1Resid_bb:",T1Resid_bb,L_bb)
      call checksum("T2resid_aa:",T2Resid_aa,L_aaaa)
      call checksum("T2resid_bb:",T2Resid_bb,L_bbbb)
      call checksum("T2Resid_ab:",T2Resid_ab,L_abab)
      Write(6,*) 
#endif

      Return
      End
