      Subroutine T1_int2(T1aa,T1bb,T2resid_aa,T2resid_bb,T2resid_ab,
     +                   W4_aa,W4_bb,W4_ab,W4_ba,W5_aa,W5_bb,
     +                   W5_ab,W5_ba,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)
 
      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension T2resid_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2resid_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2resid_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension W5_aa(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension W5_bb(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension W5_ab(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension W5_ba(Nocc_a,Nocc_b,Nvrt_a,Nocc_b)

      Dimension W4_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W4_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W4_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W4_ba(Nvrt_a,Nvrt_b,Nocc_a,Nvrt_b)

      Integer A,B,I,J,M,E

#include "ccsdlight_vars.com"

C AAAA block 
C Res(A,B,I,J) = P_(I,J)T1(E,I)*W(A,B,E,J)

      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T  = 0.0D0
         T1 = 0.0D0
         T2 = 0.0D0
      Do E = 1, Nvrt_a
C         C = (1.0D0-Ocn_va(E))
         C  = Ocn_oa(J)*(1.0D0-Ocn_va(A))*(1.0D0-Ocn_va(B))
         C1 = Ocn_oa(J)*(1.0D0-Ocn_va(A))*(1.0D0-Ocn_va(B))
         C2 = Ocn_oa(I)*(1.0D0-Ocn_va(A))*(1.0D0-Ocn_va(B))
         T  = T  + T1aa(E,I)*W4_aa(A,B,E,J)*C
         T1 = T1 + T1aa(E,I)*W4_aa(A,B,E,J)*C1
         T2 = T2 + T1aa(E,J)*W4_aa(A,B,E,I)*C2
      Enddo
         T2resid_aa(A,B,I,J) = T2resid_aa(A,B,I,J) + T1
         T2resid_aa(A,B,J,I) = T2resid_aa(A,B,J,I) - T2
      Enddo
      Enddo
      Enddo
      Enddo

C Res(A,B,I,J) = -P_(AB)T1(A,M)*W(I,J,M,B)

      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T  = 0.0D0
         T1 = 0.0D0
         T2 = 0.0D0
      Do M = 1, Nocc_a
CSSS         C = Ocn_oa(M)
         C  = Ocn_oa(I)*Ocn_oa(J)*(1.0D0-Ocn_va(B))
         C1 = Ocn_oa(I)*Ocn_oa(J)*(1.0D0-Ocn_va(B))
         C2 = Ocn_oa(I)*Ocn_oa(J)*(1.0D0-Ocn_va(A))
         T  = T  + T1aa(A,M)*W5_aa(I,J,M,B)*C
         T1 = T1 + T1aa(A,M)*W5_aa(I,J,M,B)*C1
         T2 = T2 + T1aa(B,M)*W5_aa(I,J,M,A)*C2
      Enddo
         T2resid_aa(A,B,I,J) = T2resid_aa(A,B,I,J) - T1
         T2resid_aa(B,A,I,J) = T2resid_aa(B,A,I,J) + T2
      Enddo
      Enddo
      Enddo
      Enddo

C BBBB block 
C Res(a,b,i,j) = P_(i,j)T1(e,i)*W(a,b,e,j)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         T   = 0.0D0
         T1  = 0.0D0
         T2  = 0.0D0
      Do E = 1, Nvrt_b
CSSS         C = (1.0D0-Ocn_vb(e))
         C  = Ocn_ob(j)*(1.0D0-Ocn_vb(a))*(1.0D0-Ocn_vb(b))
         C1 = Ocn_ob(j)*(1.0D0-Ocn_vb(a))*(1.0D0-Ocn_vb(b))
         C2 = Ocn_ob(i)*(1.0D0-Ocn_vb(a))*(1.0D0-Ocn_vb(b))
         T  = T  + T1bb(e,i)*W4_bb(a,b,e,j)*C
         T1 = T1 + T1bb(e,i)*W4_bb(a,b,e,j)*C1
         T2 = T2 + T1bb(e,j)*W4_bb(a,b,e,i)*C2
      Enddo
         T2resid_bb(a,b,i,j) = T2resid_bb(a,b,i,j) + T1
         T2resid_bb(a,b,j,i) = T2resid_bb(a,b,j,i) - T2
      Enddo 
      Enddo 
      Enddo 
      Enddo 

C Res(a,b,i,j) = -P_(ab)T1(a,m)*W(i,j,m,b)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         T   = 0.0D0
         T1  = 0.0D0
         T2  = 0.0D0
      Do M = 1, Nocc_b
CSSS         C = Ocn_ob(m)
         C  = Ocn_ob(i)*Ocn_ob(j)*(1.0D0-Ocn_vb(b))
         C1 = Ocn_ob(i)*Ocn_ob(j)*(1.0D0-Ocn_vb(b))
         C2 = Ocn_ob(i)*Ocn_ob(j)*(1.0D0-Ocn_vb(a))
         T  = T  + T1bb(a,m)*W5_bb(i,j,m,b)*C
         T1 = T1 + T1bb(a,m)*W5_bb(i,j,m,b)*C1
         T2 = T2 + T1bb(b,m)*W5_bb(i,j,m,a)*C2
      Enddo 
         T2resid_bb(a,b,i,j) = T2resid_bb(a,b,i,j) - T1
         T2resid_bb(b,a,i,j) = T2resid_bb(b,a,i,j) + T2
      Enddo 
      Enddo 
      Enddo 
      Enddo 
      
C ABAB block 
C Res(A,b,I,j) = T1(E,I)*W(A,b,E,j) - T(e,j)*W(A,b,e,I)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do E = 1, Nvrt_a
CSSS         C = (1.0D0-Ocn_va(E))
         C = Ocn_ob(j)*(1.0D0-Ocn_va(A))*(1.0D0-Ocn_vb(b))
         T = T + T1aa(E,I)*W4_ab(A,b,E,j)*C
      Enddo
         T2resid_ab(A,b,I,j) = T2resid_ab(A,b,I,j) + T 
      Enddo
      Enddo
      Enddo
      Enddo 

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do E = 1, Nvrt_b
CSSS         C = (1.0D0-Ocn_vb(e))
         C = Ocn_oa(I)*(1.0D0-Ocn_va(A))*(1.0D0-Ocn_vb(b))
         T = T + T1bb(e,j)*W4_ba(A,b,I,e)*C
      Enddo
         T2resid_ab(A,b,I,j) = T2resid_ab(A,b,I,j) + T
      Enddo
      Enddo
      Enddo
      Enddo

C Res(A,b,I,j) = -T1(A,M)*W(I,j,M,b) + T1(b,m)*W(I,j,m,A)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
C         C = Ocn_oa(M)
         C = Ocn_oa(I)*Ocn_ob(j)*(1.0D0-Ocn_vb(b))
         T = T + T1aa(A,M)*W5_ab(I,j,M,b)*C
      Enddo
         T2resid_ab(A,b,I,j) = T2resid_ab(A,b,I,j) - T
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_b
CSSS         C = Ocn_ob(m)
         C =  Ocn_oa(I)*Ocn_ob(j)*(1.0D0-Ocn_va(A))
         T = T + T1bb(b,m)*W5_ba(I,j,A,m)*C
      Enddo
         T2resid_ab(A,b,I,j) = T2resid_ab(A,b,I,j) - T
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " T1_in_t2"
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("T2resid_aa:",T2resid_aa,L_aaaa)
      call checksum("T2resid_bb:",T2resid_bb,L_bbbb)
      call checksum("T2resid_ab:",T2resid_ab,L_abab)
#endif

      Return
      End 
