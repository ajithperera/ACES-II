      Subroutine T1_inhtabij(T1aa,T1bb,Htabij_aa,Htabij_bb,Htabij_ab,
     +                       W4_aa,W4_bb,W4_ab,W4_ba,W5_aa,W5_bb,
     +                       W5_ab,W5_ba,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)
 
      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension Htabij_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Htabij_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Htabij_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension W5_aa(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension W5_bb(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension W5_ab(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension W5_ba(Nocc_a,Nocc_b,Nvrt_a,Nocc_b)

      Dimension W4_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W4_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W4_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W4_ba(Nvrt_a,Nvrt_b,Nocc_a,Nvrt_b)

      Integer A,B,I,J,M,E

#include "cse12_vars.com"

C AAAA block 
C Htabij_aa(A,B,I,J) <- P_(I,J)T1(E,I)*W(A,B,E,J)

      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do E = 1, Nvrt_a
         C = Ocn_oa(J)*(1.0D0-Ocn_va(A))*(1.0D0-Ocn_va(B))
         T = T + T1aa(E,I)*W4_aa(A,B,E,J)*C
      Enddo
         Htabij_aa(A,B,I,J) = Htabij_aa(A,B,I,J) + T
         Htabij_aa(A,B,J,I) = Htabij_aa(A,B,J,I) - T
      Enddo
      Enddo
      Enddo
      Enddo

C Htabij_aa(A,B,I,J) <- -P_(AB)T1(A,M)*W(I,J,M,B)

      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
         C = Ocn_oa(I)*Ocn_oa(J)*(1.0D0-Ocn_va(B))
         T = T + T1aa(A,M)*W5_aa(I,J,M,B)*C
      Enddo
         Htabij_aa(A,B,I,J) = Htabij_aa(A,B,I,J) - T
         Htabij_aa(B,A,I,J) = Htabij_aa(B,A,I,J) + T
      Enddo
      Enddo
      Enddo
      Enddo

C BBBB block 
C Htabij_bb(a,b,i,j) <- P_(i,j)T1(e,i)*W(a,b,e,j)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         T  = 0.0D0
      Do E = 1, Nvrt_b
         C = Ocn_ob(j)*(1.0D0-Ocn_vb(a))*(1.0D0-Ocn_vb(b))
         T = T + T1bb(e,i)*W4_bb(a,b,e,j)*C
      Enddo
         Htabij_bb(a,b,i,j) = Htabij_bb(a,b,i,j) + T
         Htabij_bb(a,b,j,i) = Htabij_bb(a,b,j,i) - T
      Enddo 
      Enddo 
      Enddo 
      Enddo 

C Htabij_bb(a,b,i,j) <- -P_(ab)T1(a,m)*W(i,j,m,b)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         T  = 0.0D0
      Do M = 1, Nocc_b
         C = Ocn_ob(i)*Ocn_ob(j)*(1.0D0-Ocn_vb(b))
         T = T + T1bb(a,m)*W5_bb(i,j,m,b)*C
      Enddo 
         Htabij_bb(a,b,i,j) = Htabij_bb(a,b,i,j) - T
         Htabij_bb(b,a,i,j) = Htabij_bb(b,a,i,j) + T
      Enddo 
      Enddo 
      Enddo 
      Enddo 
      
C ABAB block 
C Htabij_ab(A,b,I,j) <- T1(E,I)*W(A,b,E,j) - T(e,j)*W(A,b,e,I)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do E = 1, Nvrt_a
         C = Ocn_ob(j)*(1.0D0-Ocn_va(A))*(1.0D0-Ocn_vb(b))
         T = T + T1aa(E,I)*W4_ab(A,b,E,j)*C
      Enddo
         Htabij_ab(A,b,I,j) = Htabij_ab(A,b,I,j) + T 
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
         C = Ocn_oa(I)*(1.0D0-Ocn_va(A))*(1.0D0-Ocn_vb(b))
         T = T + T1bb(e,j)*W4_ba(A,b,I,e)*C
      Enddo
         Htabij_ab(A,b,I,j) = Htabij_ab(A,b,I,j) + T
      Enddo
      Enddo
      Enddo
      Enddo

C Htabij_ab(A,b,I,j) <- -T1(A,M)*W(I,j,M,b) + T1(b,m)*W(I,j,m,A)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
         C = Ocn_oa(I)*Ocn_ob(j)*(1.0D0-Ocn_vb(b))
         T = T + T1aa(A,M)*W5_ab(I,j,M,b)*C
      Enddo
         Htabij_ab(A,b,I,j) = Htabij_ab(A,b,I,j) - T
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
         C =  Ocn_oa(I)*Ocn_ob(j)*(1.0D0-Ocn_va(A))
         T = T + T1bb(b,m)*W5_ba(I,j,A,m)*C
      Enddo
         Htabij_ab(A,b,I,j) = Htabij_ab(A,b,I,j) - T
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " @-T1_inhtabij"
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("Htabij_aa :",Htabij_aa,L_aaaa)
      call checksum("Htabij_bb :",Htabij_bb,L_bbbb)
      call checksum("Htabij_ab :",Htabij_ab,L_abab)
#endif

      Return
      End 
