      Subroutine T1_int2_modf(T1aa,T1bb,T2resid_aa,T2resid_bb,
     +                        T2resid_ab,W4_aa,W4_bb,W4_ab,
     +                        W4_ba,W5_aa,W5_bb,W5_ab,W5_ba,
     +                        Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)
 
      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension T2resid_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2resid_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2resid_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension W5_aa(Nocc_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension W5_bb(Nocc_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension W5_ab(Nocc_a,Nvrt_b,Nocc_a,Nocc_b)
      Dimension W5_ba(Nocc_b,Nvrt_a,Nocc_b,Nocc_a)

      Dimension W4_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W4_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W4_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W4_ba(Nvrt_b,Nvrt_a,Nvrt_b,Nocc_a)

      Integer A,B,I,J,M,E

#include "ccsdlight_vars.com"

C AAAA block 
C Res(A,B,I,J) = P_(I,J)T1(E,I)*W(A,B,E,J)

      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do E = 1, Nvrt_a
         T = T + T1aa(E,I)*W4_aa(A,B,E,J)
      Enddo
         T2resid_aa(A,B,I,J) = T2resid_aa(A,B,I,J) + T
         T2resid_aa(A,B,J,I) = T2resid_aa(A,B,J,I) - T
      Enddo
      Enddo
      Enddo
      Enddo

C Res(A,B,I,J) = -P_(AB)T1(A,M)*W(M,B,I,J)

      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
         T = T + T1aa(A,M)*W5_aa(M,B,I,J)
      Enddo
         T2resid_aa(A,B,I,J) = T2resid_aa(A,B,I,J) - T
         T2resid_aa(B,A,I,J) = T2resid_aa(B,A,I,J) + T
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
         T  = 0.0D0
      Do E = 1, Nvrt_b
         T = T + T1bb(e,i)*W4_bb(a,b,e,j)
      Enddo
         T2resid_bb(a,b,i,j) = T2resid_bb(a,b,i,j) + T
         T2resid_bb(a,b,j,i) = T2resid_bb(a,b,j,i) - T
      Enddo 
      Enddo 
      Enddo 
      Enddo 

C Res(a,b,i,j) = -P_(ab)T1(a,m)*W(m,b,i,j)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         T  = 0.0D0
      Do M = 1, Nocc_b
         T = T + T1bb(a,m)*W5_bb(m,b,i,j)
      Enddo 
         T2resid_bb(a,b,i,j) = T2resid_bb(a,b,i,j) - T
         T2resid_bb(b,a,i,j) = T2resid_bb(b,a,i,j) + T
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
         T = T + T1aa(E,I)*W4_ab(A,b,E,j)
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
         T = T + T1bb(e,j)*W4_ba(b,A,e,I)
      Enddo
         T2resid_ab(A,b,I,j) = T2resid_ab(A,b,I,j) + T
      Enddo
      Enddo
      Enddo
      Enddo

C Res(A,b,I,j) = -T1(A,M)*W(M,b,I,j) + T1(b,m)*W(m,A,I,j)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
         T = T + T1aa(A,M)*W5_ab(M,b,I,j)
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
         T = T + T1bb(b,m)*W5_ba(m,A,j,I)
      Enddo
         T2resid_ab(A,b,I,j) = T2resid_ab(A,b,I,j) - T
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("T2resid_aa:",T2resid_aa,L_aaaa)
      call checksum("T2resid_bb:",T2resid_bb,L_bbbb)
      call checksum("T2resid_ab:",T2resid_ab,L_abab)
#endif

      Return
      End 
