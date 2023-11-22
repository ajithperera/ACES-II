      Subroutine T1t1_inhtabij(T1aa,T1bb,Htabij_aa,Htabij_bb,Htabij_ab,
     +                         W3_aaaa,W3_bbbb,W3_abab,W3_baba,W3_abba,
     +                         W3_baab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)
 
      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension Htabij_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Htabij_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Htabij_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension W3_aaaa(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W3_bbbb(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W3_abab(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W3_baba(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension W3_abba(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension W3_baab(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Integer A,B,I,J,M,E

#include "cse12_vars.com"

C AAAA block 
C Htabij_aa(A,B,I,J) <- -P_(AB)P_(IJ)T1(A,M)*T1(E,I)*W(M,B,E,J)
      
      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a 
      Do E = 1, Nvrt_a 
         C = Ocn_oa(J)*(1.0D0-Ocn_va(B))
         T = T + T1aa(E,I)*T1aa(A,M)*W3_aaaa(M,B,E,J)*C 
      Enddo 
      Enddo 

         Htabij_aa(A,B,I,J) = Htabij_aa(A,B,I,J) - T
         Htabij_aa(B,A,I,J) = Htabij_aa(B,A,I,J) + T
         Htabij_aa(A,B,J,I) = Htabij_aa(A,B,J,I) + T
         Htabij_aa(B,A,J,I) = Htabij_aa(B,A,J,I) - T

      Enddo 
      Enddo 
      Enddo 
      Enddo 

C BBBB block 
C Htijab_bb(a,b,i,j) <- -P_(ab)P_(ij)T1(a,m)*T1(e,i)*W(m,b,e,j)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         T  = 0.0D0
      Do M = 1, Nocc_b
      Do E = 1, Nvrt_b
         C = Ocn_ob(j)*(1.0D0-Ocn_vb(b))
         T = T + T1bb(a,m)*T1bb(e,i)*W3_bbbb(m,b,e,j)*C
      Enddo
      Enddo

         Htabij_bb(a,b,i,j) = Htabij_bb(a,b,i,j) - T
         Htabij_bb(b,a,i,j) = Htabij_bb(b,a,i,j) + T
         Htabij_bb(a,b,j,i) = Htabij_bb(a,b,j,i) + T
         Htabij_bb(b,a,j,i) = Htabij_bb(b,a,j,i) - T

      Enddo
      Enddo
      Enddo
      Enddo 

C ABAB block 
C Htabij_ab(A,b,I,j) <- -T1(A,M)*T1(E,I)*W(M,b,E,j)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
      Do E = 1, Nvrt_a
         C = Ocn_ob(j)*(1.0D0-Ocn_vb(b))
         T = T + T1aa(A,M)*T1aa(E,I)*W3_abab(M,b,E,j)*C
      Enddo
      Enddo
         Htabij_ab(A,b,I,j) = Htabij_ab(A,b,I,j) - T 
      Enddo
      Enddo
      Enddo
      Enddo 

C Htabij_ab(A,b,I,j) <- +T1(A,M)*T1(e,j)*W(M,b,I,e)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
      Do E = 1, Nvrt_b
         C = Ocn_oa(I)*(1.0D0-Ocn_vb(b))
         T = T + T1aa(A,M)*T1bb(e,j)*W3_abba(M,b,I,e)*C
      Enddo
      Enddo
         Htabij_ab(A,b,I,j) = Htabij_ab(A,b,I,j) - T 
      Enddo
      Enddo
      Enddo
      Enddo

C Htabij_ab(A,b,I,j) <- +T1(E,I)*T1(b,m)*W(m,A,j,E)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_b
      Do E = 1, Nvrt_a
         C = Ocn_ob(j)*(1.0D0-Ocn_va(A))
         T = T + T1aa(E,I)*T1bb(b,m)*W3_baab(m,A,j,E)*C
      Enddo
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
      Do E = 1, Nvrt_b
         C = Ocn_oa(I)*(1.0D0-Ocn_va(A))
         T = T + T1bb(b,m)*T1bb(e,j)*W3_baba(m,A,e,I)*C
      Enddo
      Enddo
         Htabij_ab(A,b,I,j) = Htabij_ab(A,b,I,j) - T
      Enddo
      Enddo
      Enddo
      Enddo
#ifdef _DEBUG_LVL0
      Write(6,"(a)") " @-T1T1_inhtabij"
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("Htabij_aa :",Htabij_aa,L_aaaa)
      call checksum("Htabij_bb :",Htabij_bb,L_bbbb)
      call checksum("Htabij_ab :",Htabij_ab,L_abab)
#endif 

      Return
      End 
