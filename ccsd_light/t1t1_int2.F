      Subroutine T1t1_int2(T1aa,T1bb,T2resid_aa,T2resid_bb,T2resid_ab,
     +                     W3_aaaa,W3_bbbb,W3_abab,W3_baba,W3_abba,
     +                     W3_baab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)
 
      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension T2resid_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2resid_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2resid_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension W3_aaaa(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W3_bbbb(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W3_abab(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W3_baba(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension W3_abba(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension W3_baab(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Integer A,B,I,J,M,E

#include "ccsdlight_vars.com"

#ifdef _DEBUG_LVLM
      L_aaaa = Nocc_a*Nvrt_a*Nvrt_a*Nocc_a
      L_bbbb = Nocc_b*Nvrt_b*Nvrt_b*Nocc_b
      L_abab = Nocc_a*Nvrt_b*Nvrt_a*Nocc_b
      L_baba = Nocc_b*Nvrt_a*Nvrt_b*Nocc_a
      L_abba = Nocc_a*Nvrt_b*Nocc_a*Nvrt_b
      L_baab = Nocc_b*Nvrt_a*Nocc_b*Nvrt_a
      call checksum("Form_wMBEJ:",W3_aaaa,L_aaaa)
      call checksum("Form_wmbej:",W3_bbbb,L_bbbb)
      call checksum("Form_wMbEj:",W3_abab,L_abab)
      call checksum("Form_wmBeJ:",W3_baba,L_baba)
      call checksum("Form_wmBEj:",W3_baab,L_baab)
      call checksum("Form_wMbeJ:",W3_abba,L_abba)
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("T2resid_aa:",T2resid_aa,L_aaaa)
      call checksum("T2resid_bb:",T2resid_bb,L_bbbb)
      call checksum("T2resid_ab:",T2resid_ab,L_abab)
#endif

C AAAA block 
C Res(A,B,I,J) = -P_(AB)P_(IJ)T1(A,M)*T1(E,I)*W(M,B,E,J)
      
      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a 
      Do E = 1, Nvrt_a 
CSSS         C = (1.0D0-Ocn_va(E))*Ocn_oa(M)
         C = Ocn_oa(J)*(1.0D0-Ocn_va(B))
         T = T + T1aa(E,I)*T1aa(A,M)*W3_aaaa(M,B,E,J)*C 
      Enddo 
      Enddo 

         T2resid_aa(A,B,I,J) = T2resid_aa(A,B,I,J) - T
         T2resid_aa(B,A,I,J) = T2resid_aa(B,A,I,J) + T
         T2resid_aa(A,B,J,I) = T2resid_aa(A,B,J,I) + T
         T2resid_aa(B,A,J,I) = T2resid_aa(B,A,J,I) - T

      Enddo 
      Enddo 
      Enddo 
      Enddo 

C BBBB block 
C Res(a,b,i,j) = -P_(ab)P_(ij)T1(a,m)*T1(e,i)*W(m,b,e,j)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         T  = 0.0D0
      Do M = 1, Nocc_b
      Do E = 1, Nvrt_b
CSSS         C = (1.0D0-Ocn_vb(e))*Ocn_ob(m)
         C = Ocn_ob(j)*(1.0D0-Ocn_vb(b))
         T = T + T1bb(a,m)*T1bb(e,i)*W3_bbbb(m,b,e,j)*C
      Enddo
      Enddo

         T2resid_bb(a,b,i,j) = T2resid_bb(a,b,i,j) - T
         T2resid_bb(b,a,i,j) = T2resid_bb(b,a,i,j) + T
         T2resid_bb(a,b,j,i) = T2resid_bb(a,b,j,i) + T
         T2resid_bb(b,a,j,i) = T2resid_bb(b,a,j,i) - T

      Enddo
      Enddo
      Enddo
      Enddo 

C ABAB block 
C Res(A,b,I,j) = -T1(A,M)*T1(E,I)*W(M,b,E,j)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
      Do E = 1, Nvrt_a
CSSS         C = (1.0D0-Ocn_va(E))*Ocn_oa(M)
         C = Ocn_ob(j)*(1.0D0-Ocn_vb(b))
         T = T + T1aa(A,M)*T1aa(E,I)*W3_abab(M,b,E,j)*C
      Enddo
      Enddo
         T2resid_ab(A,b,I,j) = T2resid_ab(A,b,I,j) - T 
      Enddo
      Enddo
      Enddo
      Enddo 

C Res(A,b,I,j) = +T1(A,M)*T1(e,j)*W(M,b,I,e)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
      Do E = 1, Nvrt_b
CSSS         C = (1.0D0-Ocn_vb(e))*Ocn_oa(M)
         C = Ocn_oa(I)*(1.0D0-Ocn_vb(b))
         T = T + T1aa(A,M)*T1bb(e,j)*W3_abba(M,b,I,e)*C
      Enddo
      Enddo
         T2resid_ab(A,b,I,j) = T2resid_ab(A,b,I,j) - T 
      Enddo
      Enddo
      Enddo
      Enddo

C Res(A,b,I,j) = +T1(E,I)*T1(b,m)*W(m,A,j,E)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_b
      Do E = 1, Nvrt_a
CSSS         C = (1.0D0-Ocn_va(E))*Ocn_ob(m)
         C = Ocn_ob(j)*(1.0D0-Ocn_va(A))
         T = T + T1aa(E,I)*T1bb(b,m)*W3_baab(m,A,j,E)*C
      Enddo
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
      Do E = 1, Nvrt_b
CSSS         C = (1.0D0-Ocn_vb(e))*Ocn_ob(m)
         C = Ocn_oa(I)*(1.0D0-Ocn_va(A))
         T = T + T1bb(b,m)*T1bb(e,j)*W3_baba(m,A,e,I)*C
      Enddo
      Enddo
         T2resid_ab(A,b,I,j) = T2resid_ab(A,b,I,j) - T
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " T1t1_int2"
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("T2resid_aa:",T2resid_aa,L_aaaa)
      call checksum("T2resid_bb:",T2resid_bb,L_bbbb)
      call checksum("T2resid_ab:",T2resid_ab,L_abab)
#endif

      Return
      End 
