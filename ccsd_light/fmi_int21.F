      Subroutine Fmi_int21(T2aa_old,T2bb_old,T2ab_old,T2resid_aa,
     +                     T2resid_bb,T2resid_ab,Fmi_a,Fmi_b,
     +                     Nocc_a,Nocc_b,Nvrt_a,Nvrt_b) 

      Implicit Double Precision(A-H,O-Z)

      Dimension T2aa_old(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb_old(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab_old(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)
 
      Dimension T2resid_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2resid_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2resid_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Fmi_a(Nocc_a,Nocc_a)
      Dimension Fmi_b(Nocc_b,Nocc_b)

      Integer A,B,I,J,M

#include "ccsdlight_vars.com"

#ifdef _DEBUG_LVLM
      Write(6,*)
      call checksum("Fmi_a     :",Fmi_a,Nocc_a*Nocc_a)
      call checksum("Fmi_b     :",Fmi_b,Nocc_b*Nocc_b)
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("T2aa      :",T2aa_old,L_aaaa)
      call checksum("T2bb      :",T2bb_old,L_bbbb)
      call checksum("T2ab      :",T2ab_old,L_abab)
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("T2resid_aa:",T2resid_aa,L_aaaa)
      call checksum("T2resid_bb:",T2resid_bb,L_bbbb)
      call checksum("T2resid_ab:",T2resid_ab,L_abab)
#endif

C AAAA block 

      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a 
CSSS         C = Ocn_oa(M) 
         C = 1.0D0
         T = T  + Fmi_a(M,J)*T2aa_old(A,B,I,M)*C
      Enddo 
         T2resid_aa(A,B,I,J) = T2resid_aa(A,B,I,J) - T 
C         T2resid_aa(A,B,J,I) = T2resid_aa(A,B,J,I) + T 
      Enddo 
      Enddo 
      Enddo 
      Enddo 

C BBBB block 

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         T = 0.0D0
      Do M = 1, Nocc_b
CSSS         C = Ocn_ob(m) 
         C = 1.0D0
         T = T + Fmi_b(m,j)*T2bb_old(a,b,i,m)*C
      Enddo
         T2resid_bb(a,b,i,j) = T2resid_bb(a,b,i,j) - T
C         T2resid_bb(a,b,j,i) = T2resid_bb(a,b,j,i) + T
      Enddo
      Enddo
      Enddo
      Enddo

C ABAB block 

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_b
CSSS         C = Ocn_ob(m) 
         C = 1.0D0
         T = T + Fmi_b(m,j)*T2ab_old(A,b,I,m)*C
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
      Do M = 1, Nocc_a
CSSS         C = Ocn_oa(M) 
         C = 1.0D0
         T = T + Fmi_a(M,I)*T2ab_old(A,b,M,j)*C 
      Enddo
C         T2resid_ab(A,b,I,j) = T2resid_ab(A,b,I,j)  - T
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " Fmi_int2"
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("T2resid_aa:",T2resid_aa,L_aaaa)
      call checksum("T2resid_bb:",T2resid_bb,L_bbbb)
      call checksum("T2resid_ab:",T2resid_ab,L_abab)
#endif
      Return
      End 
