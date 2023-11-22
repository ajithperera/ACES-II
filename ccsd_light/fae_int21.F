      Subroutine Fae_int21(T2aa_old,T2bb_old,T2ab_old,T2resid_aa,
     +                    T2resid_bb,T2resid_ab,Fae_a,Fae_b,Nocc_a,
     +                    Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension T2aa_old(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb_old(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab_old(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)
 
      Dimension T2resid_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2resid_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2resid_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Fae_a(Nvrt_a,Nvrt_a)
      Dimension Fae_b(Nvrt_b,Nvrt_b)

      Integer I,J,A,B,E

#include "ccsdlight_vars.com"

#ifdef _DEBUG_LVL0
      Write(6,*)
      call checksum("Fae_a     :",Fae_a,Nvrt_a*Nvrt_a)
      call checksum("Fae_b     :",Fae_b,Nvrt_b*Nvrt_b)
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
      Do E = 1, Nvrt_a 
CSSS         C = (1.0D0-Ocn_va(E))
         C = 1.0D0
         T = T + Fae_a(A,E)*T2aa_old(E,B,I,J)*C
      Enddo 
         T2resid_aa(A,B,I,J) = T2resid_aa(A,B,I,J) + T 
C         T2resid_aa(B,A,I,J) = T2resid_aa(B,A,I,J) - T 
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
      Do E = 1, Nvrt_b
CSSS         C = (1.0D0-Ocn_vb(e))
         C = 1.0D0
         T = T + Fae_b(a,e)*T2bb_old(e,b,i,j)*C
      Enddo
         T2resid_bb(a,b,i,j) = T2resid_bb(a,b,i,j) + T 
C         T2resid_bb(b,a,i,j) = T2resid_bb(b,a,i,j) - T 
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
      Do E = 1, Nvrt_a
CSSS         C = (1.0D0-Ocn_va(E))
         C = 1.0D0
         T = T + Fae_a(A,E)*T2ab_old(E,b,I,j)*C 
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
         C = 1.0D0
         T = T + Fae_b(b,e)*T2ab_old(A,e,I,j)*C 
      Enddo
C         T2resid_ab(A,b,I,j) = T2resid_ab(A,b,I,j) + T 
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " Fae_int2"
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("T2resid_aa:",T2resid_aa,L_aaaa)
      call checksum("T2resid_bb:",T2resid_bb,L_bbbb)
      call checksum("T2resid_ab:",T2resid_ab,L_abab)
#endif

      Return
      End 
