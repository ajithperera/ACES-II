      Subroutine Wabef_int2(T2aa_old,T2bb_old,T2ab_old,T2resid_aa,
     +                      T2resid_bb,T2resid_ab,Wabef_1,Wabef_2,
     +                      Wabef_3,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)
     
      Dimension T2aa_old(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb_old(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab_old(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)
 
      Dimension T2resid_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2resid_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2resid_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Wabef_1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension Wabef_2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension Wabef_3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Integer I,J,A,B,E,F

#include "ccsdlight_vars.com"

#ifdef _DEBUG_LVLM
      Write(6,*)
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
      Lenabef_1 = Nvrt_a*Nvrt_a*Nvrt_a*Nvrt_a
      Lenabef_2 = Nvrt_b*Nvrt_b*Nvrt_b*Nvrt_b
      Lenabef_3 = Nvrt_a*Nvrt_b*Nvrt_a*Nvrt_b
      call checksum("Hbar_abef1:",Wabef_1,Lenabef_1)
      call checksum("Hbar_abef2:",Wabef_2,Lenabef_2)
      call checksum("Hbar_abef3:",Wabef_3,Lenabef_3)
#endif

C AAAA block 

      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
        T = 0.0D0
      Do E = 1, Nvrt_a 
      Do F = 1, Nvrt_a 
CSSS        C = (1.0D0-Ocn_va(E))*(1.0D0-Ocn_va(F))
        C = 1.0D0
        T = T + 0.50D0*Wabef_1(A,B,E,F)*T2aa_old(E,F,I,J)*C 
      Enddo 
      Enddo 
         T2resid_aa(A,B,I,J) = T2resid_aa(A,B,I,J) + T 
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
      Do F = 1, Nvrt_b
CSSS         C = (1.0D0-Ocn_vb(e))*(1.0D0-Ocn_vb(f))
         C = 1.0D0
         T = T + 0.50D0*Wabef_2(a,b,e,f)*T2bb_old(e,f,i,j)*C 
      Enddo
      Enddo
         T2resid_bb(a,b,i,j) = T2resid_bb(a,b,i,j) + T 
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
      Do F = 1, Nvrt_b
CSSS         C = (1.0D0-Ocn_va(E))*(1.0D0-Ocn_vb(f))
         C = 1.0D0 
         T = T + Wabef_3(A,b,E,f)*T2ab_old(E,f,I,j)*C
      Enddo
      Enddo
      T2resid_ab(A,b,I,j) = T2resid_ab(A,b,I,j) + T
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " Wabef_int2"
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("T2resid_aa:",T2resid_aa,L_aaaa)
      call checksum("T2resid_bb:",T2resid_bb,L_bbbb)
      call checksum("T2resid_ab:",T2resid_ab,L_abab)
#endif

      Return
      End 
