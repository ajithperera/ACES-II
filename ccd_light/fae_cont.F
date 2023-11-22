      Subroutine Fae_cont(T2aa_old,T2bb_old,T2ab_old,Resid_aa,Resid_bb, 
     +                    Resid_ab,Fae_a,Fae_b,Nocc_a,Nocc_b,Nvrt_a,
     +                    Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension T2aa_old(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb_old(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab_old(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)
 
      Dimension Resid_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Resid_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Resid_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Fae_a(Nvrt_a,Nvrt_a)
      Dimension Fae_b(Nvrt_b,Nvrt_b)

      Integer I,J,A,B,E

#include "ccdlight_vars.com"

#ifdef _DEBUG_LVLM
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
      call checksum("Resid_aa  :",Resid_aa,L_aaaa)
      call checksum("Resid_bb  :",Resid_bb,L_bbbb)
      call checksum("Resid_ab  :",Resid_ab,L_abab)
#endif

C AAAA block 

      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T = 0.0D0 
      Do E = 1, Nvrt_a 
         C = (1.0D0-Ocn_va(E))
         T = T + Fae_a(A,E)*T2aa_old(E,B,I,J)*C 
      Enddo 
         Resid_aa(A,B,I,J) = Resid_aa(A,B,I,J) + T 
         Resid_aa(B,A,I,J) = Resid_aa(B,A,I,J) - T 
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
         C = (1.0D0-Ocn_vb(e))
         T = T + Fae_b(a,e)*T2bb_old(e,b,i,j)*C
      Enddo
         Resid_bb(a,b,i,j) = Resid_bb(a,b,i,j) + T 
         Resid_bb(b,a,i,j) = Resid_bb(b,a,i,j) - T 
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
         C = (1.0D0-Ocn_va(E))
         T = T + Fae_a(A,E)*T2ab_old(E,b,I,j)*C 
      Enddo
         Resid_ab(A,b,I,j) = Resid_ab(A,b,I,j) + T 
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVLM
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call output(Fae_b,1,Nvrt_b,1,Nvrt_b,Nvrt_b,Nvrt_b,1) 
#endif 

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do E = 1, Nvrt_b
         C = (1.0D0-Ocn_vb(e))
         T = T + Fae_b(b,e)*T2ab_old(A,e,I,j)*C 
      Enddo
         Resid_ab(A,b,I,j) = Resid_ab(A,b,I,j) + T 
      Enddo
      Enddo
      Enddo
      Enddo


#ifdef _DEBUG_LVL0
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("Resid_aa  :",Resid_aa,L_aaaa)
      call checksum("Resid_bb  :",Resid_bb,L_bbbb)
      call checksum("Resid_ab  :",Resid_ab,L_abab)
#endif

      Return
      End 
