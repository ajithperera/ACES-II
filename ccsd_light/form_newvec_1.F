      Subroutine Form_newvec_1(D1hbar_aa,D1hbar_bb,D2hbar_aa,
     +                         D2hbar_bb,D2hbar_ab,T1resid_aa,
     +                         T1resid_bb,T2resid_aa,T2resid_bb,
     +                         T2resid_ab,Ht1_aiaa,
     +                         Ht1_aibb,Ht2_abij1,Ht2_abij2,
     +                         Ht2_abij3,D1invt_aa,D1invt_bb,
     +                         D2invt_aa,D2invt_bb,D2invt_ab,
     +                         Newt10_aa,Newt10_bb,Newt20_aa,
     +                         Newt20_bb,Newt20_ab,Newt11_aa,
     +                         Newt11_bb,Newt21_aa,Newt21_bb,
     +                         Newt21_ab,Nocc_a,Nocc_b,Nvrt_a,
     +                         Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Double Precision Newt10_aa,Newt10_bb,Newt20_aa,Newt20_bb,
     +                 Newt20_ab
      Double Precision Newt11_aa,Newt11_bb,Newt21_aa,Newt21_bb,
     +                 Newt21_ab

#include "ccsdlight_vars.com"

      Dimension T1resid_aa(Nvrt_a,Nocc_a)
      Dimension T1resid_bb(Nvrt_b,Nocc_b)
      Dimension T2resid_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2resid_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2resid_ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension D1hbar_aa(Nvrt_a,Nocc_a)
      Dimension D1hbar_bb(Nvrt_b,Nocc_b)
      Dimension D2hbar_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension D2hbar_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension D2hbar_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension D1invt_aa(Nvrt_a,Nocc_a)
      Dimension D1invt_bb(Nvrt_b,Nocc_b)
      Dimension D2invt_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension D2invt_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension D2invt_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Ht1_aiaa(Nvrt_a,Nocc_a)
      Dimension Ht1_aibb(Nvrt_b,Nocc_b)
      Dimension Ht2_abij1(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Ht2_abij2(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Ht2_abij3(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Newt10_aa(Nvrt_a,Nocc_a)
      Dimension Newt10_bb(Nvrt_b,Nocc_b)
      Dimension Newt20_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Newt20_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Newt20_ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension Newt11_aa(Nvrt_a,Nocc_a)
      Dimension Newt11_bb(Nvrt_b,Nocc_b)
      Dimension Newt21_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Newt21_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Newt21_ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      L_aa   = Nocc_a*Nvrt_a
      L_bb   = Nocc_b*Nvrt_b
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b

#ifdef _DEBUG_LVLM
      Write(6,"(a)") "@-Form_newvec_1"
      call checksum("D1invt_aa :",D1invt_aa,L_aa)
      call checksum("D1invt_bb :",D1invt_bb,L_bb)
      call checksum("D2invt_aa :",D2invt_aa,L_aaaa)
      call checksum("D2invt_bb :",D2invt_bb,L_bbbb)
      call checksum("D2invt_ab :",D2invt_ab,L_abab)
      call checksum("D1hbar_aa :",D1invt_aa,L_aa)
      call checksum("D1hbar_bb :",D1invt_bb,L_bb)
      call checksum("D2hbar_aa :",D2hbar_aa,L_aaaa)
      call checksum("D2hbar_bb :",D2hbar_bb,L_bbbb)
      call checksum("D2hbar_ab :",D2hbar_ab,L_abab)
      call checksum("Ht1_aiaa  :",Ht1_aiaa,L_aa)
      call checksum("Ht1_aibb  :",Ht1_aibb,L_bb)
      call checksum("Ht2_abij1 :",Ht2_abij1,L_aaaa)
      call checksum("Ht2_abij2 :",Ht2_abij2,L_bbbb)
      call checksum("Ht2_abij3 :",Ht2_abij3,L_abab)
#endif
C Obtain the norm of the previous solution vector.

      D_norm_aa   = Ddot(L_aa,Newt10_aa,1,Newt10_aa,1)
      D_norm_bb   = Ddot(L_bb,Newt10_bb,1,Newt10_bb,1)
      D_norm_aaaa = Ddot(L_aaaa,Newt20_aa,1,Newt20_aa,1)
      D_norm_bbbb = Ddot(L_bbbb,Newt20_bb,1,Newt20_bb,1)
      D_norm_abab = Ddot(L_abab,Newt20_ab,1,Newt20_ab,1)

      D_norm_old = D_norm_aa+D_norm_bb+D_norm_aaaa+D_norm_bbbb
     +             D_norm_abab
      Write(6,"(a,F15.10)") " The norm of DeltaT(k) :", D_norm_old
C AA block 

      DO I = 1, Nocc_a
      DO A = 1, Nvrt_a
         Dinvt  = D1invt_aa(A,I)
         Dinvht = D1hbar_aa(A,I)*Ht1_aiaa(A,I)
         T1     = Newt10_aa(A,I)
         R1     = D1hbar_aa(A,I)*T1resid_aa(A,I)
         Newt11_aa(A,I) = Dinvt - (Dinvht - T1) + D_norm_old*R1
      ENDDO
      ENDDO

C BB block

      DO I = 1, Nocc_b
      DO A = 1, Nvrt_b
         Dinvt  = D1invt_bb(a,i)
         Dinvht = D1hbar_bb(a,i)*Ht1_aibb(a,i)
         T1     = Newt10_bb(a,i)
         R1     = D1hbar_bb(a,i)*T1resid_bb(a,i)
         Newt11_bb(a,i) = Dinvt - (Dinvht - T1) + D_norm_old*R1
      ENDDO
      ENDDO

C AAAA block

      DO J = 1, Nocc_a
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_a
      DO A = 1, Nvrt_a
         Dinvt  = D2invt_aa(A,B,I,J)
         Dinvht = D2hbar_aa(A,B,I,J)*Ht2_abij1(A,B,I,J)
         T2     = Newt20_aa(A,B,I,J)
         R2     = D2hbar_aa(A,B,I,J)*T2resid_aa(A,B,I,J)
         Newt21_aa(A,B,I,J) =  Dinvt - (Dinvht - T2) + D_norm_old*R2
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C BBBB block

      DO J = 1, Nocc_b
      DO I = 1, Nocc_b
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_b
         Dinvt  = D2invt_bb(a,b,i,j)
         Dinvht = D2hbar_bb(a,b,i,j)*Ht2_abij2(a,b,i,j)
         T2     = Newt20_bb(a,b,i,j)
         R2     = D2hbar_bb(a,b,i,j)*T2resid_bb(a,b,i,j)
         Newt21_bb(A,B,I,J) =  Dinvt - (Dinvht - T2) + D_norm_old*R2
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C ABAB block

      DO J = 1, Nocc_b
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_a
         Dinvt  = D2invt_ab(A,b,I,j)
         Dinvht = D2hbar_ab(A,b,I,j)*Ht2_abij3(A,b,I,j)
         T2     = Newt20_ab(A,b,I,j)
         R2     = D2hbar_ab(A,b,I,j)*T2resid_ab(A,b,I,j)
         Newt21_ab(A,B,I,J) =  Dinvt - (Dinvht - T2) + D_norm_old*R2
      ENDDO
      ENDDO
      ENDDO
      ENDDO

#ifdef _DEBUG_LVLM
      L_aa = Nocc_a*Nvrt_a
      L_bb = Nocc_b*Nvrt_b
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("Deltat1_aa:",Newt11_aa,L_aa)
      call checksum("Deltat1_bb:",Newt11_bb,L_bb)
      call checksum("Deltat2_aa:",Newt21_aa,L_aaaa)
      call checksum("Deltat2_bb:",Newt21_bb,L_bbbb)
      call checksum("Deltat2_ab:",Newt21_ab,L_abab)
#endif

      Return
      End
