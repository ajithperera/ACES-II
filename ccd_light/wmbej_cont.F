      Subroutine Wmbej_cont(T2aa_old,T2bb_old,T2ab_old,Resid_aa,
     +                      Resid_bb,Resid_ab,Wmbej_1,Wmbej_2,
     +                      Wmbej_3,Wmbej_4,Wmbej_5,Wmbej_6,
     +                      Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)


      Dimension T2aa_old(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb_old(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab_old(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)
 
      Dimension Resid_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Resid_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Resid_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Wmbej_1(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Wmbej_2(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Wmbej_3(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Wmbej_4(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension Wmbej_5(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension Wmbej_6(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Integer A,B,I,J,M,E

#include "ccdlight_vars.com"

#ifdef _DEBUG_LVLM
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("T2aa      :",T2aa_old,L_aaaa)
      call checksum("T2bb      :",T2bb_old,L_bbbb)
      call checksum("T2ab      :",T2ab_old,L_abab)
      L_aaaa = Nocc_a*Nvrt_a*Nvrt_a*Nocc_a
      L_bbbb = Nocc_b*Nvrt_b*Nvrt_b*Nocc_b
      L_abab = Nocc_a*Nvrt_b*Nvrt_a*Nocc_b
      L_baba = Nocc_b*Nvrt_a*Nvrt_b*Nocc_a
      L_abba = Nocc_a*Nvrt_b*Nocc_a*Nvrt_b
      L_baab = Nocc_b*Nvrt_a*Nocc_b*Nvrt_a
      call checksum("Form_wMBEJ:",Wmbej_1,L_aaaa)
      call checksum("Form_wmbej:",Wmbej_2,L_bbbb)
      call checksum("Form_wMbEj:",Wmbej_3,L_abab)
      call checksum("Form_wmBeJ:",Wmbej_4,L_baba)
      call checksum("Form_wmBEj:",Wmbej_5,L_abba)
      call checksum("Form_wMbeJ:",Wmbej_6,L_baab)
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
      Do M = 1, Nocc_a 
      Do E = 1, Nvrt_a 
         C = (1.0D0-Ocn_va(E))*Ocn_oa(M)
         T = T + Wmbej_1(M,B,E,J)*T2aa_old(A,E,I,M)*C 
      Enddo 
      Enddo 

      Resid_aa(A,B,I,J) = Resid_aa(A,B,I,J) + T
      Resid_aa(B,A,I,J) = Resid_aa(B,A,I,J) - T
      Resid_aa(A,B,J,I) = Resid_aa(A,B,J,I) - T
      Resid_aa(B,A,J,I) = Resid_aa(B,A,J,I) + T

      Enddo 
      Enddo 
      Enddo 
      Enddo 

      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_b
      Do E = 1, Nvrt_b
         C = (1.0D0-Ocn_vb(e))*Ocn_ob(m)
         T  = T + Wmbej_4(m,B,e,J)*T2ab_old(A,e,I,m)*C
      Enddo
      Enddo

      Resid_aa(A,B,I,J) = Resid_aa(A,B,I,J) + T 
      Resid_aa(B,A,I,J) = Resid_aa(B,A,I,J) - T 
      Resid_aa(A,B,J,I) = Resid_aa(A,B,J,I) - T 
      Resid_aa(B,A,J,I) = Resid_aa(B,A,J,I) + T 

      Enddo
      Enddo
      Enddo
      Enddo

C BBBB block 

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         T  = 0.0D0
      Do M = 1, Nocc_b
      Do E = 1, Nvrt_b
         C = (1.0D0-Ocn_vb(e))*Ocn_ob(m)
         T = T + Wmbej_2(m,b,e,j)*T2bb_old(a,e,i,m)*C
      Enddo
      Enddo

      Resid_bb(a,b,i,j) = Resid_bb(a,b,i,j) + T
      Resid_bb(b,a,i,j) = Resid_bb(b,a,i,j) - T
      Resid_bb(a,b,j,i) = Resid_bb(a,b,j,i) - T
      Resid_bb(b,a,j,i) = Resid_bb(b,a,j,i) + T

      Enddo
      Enddo
      Enddo
      Enddo 

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         T  = 0.0D0
      Do M = 1, Nocc_a
      Do E = 1, Nvrt_a
         C = (1.0D0-Ocn_va(E))*Ocn_oa(M)
         T = T + Wmbej_3(M,b,E,j)*T2ab_old(E,a,M,i)*C 
      Enddo
      Enddo

      Resid_bb(a,b,i,j) = Resid_bb(a,b,i,j) + T 
      Resid_bb(b,a,i,j) = Resid_bb(b,a,i,j) - T 
      Resid_bb(a,b,j,i) = Resid_bb(a,b,j,i) - T 
      Resid_bb(b,a,j,i) = Resid_bb(b,a,j,i) + T 

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
      Do M = 1, Nocc_a
      Do E = 1, Nvrt_a
         C = (1.0D0-Ocn_va(E))*Ocn_oa(M)
         T = T +  Wmbej_1(M,A,E,I)*T2ab_old(E,b,M,j)*C
      Enddo
      Enddo
      Resid_ab(A,b,I,j) = Resid_ab(A,b,I,j) + T 

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
        C = (1.0D0-Ocn_vb(e))*Ocn_ob(m)
        T = T +  Wmbej_2(m,b,e,j)*T2ab_old(A,e,I,m)*C 
      Enddo
      Enddo

      Resid_ab(A,b,I,j) = Resid_ab(A,b,I,j) + T 

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
      Do E = 1, Nvrt_a
         C = (1.0D0-Ocn_va(E))*Ocn_oa(M)
         T = T +  Wmbej_3(M,b,E,j)*T2aa_old(A,E,I,M)*C
      Enddo
      Enddo

      Resid_ab(A,b,I,j) = Resid_ab(A,b,I,j) + T 

      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_a
      Do I = 1, Nocc_b
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_b
         T = 0.0D0
      Do M = 1, Nocc_b
      Do E = 1, Nvrt_b
        C = (1.0D0-Ocn_vb(e))*Ocn_ob(m)
        T = T + Wmbej_4(m,B,e,J)*T2bb_old(a,e,i,m)*C 
      Enddo
      Enddo

      Resid_ab(B,a,J,i) = Resid_ab(B,a,J,i) + T 

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
      Do E = 1, Nvrt_b
         C = (1.0D0-Ocn_vb(e))*Ocn_oa(M)
         T  = T -  Wmbej_5(M,b,I,e)*T2ab_old(A,e,M,j)*C 
      Enddo
      Enddo

      Resid_ab(A,b,I,j) = Resid_ab(A,b,I,j) + T 

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
      Do E = 1, Nvrt_a
         C = (1.0D0-Ocn_va(E))*Ocn_ob(m)
         T = T - Wmbej_6(m,A,j,E)*T2ab_old(E,b,I,m)*C 
      Enddo
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
