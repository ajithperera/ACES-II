      Subroutine Form_abef(T2aa,T2bb,T2ab,W2_aa,W2_bb,W2_ab,W_aa,
     +                     W_bb,W_ab,Wabef_1,Wabef_2,Wabef_3,
     +                     Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension W2_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension W2_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension W2_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension Wabef_1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension Wabef_2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension Wabef_3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Integer M,N,A,B,E,F

#include "ccdlight_vars.com"

#ifdef _DEBUG_LVLM
      Write(6,*)
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      call checksum("T2aa      :",T2aa,L_aaaa)
      call checksum("T2bb      :",T2bb,L_bbbb)
      call checksum("T2ab      :",T2ab,L_abab)
      L_aaaa = Nvrt_a*Nvrt_a*Nvrt_a*Nvrt_a
      L_bbbb = Nvrt_b*Nvrt_b*Nvrt_b*Nvrt_b
      L_abab = Nvrt_a*Nvrt_b*Nvrt_a*Nvrt_b
      call checksum("Form_wABEF:",W2_aa,L_aaaa)
      call checksum("Form_wabef:",W2_bb,L_bbbb)
      call checksum("Form_wAbEf:",W2_ab,L_abab)
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("Form_wMNEF:",W_aa,L_aaaa)
      call checksum("Form_wmnef:",W_bb,L_bbbb)
      call checksum("Form_wMnEf:",W_ab,L_abab)
#endif
      
C Wabef_1(AB,EF) = W2(AB,EF) + 1/4T2(AB,MN)*W(MN,EF) 
 
      Do F = 1, Nvrt_a
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
      Wabef_1(A,B,E,F) =  W2_aa(A,B,E,F)
      Do M = 1, Nocc_a
      Do N = 1, Nocc_a
      C = Ocn_oa(M)*Ocn_oa(N)
      Wabef_1(A,B,E,F) = Wabef_1(A,B,E,F) + 
     +                   0.25D0*T2aa(A,B,M,N)*W_aa(M,N,E,F)*C 
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_2(ab,ef) = W2(ab,ef) + 1/4T2(ab,mn)*W(mn,ef) 

      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         Wabef_2(a,b,e,f) = W2_bb(a,b,e,f) 
      Do M = 1, Nocc_b
      Do N = 1, Nocc_b
      C = Ocn_ob(m)*Ocn_ob(n)
      Wabef_2(a,b,e,f) = Wabef_2(a,b,e,f) + 
     +                   0.25D0*T2bb(a,b,m,n)*W_bb(m,n,e,f)*C 
      Enddo
      Enddo

      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_3(Ab,Ef) = W2(Ab,Ef) + T2(Ab,Mn)*W(Mn,Ef)

      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         Wabef_3(A,b,E,F) = W2_ab(A,b,E,f)

CSSS      Do M = 1, Nocc_a
CSSS      Do N = 1, Nocc_b
CSSS      Wabef_3(A,b,E,f)= Wabef_3(A,b,E,f)
CSSS      Enddo
CSSS      Enddo

      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      L_aaaa = Nvrt_a*Nvrt_a*Nvrt_a*Nvrt_a
      L_bbbb = Nvrt_b*Nvrt_b*Nvrt_b*Nvrt_b
      L_abab = Nvrt_a*Nvrt_b*Nvrt_a*Nvrt_b
      call checksum("Form_wABEF:",Wabef_1,L_aaaa)
      call checksum("Form_wabef:",Wabef_2,L_bbbb)
      call checksum("Form_wAbEf:",Wabef_3,L_abab)
#endif
      Return
      End
