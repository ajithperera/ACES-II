      Subroutine Form_abef(T2aa,T2bb,T2ab,T1aa,T1bb,W2_aa,W2_bb,w2_ab,
     +                     W_aa,W_bb,W_ab,W4_aa,W4_bb,W4_ab,W4_ba,
     +                     Wabef_1,Wabef_2,Wabef_3,Nocc_a,Nocc_b,
     +                     Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension W2_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension W2_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension W2_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Dimension W4_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W4_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W4_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W4_ba(Nvrt_a,Nvrt_b,Nocc_a,Nvrt_a)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension Wabef_1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension Wabef_2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension Wabef_3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Integer M,N,A,B,E,F

#include "ccsdlight_vars.com"

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
         C = (1.0D0-Ocn_va(A))*(1.0D0-Ocn_va(B))
         Wabef_1(A,B,E,F) =  W2_aa(A,B,E,F)*C
      Do M = 1, Nocc_a
      Do N = 1, Nocc_a
CSSS      C = Ocn_oa(M)*Ocn_oa(N)
      C = 1.0D0
      Wabef_1(A,B,E,F) = Wabef_1(A,B,E,F) + 
     +                   0.25D0*T2aa(A,B,M,N)*W_aa(M,N,E,F)*C 
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_1(A,B,E,F) = -P_(A,B)T1(B,M)W(E,F,A,M)

      Do F = 1, Nvrt_a
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T  = 0.0D0
         T1 = 0.0D0
         T2 = 0.0D0
      Do M = 1, Nocc_a
CSSS         C =  Ocn_oa(M)
         C  = (1.0D0-Ocn_va(A))
         C1 = (1.0D0-Ocn_va(A))
         C2 = (1.0D0-Ocn_va(B))
         T  = T  + T1aa(B,M)*W4_aa(E,F,A,M)*C
         T1 = T1 + T1aa(B,M)*W4_aa(E,F,A,M)*C1
         T2 = T2 + T1aa(A,M)*W4_aa(E,B,A,M)*C2
      Enddo
         Wabef_1(A,B,E,F) =  Wabef_1(A,B,E,F) - T1
         Wabef_1(B,A,E,F) =  Wabef_1(B,A,E,F) + T2
      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_2(ab,ef) = W2(ab,ef) + 1/4T2(ab,mn)*W(mn,ef) 

      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         C = (1.0D0-Ocn_vb(a))*(1.0D0-Ocn_vb(b))
         Wabef_2(a,b,e,f) = W2_bb(a,b,e,f)*C
      Do M = 1, Nocc_b
      Do N = 1, Nocc_b
CSSS      C = Ocn_ob(m)*Ocn_ob(n)
      C = 1.0D0
      Wabef_2(a,b,e,f) = Wabef_2(a,b,e,f) + 
     +                   0.25D0*T2bb(a,b,m,n)*W_bb(m,n,e,f)*C 
      Enddo
      Enddo

      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_2(a,b,e,f) = -P_(a,b)T1(b,m)*W(e,f,a,m)

      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         T  = 0.0D0
         T1 = 0.0D0
         T2 = 0.0D0
      Do M = 1, Nocc_b
CSSS         C =  Ocn_ob(m)
         C  = (1.0D0-Ocn_vb(a))
         C1 = (1.0D0-Ocn_vb(a))
         C2 = (1.0D0-Ocn_vb(b))
         T  = T  + T1bb(b,m)*W4_bb(e,f,a,m)*C
         T1 = T1 + T1bb(b,m)*W4_bb(e,f,a,m)*C1
         T2 = T2 + T1bb(a,m)*W4_bb(e,f,b,m)*C2
      Enddo
         Wabef_2(a,b,e,f) =  Wabef_2(a,b,e,f) - T1
         Wabef_2(b,a,e,f) =  Wabef_2(b,a,e,f) + T2
      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_3(Ab,Ef) = W2(Ab,Ef) + T2(Ab,Mn)*W(Mn,Ef)

      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         C = (1.0D0-Ocn_va(A))*(1.0D0-Ocn_vb(b))
         Wabef_3(A,b,E,F) = W2_ab(A,b,E,f)*C

CSSS      Do M = 1, Nocc_a
CSSS      Do N = 1, Nocc_b
CSSS      Wabef_3(A,b,E,f)= Wabef_3(A,b,E,f)
CSSS      Enddo
CSSS      Enddo

      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_3(Ab,Ef) = -T1(b,m)*W(E,f,A,m) 
      
      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a  
         T = 0.0D0
      Do M = 1, Nocc_b
CSSS         C = Ocn_ob(m)
         C = (1.0D0-Ocn_va(A))
         T = T + T1bb(b,m)*W4_ab(E,f,A,m)*C
      Enddo
         Wabef_3(A,b,E,F) = Wabef_3(A,b,E,F) - T
      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_3(Ab,Ef) = -T1(A,M)*W(E,f,b,M)->T1(A,M)*W(E,f,M,b)

      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
CSSS         C = Ocn_oa(M)
         C = (1.0D0-Ocn_vb(b))
         T = T + T1aa(A,M)*W4_ba(E,f,M,b)*C
      Enddo
         Wabef_3(A,b,E,F) = Wabef_3(A,b,E,F) - T
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
