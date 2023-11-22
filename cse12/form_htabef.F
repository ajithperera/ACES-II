      Subroutine Form_htabef(T1aa,T1bb,T1tau_aa,T1tau_bb,T1tau_ab,
     +                       W2_aa,W2_bb,w2_ab,W_aa,W_bb,W_ab,W4_aa,
     +                       W4_bb,W4_ab,W4_ba,Wabef_1,Wabef_2,
     +                       Wabef_3,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension T1tau_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T1tau_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T1tau_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

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

#include "cse12_vars.com"


C Wabef_1(AB,EF) <- + 1/4T1tau(AB,MN)*W(MN,EF)  

      Do F = 1, Nvrt_a
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T = 0.0D0
         C = (1.0D0-Ocn_va(A))*(1.0D0-Ocn_va(B))
         Wabef_1(A,B,E,F) =  W2_aa(A,B,E,F)*C
      Do M = 1, Nocc_a
      Do N = 1, Nocc_a
         T = T + 0.25D0*T1tau_aa(A,B,M,N)*W_aa(M,N,E,F)
      Enddo
      Enddo
         Wabef_1(A,B,E,F) =  Wabef_1(A,B,E,F) + T
      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_1(A,B,E,F) <- -P_(A,B)T1(B,M)W(E,F,A,M) +  W2(AB,EF)

      Do F = 1, Nvrt_a
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
         C = (1.0D0-Ocn_va(A))
         T = T + T1aa(B,M)*W4_aa(E,F,A,M)*C
      Enddo
         Wabef_1(A,B,E,F) =  Wabef_1(A,B,E,F) - T
         Wabef_1(B,A,E,F) =  Wabef_1(B,A,E,F) + T
      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_2(ab,ef) <- + 1/4T1tau(ab,mn)*W(mn,ef) + W2(ab,ef)

      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         T = 0.0D0
         C = (1.0D0-Ocn_vb(a))*(1.0D0-Ocn_vb(b))
         Wabef_2(a,b,e,f) = W2_bb(a,b,e,f)*C
      Do M = 1, Nocc_b
      Do N = 1, Nocc_b
         T = T + 0.25D0*T1tau_bb(a,b,m,n)*W_bb(m,n,e,f)
      Enddo
      Enddo
         Wabef_2(a,b,e,f) =  Wabef_2(a,b,e,f) + T
      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_2(a,b,e,f) <- -P_(a,b)T1(b,m)*W(e,f,a,m) 

      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         T = 0.0D0
      Do M = 1, Nocc_b
         C = (1.0D0-Ocn_vb(a))
         T = T + T1bb(b,m)*W4_bb(e,f,a,m)*C
      Enddo
         Wabef_2(a,b,e,f) =  Wabef_2(a,b,e,f) - T
         Wabef_2(b,a,e,f) =  Wabef_2(b,a,e,f) + T
      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_3(Ab,Ef) = -T1(b,m)*W(E,f,A,m) + W2(Ab,eF)

      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
         C = (1.0D0-Ocn_va(A))*(1.0D0-Ocn_vb(b))
         Wabef_3(A,b,E,F) = W2_ab(A,b,E,f)*C
      Do M = 1, Nocc_b
         C = (1.0D0-Ocn_va(A))
         T = T + T1bb(b,m)*W4_ab(E,f,A,m)*C
      Enddo
         Wabef_3(A,b,E,f) = Wabef_3(A,b,E,f) - T
      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_3(Ab,Ef) <- T1(A,M)*W(E,f,M,b)

      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
         C = (1.0D0-Ocn_vb(b))
         T = T + T1aa(A,M)*W4_ba(E,f,M,b)*C
      Enddo
         Wabef_3(A,b,E,f) = Wabef_3(A,b,E,f) - T
      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_3(Ab,Ef) <- T1tau_ab(Ab,Mn)*W(Mn,Ef)

      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
      Do N = 1, Nocc_b
         T = T + 0.50D0*T1tau_ab(A,b,M,n)*W_ab(M,n,E,f)
      Enddo
      Enddo
         Wabef_3(A,b,E,f)= Wabef_3(A,b,E,f) + T
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " @-form_htabef"
      L_aaaa = Nvrt_a*Nvrt_a*Nvrt_a*Nvrt_a
      L_bbbb = Nvrt_b*Nvrt_b*Nvrt_b*Nvrt_b
      L_abab = Nvrt_a*Nvrt_b*Nvrt_a*Nvrt_b
      call checksum("Form_wABEF:",Wabef_1,L_aaaa)
      call checksum("Form_wabef:",Wabef_2,L_bbbb)
      call checksum("Form_wAbEf:",Wabef_3,L_abab)
#endif

      Return
      End
