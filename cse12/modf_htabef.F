      Subroutine Modf_htabef(T1aa,T1bb,W2_aa,W2_bb,w2_ab,W4_aa,W4_bb,
     +                       W4_ab,W4_ba,Wabef_1,Wabef_2,Wabef_3,
     +                       Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension W2_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension W2_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension W2_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Dimension W4_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W4_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W4_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W4_ba(Nvrt_a,Nvrt_b,Nocc_a,Nvrt_a)

      Dimension Wabef_1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension Wabef_2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension Wabef_3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Integer M,N,A,B,E,F

#include "cse12_vars.com"

C These simply subtract these contributions since these terms have 
C have now twice what they suppose to be 

C Wabef_1(AB,EF) <- W2_aa(A,B,E,F) 

      Do F = 1, Nvrt_a
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         C = (1.0D0-Ocn_va(A))*(1.0D0-Ocn_va(B))
         Wabef_1(A,B,E,F) = Wabef_1(A,B,E,F) - W2_aa(A,B,E,F)*C
      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_1(A,B,E,F) <- -P_(A,B)T1(B,M)W(E,F,A,M) 

      Do F = 1, Nvrt_a
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
         C = (1.0D0-Ocn_va(A))
         T = T + T1aa(B,M)*W4_aa(E,F,A,M)*C
      Enddo
         Wabef_1(A,B,E,F) =  Wabef_1(A,B,E,F) + T
         Wabef_1(B,A,E,F) =  Wabef_1(B,A,E,F) - T
      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_2(a,b,e,f) <- W2(ab,ef)

      Do F = 1, Nvrt_a
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         C = (1.0D0-Ocn_vb(a))*(1.0D0-Ocn_vb(b))
         Wabef_2(a,b,e,f) = Wabef_2(a,b,e,f) - W2_bb(a,b,e,f)*C 
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
         Wabef_2(a,b,e,f) =  Wabef_2(a,b,e,f) + T
         Wabef_2(b,a,e,f) =  Wabef_2(b,a,e,f) - T
      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_3(a,b,e,f) <- W2_ab(A,b,E,f)

      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
         C = (1.0D0-Ocn_va(A))*(1.0D0-Ocn_vb(b))
         Wabef_3(A,b,E,F) = Wabef_3(A,b,E,F) - W2_ab(A,b,E,f)*C
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
         C = (1.0D0-Ocn_va(A))
         T = T + T1bb(b,m)*W4_ab(E,f,A,m)*C
      Enddo
         Wabef_3(A,b,E,f) = Wabef_3(A,b,E,f) + T
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
         Wabef_3(A,b,E,f) = Wabef_3(A,b,E,f) + T
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " @-modf_htabef"
      L_aaaa = Nvrt_a*Nvrt_a*Nvrt_a*Nvrt_a
      L_bbbb = Nvrt_b*Nvrt_b*Nvrt_b*Nvrt_b
      L_abab = Nvrt_a*Nvrt_b*Nvrt_a*Nvrt_b
      call checksum("New_wABEF :",Wabef_1,L_aaaa)
      call checksum("New_wabef :",Wabef_2,L_bbbb)
      call checksum("New_wAbEf :",Wabef_3,L_abab)
#endif

      Return
      End
