      Subroutine Form_wabef_intms(T1aa,T1bb,W4_aa,W4_bb,W4_ab,W4_ba,
     +                            Wabef_1,Wabef_2,Wabef_3,Nocc_a,
     +                            Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension W4_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W4_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W4_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W4_ba(Nvrt_a,Nvrt_b,Nocc_a,Nvrt_a)

      Dimension Wabef_1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension Wabef_2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension Wabef_3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Integer M,N,A,B,E,F

#include "ccsdlight_vars.com"

C Wabef_1(A,B,E,F) = -P_(A,B)T1(B,M)W(E,F,A,M)

      L_aaaa = Nvrt_a*Nvrt_a*Nvrt_a*Nvrt_a
      L_bbbb = Nvrt_b*Nvrt_b*Nvrt_b*Nvrt_b
      L_abab = Nvrt_a*Nvrt_b*Nvrt_a*Nvrt_b
  
      Call Dzero(Wabef_1,L_aaaa)
      Call Dzero(Wabef_2,L_bbbb)
      Call Dzero(Wabef_3,L_abab)

      Do F = 1, Nvrt_a
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T = 0.0D0
         C = (1.0D0-Ocn_va(A))
      Do M = 1, Nocc_a
         T = T + T1aa(B,M)*W4_aa(E,F,A,M)*C
      Enddo
         Wabef_1(A,B,E,F) =  Wabef_1(A,B,E,F) - T
         Wabef_1(B,A,E,F) =  Wabef_1(B,A,E,F) + T
      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_2(a,b,e,f) = -P_(a,b)T1(b,m)*W(e,f,a,m)

      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         T = 0.0D0
         C = (1.0D0-Ocn_vb(a))
      Do M = 1, Nocc_b
         T = T + T1bb(b,m)*W4_bb(e,f,a,m)*C
      Enddo
         Wabef_2(a,b,e,f) =  Wabef_2(a,b,e,f) - T
         Wabef_2(b,a,e,f) =  Wabef_2(b,a,e,f) + T
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
         C = (1.0D0-Ocn_va(A))
      Do M = 1, Nocc_b
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
         C = (1.0D0-Ocn_vb(b))
      Do M = 1, Nocc_a
         T = T + T1aa(A,M)*W4_ba(E,f,M,b)*C
      Enddo
         Wabef_3(A,b,E,F) = Wabef_3(A,b,E,F) - T
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      call checksum("WABEF_Intm:",Wabef_1,L_aaaa)
      call checksum("Wabef_intm:",Wabef_2,L_bbbb)
      call checksum("WAbEf_intm:",Wabef_3,L_abab)
#endif

      Return
      End
