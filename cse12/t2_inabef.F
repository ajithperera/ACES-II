      Subroutine T2_inabef(T2aa,T2bb,T2ab,W_aa,W_bb,W_ab,Wabef_1,
     +                     Wabef_2,Wabef_3,Nocc_a,Nocc_b,Nvrt_a,
     +                     Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension Wabef_1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension Wabef_2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension Wabef_3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Integer M,N,A,B,E,F

#include "cse12_vars.com"
      
C Wabef_1(AB,EF) <- + 1/4T2(AB,MN)*W(MN,EF) 
 
      Do F = 1, Nvrt_a
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
      Do N = 1, Nocc_a
         T = T + 0.25D0*T2aa(A,B,M,N)*W_aa(M,N,E,F)
      Enddo
      Enddo
         Wabef_1(A,B,E,F) = T 
      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_2(ab,ef) = + 1/4T2(ab,mn)*W(mn,ef) 

      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         T = 0.0D0
      Do M = 1, Nocc_b
      Do N = 1, Nocc_b
         T = T + 0.25D0*T2bb(a,b,m,n)*W_bb(m,n,e,f)
      Enddo
      Enddo
          Wabef_2(a,b,e,f) = T
      Enddo
      Enddo
      Enddo
      Enddo

C Wabef_3(Ab,Ef) <- W2(Ab,Ef) + T2(Ab,Mn)*W(Mn,Ef)

      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
      Do N = 1, Nocc_b
         T = T + 0.50D0*T2ab(A,b,M,n)*W_ab(M,n,E,f)
      Enddo
      Enddo 
         Wabef_3(a,b,e,f) = T
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " @-T2_inabef" 
      L_aaaa = Nvrt_a*Nvrt_a*Nvrt_a*Nvrt_a
      L_bbbb = Nvrt_b*Nvrt_b*Nvrt_b*Nvrt_b
      L_abab = Nvrt_a*Nvrt_b*Nvrt_a*Nvrt_b
      call checksum("T2->ABEF  :",Wabef_1,L_aaaa)
      call checksum("T2->abef  :",Wabef_2,L_bbbb)
      call checksum("T2->AbEf  :",Wabef_3,L_abab)
#endif
      Return
      End
