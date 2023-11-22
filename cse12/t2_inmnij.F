      Subroutine T2_inmnij(T2aa,T2bb,T2ab,W_aa,W_bb,W_ab,Wmnij_1,
     +                     Wmnij_2,Wmnij_3,Nocc_a,Nocc_b,Nvrt_a,
     +                     Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_b)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension Wmnij_1(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension Wmnij_2(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension Wmnij_3(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Integer M,N,I,J,E,F

#include "cse12_vars.com" 

C Wmnij_1(MN,IJ) <-  1/4T2(IJ,EF)*W(MN,EF) 
 
      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do N = 1, Nocc_a
      Do M = 1, Nocc_a
         T = 0.0D0
      Do E = 1, Nvrt_a
      Do F = 1, Nvrt_a
         T = T + 0.25D0*T2aa(E,F,I,J)*W_aa(M,N,E,F)
      Enddo
      Enddo
          Wmnij_1(M,N,I,J) = T 
      Enddo
      Enddo
      Enddo
      Enddo

C Wmnij_2(mn,ij) <-  1/4T2(ij,ef)*W2(mn,ef) 

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
      Do N = 1, Nocc_b
      Do M = 1, Nocc_b
         T = 0.0D0
      Do E = 1, Nvrt_b
      Do F = 1, Nvrt_b
         T = T + 0.25D0*T2bb(e,f,i,j)*W_bb(m,n,e,f)
      Enddo
      Enddo
         Wmnij_2(m,n,i,j) = T
      Enddo
      Enddo
      Enddo
      Enddo

C Wmnij_2(Mn,Ij) <- T2(Ij,Ef)*W2(Mn,Ef)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do N = 1, Nocc_b
      Do M = 1, Nocc_a
         T = 0.0D0
      Do E = 1, Nvrt_a
      Do F = 1, Nvrt_b
         T  = T + 0.50D0*T2ab(E,f,I,j)*W_ab(M,n,E,f)
      Enddo
      Enddo
        Wmnij_3(M,n,I,j) = T
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " @-T2_inmnij" 
      L_aaaa = Nocc_a*Nocc_a*Nocc_a*Nocc_a
      L_bbbb = Nocc_b*Nocc_b*Nocc_b*Nocc_b
      L_abab = Nocc_a*Nocc_b*Nocc_a*Nocc_b
      call checksum("T2->MNIJ  :",Wmnij_1,L_aaaa)
      call checksum("T2->mnij  :",Wmnij_2,L_bbbb)
      call checksum("T2->MnIj  :",Wmnij_3,L_abab)
#endif

      Return
      End
