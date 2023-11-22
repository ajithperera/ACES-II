      Subroutine T2mnij_inrabij(T2aa_old,T2bb_old,T2ab_old,T2resid_aa,
     +                          T2resid_bb,T2resid_ab,Wmnij_1,
     +                          Wmnij_2,Wmnij_3,Nocc_a,Nocc_b,Nvrt_a,
     +                          Nvrt_b)

      Implicit Double Precision(A-H,O-Z)
    
      Dimension T2resid_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2resid_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2resid_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension T2aa_old(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb_old(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab_old(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Wmnij_1(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension Wmnij_2(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension Wmnij_3(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Integer I,J,A,B,M,N

#include "cse12_vars.com"

C AAAA block 

      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a 
      Do N = 1, Nocc_a 
         T = T + 0.50D0*Wmnij_1(M,N,I,J)*T2aa_old(A,B,M,N)
      Enddo 
      Enddo 
         T2resid_aa(A,B,I,J) = T2resid_aa(A,B,I,J) + T 
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
      Do M = 1, Nocc_b
      Do N = 1, Nocc_b
         T = T + 0.50D0*Wmnij_2(m,n,i,j)*T2bb_old(a,b,m,n)
      Enddo
      Enddo
         T2resid_bb(a,b,i,j) = T2resid_bb(a,b,i,j) + T 
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
      Do N = 1, Nocc_b
         T = T + Wmnij_3(M,n,I,j)*T2ab_old(A,b,M,n)
      Enddo
      Enddo
         T2resid_ab(A,b,I,j) = T2resid_ab(A,b,I,j) + T 
      Enddo
      Enddo
      Enddo
      Enddo 

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " @-T2mnij_inrabij"
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("T2resid_aa:",T2resid_aa,L_aaaa)
      call checksum("T2resid_bb:",T2resid_bb,L_bbbb)
      call checksum("T2resid_ab:",T2resid_ab,L_abab)
#endif

      Return
      End 
