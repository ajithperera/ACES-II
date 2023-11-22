      Subroutine Htmi_inrabij(T2aa_old,T2bb_old,T2ab_old,T2resid_aa,
     +                        T2resid_bb,T2resid_ab,Htmi_aa,Htmi_bb,
     +                        Nocc_a,Nocc_b,Nvrt_a,Nvrt_b) 

      Implicit Double Precision(A-H,O-Z)

      Dimension T2resid_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2resid_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2resid_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension T2aa_old(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb_old(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab_old(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Htmi_aa(Nocc_a,Nocc_a)
      Dimension Htmi_bb(Nocc_b,Nocc_b)

      Integer A,B,I,J,M

#include "cse12_vars.com"

C Rijab_aa(ab,ij) <- Htmi_aa(M,J)*T2aa(A,B,I,M)

      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a 
         T = T  + Htmi_aa(M,J)*T2aa_old(A,B,I,M)
      Enddo 
         T2resid_aa(A,B,I,J) = T2resid_aa(A,B,I,J) - T 
         T2resid_aa(A,B,J,I) = T2resid_aa(A,B,J,I) + T 
      Enddo 
      Enddo 
      Enddo 
      Enddo 

C Rijab_bb(ab,ij) <- Htmi_bb(M,J)*T2bb(A,B,I,M)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         T = 0.0D0
      Do M = 1, Nocc_b
         T = T + Htmi_bb(m,j)*T2bb_old(a,b,i,m)
      Enddo
         T2resid_bb(a,b,i,j) = T2resid_bb(a,b,i,j) - T
         T2resid_bb(a,b,j,i) = T2resid_bb(a,b,j,i) + T
      Enddo
      Enddo
      Enddo
      Enddo

C Rijab_ab(ab,ij) <- Htmi(a,b)(M,J)*T2ab(A,B,I,M)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_b
         T = T + Htmi_bb(m,j)*T2ab_old(A,b,I,m)
      Enddo
         T2resid_ab(A,b,I,j) = T2resid_ab(A,b,I,j) - T
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
         T = T + Htmi_aa(M,I)*T2ab_old(A,b,M,j)
      Enddo
         T2resid_ab(A,b,I,j) = T2resid_ab(A,b,I,j)  - T
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " Htmi_inrabij"
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("T2resid_aa:",T2resid_aa,L_aaaa)
      call checksum("T2resid_bb:",T2resid_bb,L_bbbb)
      call checksum("T2resid_ab:",T2resid_ab,L_abab)
#endif
      Return
      End 
