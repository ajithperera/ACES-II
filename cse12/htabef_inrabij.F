      Subroutine Htabef_inrabij(T2aa_old,T2bb_old,T2ab_old,T2resid_aa,
     +                          T2resid_bb,T2resid_ab,Htabef_1,
     +                          Htabef_2,Htabef_3,Nocc_a,Nocc_b,Nvrt_a,
     +                          Nvrt_b)

      Implicit Double Precision(A-H,O-Z)
     
      Dimension T2aa_old(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb_old(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab_old(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension T2resid_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2resid_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2resid_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Htabef_1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension Htabef_2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension Htabef_3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Integer I,J,A,B,E,F

#include "cse12_vars.com"

C AAAA block 

      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
        T = 0.0D0
      Do E = 1, Nvrt_a 
      Do F = 1, Nvrt_a 
        T = T + 0.50D0*Htabef_1(A,B,E,F)*T2aa_old(E,F,I,J)
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
      Do E = 1, Nvrt_b
      Do F = 1, Nvrt_b
         T = T + 0.50D0*Htabef_2(a,b,e,f)*T2bb_old(e,f,i,j)
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
      Do E = 1, Nvrt_a
      Do F = 1, Nvrt_b
         T = T + Htabef_3(A,b,E,f)*T2ab_old(E,f,I,j)
      Enddo
      Enddo
         T2resid_ab(A,b,I,j) = T2resid_ab(A,b,I,j) + T
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " Htabef_inRabij"
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("T2resid_aa:",T2resid_aa,L_aaaa)
      call checksum("T2resid_bb:",T2resid_bb,L_bbbb)
      call checksum("T2resid_ab:",T2resid_ab,L_abab)
      Write(6,*)
#endif
      Return
      End 
