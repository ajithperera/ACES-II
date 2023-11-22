      Subroutine Htabef_inhtabij(T2aa_old,T2bb_old,T2ab_old,Htabij_1,
     +                           Htabij_2,Htabij_3,Htabef_1,Htabef_2,
     +                           Htabef_3,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)
     
      Dimension T2aa_old(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb_old(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab_old(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)
 
      Dimension Htabij_1(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Htabij_2(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Htabij_3(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

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
         Htabij_1(A,B,I,J) = Htabij_1(A,B,I,J) + T 
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
         Htabij_2(a,b,i,j) = Htabij_2(a,b,i,j) + T 
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
         Htabij_3(A,b,I,j) = Htabij_3(A,b,I,j) + T
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " @-Htabef_inhtabij"
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("Htabij_aa :",Htabij_1,L_aaaa)
      call checksum("Htabij_bb :",Htabij_2,L_bbbb)
      call checksum("Htabij_ab :",Htabij_3,L_abab)
      Write(6,*)
#endif
      Return
      End 
