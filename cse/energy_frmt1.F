      Subroutine Energy_frmt1(T1aa_old,T1bb_old,W0_aa,W0_bb,W0_ab,
     +                        Fockov_a,Fockov_b,Nocc_a,Nocc_b,Nvrt_a,
     +                        Nvrt_b,E1_aa,E1_bb,E1_ab) 

      Implicit Double Precision(A-H,O-Z)
      
      Dimension T1aa_old(Nvrt_a,Nocc_a)
      Dimension T1bb_old(Nvrt_b,Nocc_b)

      Dimension W0_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension w0_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension W0_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Fockov_a(Nocc_a,Nvrt_a)
      Dimension Fockov_b(Nocc_b,Nvrt_b)
 
      Integer A,B,I,J

#include "cse_vars.com"

#ifdef _DEBUG_LVLM
      L_aa = Nocc_a*Nvrt_a
      L_bb = Nocc_b*Nvrt_b
      call checksum("T1aa      :",T1aa_old,L_aa)
      call checksum("T1bb      :",T1bb_old,L_bb)
#endif
C AA block 

      E1_aa = 0.0D0
      DO J = 1, Nocc_a
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_a 
      DO A = 1, Nvrt_a
         TC = T1aa_old(A,I)*T1aa_old(B,J)
         TE = T1aa_old(B,I)*T1aa_old(A,J) 
         E1_aa = E1_aa + (TC-TE)*W0_aa(A,B,I,J)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C BB block 

      E1_bb = 0.0D0
      DO J = 1, Nocc_b
      DO I = 1, Nocc_b
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_b
         TC = T1bb_old(a,i)*T1bb_old(b,j)
         TE = T1bb_old(b,i)*T1bb_old(a,j)
         E1_bb = E1_bb + (TC-TE)*W0_bb(a,b,i,j)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C AB block

      E1_ab = 0.0D0
      DO J = 1, Nocc_b
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_a
         E1_ab = E1_ab + T1aa_old(A,I)*T1bb_old(b,j)*W0_ab(A,b,I,j)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

#ifdef _DEBUG_LVLM
      Write(6,"(a,3(1X,F15.10))") "E1aa,E1bb,E1ab       :",E1_aa,E1_bb,
     +                             E1_ab
#endif 
  
      Enhf_aa = 0.0D0
      DO I = 1, Nocc_a
      DO A = 1, Nvrt_a
         Enhf_aa = Enhf_aa + T1aa_old(A,I)*Fockov_a(A,I)
      ENDDo
      ENDDO 
     
      Enhf_bb = 0.0D0
      DO I = 1, Nocc_b
      DO A = 1, Nvrt_b
         Enhf_bb = Enhf_bb + T1aa_old(A,I)*Fockov_a(A,I)
      ENDDo
      ENDDO 

#ifdef _DEBUG_LVL0
      Write(6,"(a,2(1X,F15.10))") "Enhf_aa,Enhf_bb      :",Enhf_aa,
     +                                                     Enhf_bb
#endif
      E1_aa = Enhf_aa
      E1_bb = Enhf_bb
      E1_ab = 0.0D0 

      Return 
      End

