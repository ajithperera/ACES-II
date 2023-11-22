      Subroutine Hbar_diagonals(Hbar_mia,Hbar_mib,Hbar_aea,Hbar_aeb,
     +                          Hbar_ijkl1,Hbar_ijkl2,Hbar_ijkl3,
     +                          Hbar_abef1,Hbar_abef2,Hbar_abef3,
     +                          Hbar_mbej1,Hbar_mbej2,Hbar_mbej3,
     +                          Hbar_mbej4,Hbar_mbej5,Hbar_mbej6,
     +                          D1hbar_aa,D1hbar_bb,D2hbar_aa,
     +                          D2hbar_bb,D2hbar_ab,Work,Maxcor,
     +                          Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)
  
      Implicit Double Precision(A-H,O-Z)

#include "ccsdlight_vars.com"
   
      Dimension Work(Maxcor)

      Dimension Hbar_abef1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension Hbar_abef2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension Hbar_abef3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Dimension Hbar_ijkl1(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension Hbar_ijkl2(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension Hbar_ijkl3(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Dimension Hbar_mbej1(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Hbar_mbej2(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Hbar_mbej3(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Hbar_mbej4(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension Hbar_mbej5(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension Hbar_mbej6(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Dimension Hbar_mia(Nocc_a,Nocc_a)
      Dimension Hbar_mib(Nocc_b,Nocc_b)
      Dimension Hbar_aea(Nvrt_a,Nvrt_a)
      Dimension Hbar_aeb(Nvrt_b,Nvrt_b)
      Dimension Hbar_mea(Nocc_a,Nvrt_a)
      Dimension Hbar_meb(Nocc_b,Nvrt_b)

      Dimension D1hbar_aa(Nvrt_a,Nocc_a)
      Dimension D1hbar_bb(Nvrt_b,Nocc_b)
      Dimension D2hbar_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension D2hbar_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension D2hbar_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Data Zero /0.0D0/

      Integer I,J,A,B

c AA diagonal block

      Do I = 1, Nocc_a
      Do A = 1, Nvrt_a
         D1hbar_aa(A,I) = Hbar_aea(A,A) - Hbar_mia(I,I) +
     +                    Hbar_mbej1(I,A,A,I) 
      Enddo
      Enddo

c BB diagonal block

      Do I = 1, Nocc_b
      Do A = 1, Nvrt_b
         D1hbar_bb(a,i) = Hbar_aeb(a,a) - Hbar_mib(i,i) +
     +                    Hbar_mbej2(i,a,a,i)
      Enddo
      Enddo

C AAAA diagonal block (The diagonals of these are zero due to antisymmetry).

      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
         Do B = 1, Nvrt_a
         Do A = 1, Nvrt_a
            If (A .NE. B .AND. I .NE. J) Then
            D2hbar_aa(A,B,I,J) = Hbar_aea(A,A) + Hbar_aea(B,B) -
     +                           Hbar_mia(I,I) - Hbar_mia(J,J) +
     +                           Hbar_mbej1(I,A,A,I)           +
     +                           Hbar_mbej1(J,B,B,J)           +
     +                           Hbar_mbej1(I,B,B,I)           +
     +                           Hbar_mbej1(J,A,A,J)           +
     +                           Hbar_ijkl1(I,J,I,J)           +
     +                           Hbar_abef1(A,B,A,B)
            Else
            D2hbar_aa(A,B,I,J) = Zero
            Endif 
         Enddo                 
         Enddo 
      Enddo
      Enddo


C BBBB diagonal block

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
         Do B = 1, Nvrt_b
         Do A = 1, Nvrt_b
            If (a .NE. b .AND. i .NE. j) Then
            D2hbar_bb(a,b,i,j) = Hbar_aeb(a,a) + Hbar_aeb(b,b) -
     +                           Hbar_mib(i,i) - Hbar_mib(j,j) +
     +                           Hbar_mbej2(i,a,a,i)           +
     +                           Hbar_mbej2(j,b,b,j)           +
     +                           Hbar_mbej2(i,b,b,i)           +
     +                           Hbar_mbej2(j,a,a,j)           +
     +                           Hbar_ijkl2(i,j,i,j)           +
     +                           Hbar_abef2(a,b,a,b) 
            Else
            D2hbar_bb(a,b,i,j) = Zero 
            Endif 
         Enddo
         Enddo
      Enddo
      Enddo

C ABAB diagonal block

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
         Do B = 1, Nvrt_b
         Do A = 1, Nvrt_a
            D2hbar_ab(A,b,I,j) = Hbar_aea(A,A) + Hbar_aeb(b,b) -
     +                           Hbar_mia(I,I) - Hbar_mib(j,j) +
     +                           Hbar_ijkl3(I,j,I,j)           +
     +                           Hbar_abef3(A,b,A,b)           -
     +                           Hbar_mbej5(I,b,I,b)           -
     +                           Hbar_mbej6(j,A,j,A)           + 
     +                           Hbar_mbej1(I,A,A,I)           + 
     +                           Hbar_mbej2(j,b,b,j) 

         Enddo
         Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      L_aa   = Nocc_a*Nvrt_a
      L_bb   = Nocc_b*Nvrt_b
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("D1hbar_aa :",D1hbar_aa,L_aa)
      call checksum("D1hbar_bb :",D1hbar_bb,L_bb)
      call checksum("D2hbar_aa :",D2hbar_aa,L_aaaa)
      call checksum("D2hbar_bb :",D2hbar_bb,L_bbbb)
      call checksum("D2hbar_ab :",D2hbar_ab,L_abab)
#endif

      Return
      End 



