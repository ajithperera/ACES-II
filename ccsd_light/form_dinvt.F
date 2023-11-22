      Subroutine Form_dinvt(D1hbar_aa,D1hbar_bb,D2hbar_aa,D2hbar_bb,
     +                      D2hbar_ab,D1invt_aa,D1invt_bb,D2invt_aa,
     +                      D2invt_bb,D2invt_ab,T1aa,T1bb,T2aa,T2bb,
     +                      T2ab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

#include "ccsdlight_vars.com"

      Dimension D1hbar_aa(Nvrt_a,Nocc_a)
      Dimension D1hbar_bb(Nvrt_b,Nocc_b)
      Dimension D2hbar_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension D2hbar_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension D2hbar_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension D1invt_aa(Nvrt_a,Nocc_a)
      Dimension D1invt_bb(Nvrt_b,Nocc_b)
      Dimension D2invt_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension D2invt_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension D2invt_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)
      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)
      
      Integer I,J,A,B

#ifdef _DEBUG_LVLM
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

C AA block

      DO I = 1, Nocc_a
      DO A = 1, Nvrt_a
         Daa = D1hbar_aa(A,I)
         If (Dabs(Daa) .Gt. Denom_tol) Then
             D1invt_aa(A,I) = T1aa(A,I)/Daa
         ELse
            If (Regular) Then
                 Daa2 = Daa*Daa
                 Daa = (Daa2 + Rfac)/Daa
                 D1invt_aa(A,I) = T1aa(A,I)/Daa
             Else
                 D1invt_aa(A,I) = 0.0D0
             Endif
             D1invt_aa(A,I) = 0.0D0
         Endif
      ENDDO
      ENDDO

C BB block

      DO I = 1, Nocc_b
      DO A = 1, Nvrt_b
         Dbb = D1hbar_bb(a,i) 
         If (Dabs(Dbb) .Gt. Denom_tol) Then
            D1invt_bb(a,i) = T1bb(a,i)/Dbb
         Else
            If (Regular) Then
                 Dbb2 = Dbb*Dbb
                 Dbb = (Dbb2 + Rfac)/Dbb
                 D1invt_bb(A,I) = T1bb(A,I)/Dbb
             Else
                 D1invt_bb(A,I) = 0.0D0
             Endif
             D1invt_bb(A,I) = 0.0D0
         Endif
      ENDDO
      ENDDO

C AAAA block (also taking account the antisymmetry)

      DO J = 1, Nocc_a
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_a
      DO A = 1, Nvrt_a
         If (A .NE. B .AND. I .NE. J) Then
         Daaaa = D2hbar_aa(A,B,I,J)
         If (Dabs(Daaaa) .Gt. Denom_tol) Then
             D2invt_aa(A,B,I,J) = T2aa(A,B,I,J)/Daaaa
         ELse
            If (Regular) Then
                Daaaa2 = Daaaa*Daaaa
                Daaaa  = (Daaaa2 + Rfac)/Daaaa
                D2invt_aa(A,B,I,J) = T2aa(A,B,I,J)/Daaaa
            Else
               D2invt_aa(A,B,I,J) = 0.0D0
            Endif
         Endif

         Else
            D2invt_aa(A,B,I,J) = 0.0D0
         Endif 
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C BBBB block

      DO J = 1, Nocc_b
      DO I = 1, Nocc_b
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_b
         If (A .NE. B .AND. I .NE. J) Then
         Dbbbb = D2hbar_bb(a,b,i,j)
         If (Dabs(Dbbbb) .Gt. Denom_tol) Then
            D2invt_bb(a,b,i,j) = T2bb(a,b,i,j)/Dbbbb
         Else
            If (Regular) Then
                Dbbbb2 = Dbbbb*Dbbbb
                Dbbbb  = (Dbbbb2 + Rfac)/Dbbbb
                D2invt_bb(a,b,i,j) = T2bb(a,b,i,j)/Dbbbb
            Else
                D2invt_bb(a,b,i,j) = 0.0D0
            Endif
         Endif
         Else 
            D2invt_bb(a,b,i,j) = 0.0D0
         Endif 
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C ABAB block

      DO J = 1, Nocc_b
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_a
         Dabab = D2hbar_ab(A,b,I,j)
         If (Dabs(Dabab) .Gt. Denom_tol) Then
            D2invt_ab(A,b,I,j) = T2ab(A,b,I,j)/Dabab
         Else
            If (Regular) Then
                Dabab2 = Dabab*Dabab
                Dabab  = (Dabab2 + Rfac)/Dabab
                D2invt_ab(A,b,I,j) = T2ab(A,b,I,j)/Dabab
            Else
                D2invt_ab(A,b,I,j) = 0.0D0
            Endif
            D2invt_ab(A,b,I,j) = 0.0D0
         Endif
      ENDDO
      ENDDO
      ENDDO
      ENDDO

#ifdef _DEBUG_LVL0
      L_aa = Nocc_a*Nvrt_a
      L_bb = Nocc_b*Nvrt_b
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("D1invt_aa :",D1invt_aa,L_aa)
      call checksum("D1invt_bb :",D1invt_bb,L_bb)
      call checksum("D2invt_aa :",D2invt_aa,L_aaaa)
      call checksum("D2invt_bb :",D2invt_bb,L_bbbb)
      call checksum("D2invt_ab :",D2invt_ab,L_abab)
#endif

      Return
      End 

