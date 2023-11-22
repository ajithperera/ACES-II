      Subroutine T2abij_inrabij(T2aa,T2bb,T2ab,T2resid_aa,T2resid_bb,
     +                          T2resid_ab,fockoo_a,fockoo_b,fockvv_a,
     +                          fockvv_b,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)
      
      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension T2resid_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2resid_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2resid_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension fockoo_a(Nocc_a,Nocc_a)
      Dimension fockoo_b(Nocc_b,Nocc_b)
      Dimension fockvv_a(Nvrt_a,Nvrt_a)
      Dimension fockvv_b(Nvrt_b,Nvrt_b)

      Integer I,J,A,B,K,E

#include "cse12_vars.com"
      
C AAAA block 

      DO J = 1, Nocc_a
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_a
      DO A = 1, Nvrt_a
         T = 0.0D0
         C = Ocn_oa(J)
      DO K = 1, Nocc_a
         IF (K .EQ. J) Then
             T = T + T2aa(A,B,I,K)*Fockoo_a(K,J)*C 
         Endif 
      ENDDO 
         T2resid_aa(A,B,I,J) =  T2resid_aa(A,B,I,J) - T
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO J = 1, Nocc_a
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_a
      DO A = 1, Nvrt_a
         T = 0.0D0
         C = Ocn_oa(I)
      DO K = 1, Nocc_a
         IF (K .EQ. I) Then
             T = T + T2aa(A,B,K,J)*Fockoo_a(K,I)*C
         Endif 
      ENDDO 
         T2resid_aa(A,B,I,J) =  T2resid_aa(A,B,I,J) - T
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO J = 1, Nocc_a
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_a
      DO A = 1, Nvrt_a
         T = 0.0D0
         C = Ocn_va(A)
      DO E = 1, Nvrt_a
         IF (E .EQ. A) Then
             T = T + T2aa(E,B,I,J)*Fockvv_a(A,E)*C
         Endif
      ENDDO
         T2resid_aa(A,B,I,J) =  T2resid_aa(A,B,I,J) + T
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO J = 1, Nocc_a
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_a
      DO A = 1, Nvrt_a
         T = 0.0D0
         C = Ocn_va(B)
      DO E = 1, Nvrt_a
         IF (E .EQ. B) Then
             T = T + T2aa(A,E,K,J)*Fockvv_a(B,E)*C 
         Endif
      ENDDO
         T2resid_aa(A,B,I,J) =  T2resid_aa(A,B,I,J) + T
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C BBBB block 

      DO J = 1, Nocc_b
      DO I = 1, Nocc_b
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_b
         T = 0.0D0
         C = Ocn_ob(j)
      DO K = 1, Nocc_a
         IF (K .EQ. J) Then
             T = T + T2bb(a,b,i,k)*Fockoo_b(k,j)*C
         Endif
      ENDDO
         T2resid_bb(a,b,i,j) =  T2resid_bb(a,b,i,j) - T
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO J = 1, Nocc_b
      DO I = 1, Nocc_b
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_b
         T = 0.0D0
         C = Ocn_ob(i)
      DO K = 1, Nocc_a
         IF (K .EQ. I) Then
             T = T + T2bb(a,b,k,j)*Fockoo_b(k,i)*C
         Endif
      ENDDO
         T2resid_bb(a,b,i,j) =  T2resid_bb(a,b,i,j) - T
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO J = 1, Nocc_b
      DO I = 1, Nocc_b
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_b
         T = 0.0D0
         C = Ocn_vb(a)
      DO E = 1, Nvrt_b
         IF (E .EQ. A) Then
             T = T + T2aa(e,b,i,j)*Fockvv_b(a,e)*C
         Endif
      ENDDO
         T2resid_bb(a,b,i,j) =  T2resid_bb(a,b,i,j) + T
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO J = 1, Nocc_b
      DO I = 1, Nocc_b
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_b
         T = 0.0D0
         C = Ocn_vb(b)
      DO E = 1, Nvrt_b
         IF (E .EQ. B) Then
             T = T + T2bb(a,e,k,j)*Fockvv_b(b,e)*C 
         Endif
      ENDDO
         T2resid_bb(a,b,i,j) =  T2resid_bb(a,b,i,j) + T
      ENDDO
      ENDDO
      ENDDO
      ENDDO


C ABAB block 

      DO J = 1, Nocc_b
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_a
         T = 0.0D0
         C = Ocn_ob(j)
      DO K = 1, Nocc_b
         IF (K .EQ. J) Then
             T = T + T2ab(A,b,I,k)*Fockoo_b(k,j)*C
         Endif
      ENDDO
         T2resid_ab(A,b,I,j) =  T2resid_ab(A,b,i,j) - T
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO J = 1, Nocc_b
      DO I = 1, Nocc_b
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_b
         T = 0.0D0
         C = Ocn_oa(I)
      DO K = 1, Nocc_a
         IF (K .EQ. I) Then
             T = T + T2ab(A,b,K,j)*FocKoo_a(K,I)*C
         Endif
      ENDDO
         T2resid_ab(A,b,I,j) =  T2resid_ab(A,b,I,j) - T
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO J = 1, Nocc_b
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_a
         T = 0.0D0
         C = Ocn_va(a)
      DO E = 1, Nvrt_a
         IF (E .EQ. A) Then
             T = T + T2ab(E,b,I,j)*Fockvv_a(A,E)*C
         Endif
      ENDDO
         T2resid_ab(A,b,I,j) =  T2resid_ab(A,b,I,j) + T
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO J = 1, Nocc_b
      DO I = 1, Nocc_b
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_b
         T = 0.0D0
         C = Ocn_vb(b)
      DO E = 1, Nvrt_b
         IF (E .EQ. B) Then
             T = T + T2bb(A,e,I,j)*Fockvv_b(b,e)*C
         Endif
      ENDDO
         T2resid_ab(A,b,I,j) =  T2resid_ab(A,b,I,j) + T
      ENDDO
      ENDDO
      ENDDO
      ENDDO

#ifdef _DEBUG_LVL0
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("T2aa      :",T2aa_new,L_aaaa)
      call checksum("T2bb      :",T2bb_new,L_bbbb)
      call checksum("T2ab      :",T2ab_new,L_abab)
#endif

      Return
      End
