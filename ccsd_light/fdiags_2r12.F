      Subroutine Fdiags_2r12(T1aa,T1bb,T2aa,T2bb,T2ab,T1resid_aa,
     +                       T1resid_bb,T2resid_aa,T2resid_bb,
     +                       T2resid_ab,Fockoo_a,Fockoo_b,Fockvv_a,
     +                       Fockvv_b,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)
      
      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)
      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension T1resid_aa(Nvrt_a,Nocc_a)
      Dimension T1resid_bb(Nvrt_b,Nocc_b)
      Dimension T2resid_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2resid_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2resid_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Fockoo_a(Nocc_a,Nocc_a)
      Dimension Fockoo_b(Nocc_b,Nocc_b)
      Dimension Fockvv_a(Nvrt_a,Nvrt_a)
      Dimension Fockvv_b(Nvrt_b,Nvrt_b)

      Integer I,J,A,B,E

#include "ccsdlight_vars.com"

C AA block 

      DO I = 1, Nocc_a
      DO A = 1, Nvrt_a
         T = 0.0D0
         C = Ocn_oa(I)
      DO K = 1, Nocc_a
         If (K .EQ. I) Then
             T = T + T1aa(A,K)*Fockoo_a(K,I)
         Endif 
      ENDDO 
        T1resid_aa(A,I) = T1resid_aa(A,I) - T*C
      ENDDO
      ENDDO

      DO I = 1, Nocc_a
      DO A = 1, Nvrt_a
         T = 0.0D0
         C = (1.0D0-Ocn_va(A))
      DO E = 1, Nvrt_a
         If (A .EQ. E) Then
             T = T + T1aa(E,I)*Fockvv_a(E,A)
         Endif 
      ENDDO 
         T1resid_aa(A,I) = T1resid_aa(A,I) + T*C 
      ENDDO
      ENDDO

C BB block 

      DO I = 1, Nocc_b
      DO A = 1, Nvrt_b
         T = 0.0D0
         C = Ocn_ob(i)
      DO K = 1, Nocc_b
         If (k .EQ. i) Then
             T = T + T1bb(a,k)*Fockoo_b(k,i)
         Endif 
      ENDDO
         T1resid_bb(a,i) = T1resid_bb(a,i) - T*C
      ENDDO
      ENDDO

      DO I = 1, Nocc_b
      DO A = 1, Nvrt_b
         T = 0.0D0
         C = (1.0D0-Ocn_vb(a))
      DO E = 1, Nvrt_b
         If (A .EQ.  E) Then
             T = T + T1bb(e,i)*Fockvv_b(e,a)
         Endif 
      ENDDO
         T1resid_bb(a,i) = T1resid_bb(a,i) + T*C
      ENDDO
      ENDDO

C AAAA block

      DO J = 1, Nocc_a
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_a
      DO A = 1, Nvrt_a
         T = 0.0D0 
         C = Ocn_oa(J)
      DO K = 1, Nocc_a
         If (J .EQ. K) Then
            T = T + T2aa(A,B,I,K)*Fockoo_a(J,K)
         Endif 
      ENDDO
         T2resid_aa(A,B,I,J) = T2resid_aa(A,B,I,J) - T*C
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
         If (I .EQ. K) Then
            T = T + T2aa(A,B,K,J)*Fockoo_a(I,K)
         Endif
      ENDDO
         T2resid_aa(A,B,I,J) = T2resid_aa(A,B,I,J) - T*C
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO J = 1, Nocc_a
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_a
      DO A = 1, Nvrt_a
         T = 0.0D0
         C = (1.0D0 - Ocn_va(B))
      DO E = 1, Nvrt_a
         If (B .EQ. E) Then
            T = T + T2aa(A,E,I,J)*Fockvv_a(E,B)
         Endif
      ENDDO
         T2resid_aa(A,B,I,J) = T2resid_aa(A,B,I,J) + T*C
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO J = 1, Nocc_a
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_a
      DO A = 1, Nvrt_a
         T = 0.0D0
         C = (1.0D0 - Ocn_va(A))
      DO E = 1, Nvrt_a
         If (A .EQ. E) Then
            T = T + T2aa(E,B,I,J)*Fockvv_a(E,A)
         Endif
      ENDDO
         T2resid_aa(A,B,I,J) = T2resid_aa(A,B,I,J) + T*C
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C BBBB Block

      DO J = 1, Nocc_b
      DO I = 1, Nocc_b
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_b
         T = 0.0D0
         C =  Ocn_ob(j)
      DO K = 1, Nocc_b
         If (j .EQ. k) Then
            T = T + T2bb(a,b,i,k)*Fockoo_b(k,j)
         Endif
      ENDDO
         T2resid_bb(a,b,i,j) = T2resid_bb(a,b,i,j) - T*C
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO J = 1, Nocc_b
      DO I = 1, Nocc_b
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_b
         T = 0.0D0
         C =  Ocn_ob(i)
      DO K = 1, Nocc_b
         If (i .EQ. k) Then
            T = T + T2bb(a,b,k,j)*Fockoo_b(k,i)
         Endif
      ENDDO
         T2resid_bb(a,b,i,j) = T2resid_bb(a,b,i,j) - T*C
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO J = 1, Nocc_b
      DO I = 1, Nocc_b
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_b
         T = 0.0D0
         C = (1.0D0 - Ocn_vb(b))
      DO E = 1, Nvrt_b
         If (b .EQ. e) Then
            T = T + T2bb(a,e,i,j)*Fockvv_b(e,b)
         Endif
      ENDDO
         T2resid_bb(a,b,i,j) = T2resid_bb(a,b,i,j) + T*C
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO J = 1, Nocc_b
      DO I = 1, Nocc_b
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_b
         T = 0.0D0
         C = (1.0D0 - Ocn_vb(a))
      DO E = 1, Nvrt_b
         If (a .EQ. e) Then
            T = T + T2bb(e,b,i,j)*Fockvv_b(e,a)
         Endif
      ENDDO
         T2resid_bb(a,b,i,j) = T2resid_bb(a,b,i,j) + T*C
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
         If (j .EQ. k) Then
            T = T + T2ab(A,b,I,k)*Fockoo_b(k,j)
         Endif
      ENDDO
         T2resid_ab(A,b,I,j) = T2resid_ab(A,b,I,j) - T*C
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO J = 1, Nocc_b
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_a
         T = 0.0D0
         C = Ocn_oa(I)
      DO K = 1, Nocc_a 
         If (I .EQ. K) Then
            T = T + T2ab(A,b,K,j)*Fockoo_a(K,I)
         Endif
      ENDDO
         T2resid_ab(A,b,I,j) = T2resid_ab(A,b,I,j) - T*C
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO J = 1, Nocc_b
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_a
         T = 0.0D0
         C = (1.0D0 - Ocn_vb(b))
      DO E = 1, Nvrt_b
         If (b .EQ. e) Then
            T = T + T2ab(A,e,I,j)*Fockvv_b(e,b)
         Endif
      ENDDO
         T2resid_ab(A,b,I,j) = T2resid_ab(A,b,I,j) + T*C
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO J = 1, Nocc_b
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_a
         T = 0.0D0
         C = (1.0D0 - Ocn_va(A))
      DO E = 1, Nvrt_a
         If (A .EQ. E) Then
            T = T + T2ab(E,b,I,j)*Fockvv_a(E,A)
         Endif
      ENDDO
         T2resid_ab(A,b,I,j) = T2resid_ab(A,b,I,j) + T*C
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
      call checksum("T1resid_aa:",T1resid_aa,L_aa)
      call checksum("T1resid_bb:",T1resid_bb,L_bb)
      call checksum("T2resid_aa:",T2resid_aa,L_aaaa)
      call checksum("T2resid_bb:",T2resid_bb,L_bbbb)
      call checksum("T2resid_ab:",T2resid_ab,L_abab)
#endif

      Return
      End
