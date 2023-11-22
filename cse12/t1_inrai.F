      Subroutine T1_inrai(T1aa,T1bb,T1resid_aa,T1resid_bb,Fockoo_a,
     +                    Fockoo_b,Fockvv_a,Fockvv_b,Nocc_a,Nocc_b,
     +                    Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension Fockoo_a(Nocc_a,Nocc_a)
      Dimension Fockoo_b(Nocc_b,Nocc_b)
      Dimension Fockvv_a(Nvrt_a,Nvrt_a)
      Dimension Fockvv_b(Nvrt_b,Nvrt_b)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension T1resid_aa(Nvrt_a,Nocc_a)
      Dimension T1resid_bb(Nvrt_b,Nocc_b)

      Integer A,B,I,K

#include "cse12_vars.com"

C Rai_a(A,I) <- -T1(A,K)*Fockoo_a(K,I)

      DO I = 1, Nocc_a
      DO A = 1, Nvrt_a
         T = 0.0D0
         C = Ocn_oa(I)
      DO K = 1, Nocc_a
         If (K .EQ. I) Then
             T = T + T1aa(A,K)*Fockoo_a(K,I)*C
         Endif
      ENDDO
         T1resid_aa(A,I) = T1resid_aa(A,I) - T
      ENDDO
      ENDDO

C Rai_a(A,I) <- +T1(B,I)*Fockvv_a(B,A)

      DO I = 1, Nocc_a
      DO A = 1, Nvrt_a
         T = 0.0D0
         C = (1.0D0-Ocn_va(A))
      DO B = 1, Nvrt_a
         If (A .EQ. B) Then
             T = T + T1aa(B,I)*Fockvv_a(B,A)*C
         Endif
      ENDDO
         T1resid_aa(A,I) = T1resid_aa(A,I) + T
      ENDDO
      ENDDO

C Rai_b(a,i) <- -T1(a,k)*Fockoo_b(k,i)

      DO I = 1, Nocc_b
      DO A = 1, Nvrt_b
         T = 0.0D0
         C = Ocn_ob(i)
      DO K = 1, Nocc_b
         If (k .EQ. i) Then
             T = T + T1bb(a,k)*Fockoo_b(k,i)*C
         Endif
      ENDDO
         T1resid_bb(a,i) = T1resid_bb(a,i) - T
      ENDDO
      ENDDO

C Rai_b(a,i) <- +T1(b,i)*Fockvv_a(b,a)

      DO I = 1, Nocc_b
      DO A = 1, Nvrt_b
         T = 0.0D0
         C = (1.0D0-Ocn_vb(a))
      DO B = 1, Nvrt_b
         If (a .EQ. b) Then
             T = T + T1aa(b,i)*Fockvv_b(b,a)*C
         Endif
      ENDDO
         T1resid_bb(a,i) = T1resid_bb(a,i) + T
      ENDDO
      ENDDO

#ifdef _DEBUG_LVL0
      L_aa = Nocc_a*Nvrt_a
      L_bb = Nocc_b*Nvrt_b
      call checksum("t1_inrai  :",T1resid_aa,L_aa)
      call checksum("t1_inrai  :",T1resid_bb,L_bb)
#endif
      Return 
      End

