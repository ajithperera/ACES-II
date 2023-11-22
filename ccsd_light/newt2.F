      Subroutine Newt2(T2aa_new,T2bb_new,T2ab_new,T2resid_aa,T2resid_bb,
     +                 T2resid_ab,fockoo_a,fockoo_b,fockvv_a,fockvv_b,
     +                 Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)
      
      Dimension T2aa_new(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb_new(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab_new(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension T2resid_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2resid_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2resid_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension fockoo_a(Nocc_a,Nocc_a)
      Dimension fockoo_b(Nocc_b,Nocc_b)
      Dimension fockvv_a(Nvrt_a,Nvrt_a)
      Dimension fockvv_b(Nvrt_b,Nvrt_b)

      Integer I,J,A,B

#include "ccsdlight_vars.com"
      
C AAAA block 

      DO J = 1, Nocc_a
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_a
      DO A = 1, Nvrt_a
         E_a = (1.0D0-Ocn_va(A))
         E_b = (1.0D0-Ocn_va(B))
         E_i = Ocn_oa(I)
         E_j = Ocn_oa(J)
         C   = E_a*E_b*E_i*E_j
         Daaaa = (Fockoo_a(I,I)*E_i+Fockoo_a(J,J)*E_j-
     +            Fockvv_a(A,A)*E_a-Fockvv_a(B,B)*E_b)
         If (Dabs(Daaaa) .Gt. Denom_tol) Then  
             T2aa_new(A,B,I,J) = T2resid_aa(A,B,I,J)/Daaaa
         ELse
            If (Regular) Then
                Daaaa2 = Daaaa*Daaaa
                Daaaa  = (Daaaa2 + Rfac)/Daaaa
                T2aa_new(A,B,I,J) = T2resid_aa(A,B,I,J)/Daaaa
            Else
                T2aa_new(A,B,I,J) = 0.0D0
            Endif 
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
         E_a = (1.0D0-Ocn_vb(a))
         E_b = (1.0D0-Ocn_vb(b))
         E_i = Ocn_ob(i)
         E_j = Ocn_ob(j)
         C   = E_a*E_b*E_i*E_j
         Dbbbb = (Fockoo_b(i,i)*E_i+Fockoo_b(j,j)*E_j-
     +            Fockvv_b(a,a)*E_a-Fockvv_b(b,b)*E_b)
         If (Dabs(Dbbbb) .Gt. Denom_tol) Then  
            T2bb_new(a,b,i,j) = T2resid_bb(a,b,i,j)/Dbbbb
         Else 
            If (Regular) Then
                Dbbbb2 = Dbbbb*Dbbbb
                Dbbbb  = (Dbbbb2 + Rfac)/Dbbbb
                T2bb_new(a,b,i,j) = T2resid_bb(a,b,i,j)/Dbbbb
            Else
                T2bb_new(a,b,i,j) = 0.0D0
            Endif 
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
         E_a = (1.0D0-Ocn_va(A))
         E_b = (1.0D0-Ocn_vb(b))
         E_i = Ocn_oa(I)
         E_j = Ocn_ob(j) 
         C   = E_a*E_b*E_i*E_j
         Dabab = (Fockoo_a(I,I)*E_i+Fockoo_b(j,j)*E_j-
     +            Fockvv_a(A,A)*E_a-Fockvv_b(b,b)*E_b)
         If (Dabs(Dabab) .Gt. Denom_tol) Then  
            T2ab_new(A,b,I,j) = T2resid_ab(A,b,I,j)/Dabab
#ifdef _DEBUG_LVLM
            Write(6,"(2a,4(1x,i2),2(1x,F15.10))") "A,B,I,J,Dabab,",
     +                                           "T2resid,T2ab_new: ",
     +                                            A,B,I,J,
     +                                            Dabab,
     +                                            T2resid_ab(A,b,I,j)
#endif 
         Else 
            If (Regular) Then
                Dabab2 = Dabab*Dabab
                Dabab  = (Dabab2 + Rfac)/Dabab
                T2ab_new(A,b,I,j) = T2resid_ab(A,b,I,j)/Dabab
            Else
               T2ab_new(A,b,I,j) = 0.0D0
            Endif 
         Endif 
           
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
