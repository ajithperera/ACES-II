      Subroutine Form_rabij(T2aa_new,T2bb_new,T2ab_new,R2abij_aa,
     +                      R2abij_bb,R2abij_ab,Fockoo_a,fockoo_b,
     +                      Fockvv_a,fockvv_b,Nocc_a,Nocc_b,Nvrt_a,
     +                      Nvrt_b,Scale)

      Implicit Double Precision(A-H,O-Z)
      
      Dimension T2aa_new(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb_new(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab_new(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension R2abij_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension R2abij_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension R2abij_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension fockoo_a(Nocc_a,Nocc_a)
      Dimension fockoo_b(Nocc_b,Nocc_b)
      Dimension fockvv_a(Nvrt_a,Nvrt_a)
      Dimension fockvv_b(Nvrt_b,Nvrt_b)

      Integer I,J,A,B

#include "cse12_vars.com"
      
C AAAA block 

      DO J = 1, Nocc_a
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_a
      DO A = 1, Nvrt_a
         E_a = (1.0D0-Ocn_va(A))
         E_b = (1.0D0-Ocn_va(B))
         E_i = Ocn_oa(I)
         E_j = Ocn_oa(J)
         Daaaa = (Fockoo_a(I,I)*E_i+Fockoo_a(J,J)*E_j-
     +            Fockvv_a(A,A)*E_a-Fockvv_a(B,B)*E_b)
         R2abij_aa(A,B,I,J) = R2abij_aa(A,B,I,J) - 
     +                        Daaaa*T2aa_new(A,B,I,J)*Scale 
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
         Dbbbb = (Fockoo_b(i,i)*E_i+Fockoo_b(j,j)*E_j-
     +            Fockvv_b(a,a)*E_a-Fockvv_b(b,b)*E_b)
         R2abij_bb(a,b,i,j) = R2abij_bb(a,b,i,j) - 
     +                        Dbbbb*T2bb_new(a,b,i,j)*Scale 
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
         Dabab = (Fockoo_a(I,I)*E_i+Fockoo_b(j,j)*E_j-
     +            Fockvv_a(A,A)*E_a-Fockvv_b(b,b)*E_b)
         R2abij_ab(A,b,I,j) = R2abij_ab(A,b,I,j) - 
     +                        Dabab*T2ab_new(A,b,I,j)*Scale 
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      Return
      End
