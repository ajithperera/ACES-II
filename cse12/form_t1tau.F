      Subroutine Form_t1tau(T1tau_aa,T1tau_bb,T1tau_ab,T1aa,T1bb,
     +                      Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)
      
      Dimension T1tau_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T1tau_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T1tau_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Integer I,J,A,B

#include "cse12_vars.com"

C AAAA block 

      DO J = 1, Nocc_a
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_a
      DO A = 1, Nvrt_a
         T1tau_aa(A,B,I,J) = (T1aa(A,I)*T1aa(B,J)  -
     +                        T1aa(B,I)*T1aa(A,J))
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C BBBB block 

      DO J = 1, Nocc_b
      DO I = 1, Nocc_b
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_b
         T1tau_bb(a,b,i,j) = (T1bb(a,i)*T1bb(b,j) - 
     +                        T1bb(b,i)*T1bb(a,j))
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C ABAB block 


      DO J = 1, Nocc_b
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_a
         T1tau_ab(A,b,I,J) = T1aa(A,I)*T1bb(b,j)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

#ifdef _DEBUG_LVL0
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("T1tau_aa  :",T1tau_aa,L_aaaa)
      call checksum("T1tau_bb  :",T1tau_bb,L_bbbb)
      call checksum("T1tau_ab  :",T1tau_ab,L_abab)
      Write(6,*) 
#endif

      Return
      End
