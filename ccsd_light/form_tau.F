      Subroutine Form_tau(T2aa,T2bb,T2ab,Tau_aa,Tau_bb,Tau_ab,
     +                    Tap_aa,Tap_bb,Tap_ab,T1aa,T1bb,Nocc_a,
     +                    Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)
      
      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Tau_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Tau_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Tau_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Tap_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Tap_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Tap_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Integer I,J,A,B

#include "ccsdlight_vars.com"

C AAAA block 

      DO J = 1, Nocc_a
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_a
      DO A = 1, Nvrt_a
         Tau_aa(A,B,I,J) = T2aa(A,B,I,J) + T1aa(A,I)*T1aa(B,J)  -
     +                     T1aa(B,I)*T1aa(A,J)
         Tap_aa(A,B,I,J) = T2aa(A,B,I,J) + (T1aa(A,I)*T1aa(B,J)  -
     +                     T1aa(B,I)*T1aa(A,J))*0.50D0
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C BBBB block 

      DO J = 1, Nocc_b
      DO I = 1, Nocc_b
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_b
         Tau_bb(a,b,i,j) = T2bb(a,b,i,j) + T1bb(a,i)*T1bb(b,j) - 
     +                     T1bb(b,i)*T1bb(a,j)
         Tap_bb(a,b,i,j) = T2bb(a,b,i,j) + (T1bb(a,i)*T1bb(b,j) - 
     +                     T1bb(b,i)*T1bb(a,j))*0.50D0
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C ABAB block 


      DO J = 1, Nocc_b
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_a
         Tau_ab(A,b,I,J) = T2ab(A,b,I,j) + T1aa(A,I)*T1bb(b,j) 
         Tap_ab(A,b,I,J) = T2ab(A,b,I,j) + T1aa(A,I)*T1bb(b,j)*0.50D0
      ENDDO
      ENDDO
      ENDDO
      ENDDO

#ifdef _DEBUG_LVL0
      Write(6,*) 
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("Tau_aa    :",Tau_aa,L_aaaa)
      call checksum("Tau_bb    :",Tau_bb,L_bbbb)
      call checksum("Tau_ab    :",Tau_ab,L_abab)
      call checksum("Tap_aa    :",Tap_aa,L_aaaa)
      call checksum("Tap_bb    :",Tap_bb,L_bbbb)
      call checksum("Tap_ab    :",Tap_ab,L_abab)
#endif

      Return
      End
