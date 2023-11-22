      Subroutine Get_ref_energy(W0_aa,W0_bb,W0_ab,Fock_a,Fock_b,
     +                          E1_aa,E1_bb,Nbasis,Naobfns,Nocc_a,
     +                          Nocc_b,Nvrt_a,Nvrt_b,
     +                          Work,Memleft,Repuls,Ref_energy)

      Implicit Double Precision(A-H,O-Z)
      
      Dimension W0_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension w0_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension W0_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Fock_a(Nbasis,Nbasis),Fock_b(Nbasis,Nbasis)

      Dimension Work(Memleft)
 
      Integer A,B,I,J

#include "ccsdlight_vars.com"
      
C AAAA block 

      E2_aa = 0.0D0
      DO J = 1, Nocc_a
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_a 
      DO A = 1, Nvrt_a
         If (A .EQ. I .And. B .EQ. J) Then
             C     = Ocn_oa(I)*Ocn_oa(J)
             E2_aa = E2_aa + W0_aa(A,B,I,J)*C
         Endif 
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C BBBB block 

      E2_bb = 0.0D0
      DO J = 1, Nocc_b
      DO I = 1, Nocc_b
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_b
         If (A .EQ. I .And. B .EQ. J) Then
            C     =  Ocn_ob(I)*Ocn_ob(J) 
            E2_bb = E2_bb+ W0_bb(a,b,i,j)*C
         Endif 
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C ABAB block 

      E2_ab = 0.0D0
      DO J = 1, Nocc_b
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_a
         If (A .EQ. I .And. B .EQ. J) Then
            C     = Ocn_oa(I)*Ocn_ob(J) 
            E2_ab = E2_ab + W0_ab(A,b,I,j)*C
         Endif 
      ENDDO
      ENDDO
      ENDDO
      ENDDO

#ifdef _DEBUG_LVL0
      Write(6,*)
      E2_a = E2_aa*0.25
      E2_b = E2_bb*0.25
      Write(6,"(a,3(1X,F15.10))") "E1aa,E1bb            :",E1_aa,E1_bb
      Write(6,"(a,3(1X,F15.10))") "E2aa,E2bb,E2ab       :",E2_aa,E2_bb,
     +                             E2_ab
#endif 
      E_tot      = E1_aa + E1_bb + (E2_aa + E2_bb)*0.50D0 + E2_ab
      Ref_energy = E_tot + Repuls 
      
#ifdef _DEBUG_LVL0
      Write(6,"(a,3(1X,F15.10))") "Reference E_tot      :",Ref_energy  
#endif 
     
      Return 
      End

