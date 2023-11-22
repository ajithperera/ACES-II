      subroutine Built_frac_fock_mo(W0_aa,W0_bb,W0_ab,Fock_aa,
     +                              Fock_bb,Nocc_a,Nocc_b,Nvrt_a,
     +                              Nvrt_b,Nbasis)

      Implicit Double Precision(A-H,O-Z)

      Dimension W0_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension w0_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension W0_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Fock_aa(Nbasis,Nbasis) 
      Dimension Fock_bb(Nbasis,Nbasis) 

      Integer P,Q,R,S

#include "cse_vars.com"

      Data Half /0.50D0/

C Note that for fractional occupations Nocc_a,Nvrta_a,Nocc_b,Nvrt_b = Nbasis.
C It is a matter of convenince that 

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "@-built_frac_fock: one Hamiltonian matrices"
      call output(Fock_aa,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
      call output(Fock_bb,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
#endif 

C built AA MO fock matrix

      Do Q = 1, Nocc_a
      Do P = 1, Nocc_a
         T = 0.0D0
         DO R = 1, Nocc_a 
         DO S = 1, Nocc_a 
            If (R. EQ. S) Then
               C = Ocn_oa(R)
               T = T + W0_aa(P,S,Q,R)*C*Half 
            Endif
         Enddo
         Enddo 
         Fock_aa(P,Q) = Fock_aa(P,Q) + T
      Enddo
      Enddo 

      Do Q = 1, Nocc_a
      Do P = 1, Nocc_a
         T = 0.0D0
         DO R = 1, Nocc_a
         DO S = 1, Nocc_a
            If (R. EQ. S) Then
               C = Ocn_oa(R)
               T = T + W0_aa(R,P,S,Q)*C*Half
            Endif
         Enddo
         Enddo
         Fock_aa(P,Q) = Fock_aa(P,Q) + T
      Enddo
      Enddo

      Do Q = 1, Nocc_a
      Do P = 1, Nocc_a
         T = 0.0D0
         DO R = 1, Nocc_b
         DO S = 1, Nocc_b
            If (R. EQ. S) Then
               C = Ocn_ob(r)
               T = T + W0_ab(P,s,Q,r)*C 
            Endif
         Enddo
         Enddo 
         Fock_aa(P,Q) = Fock_aa(P,Q) + T 
      Enddo
      Enddo 

C Built BB MO fock matrix 

      Do Q = 1, Nocc_b
      Do P = 1, Nocc_b
         T = 0.0D0
         DO R = 1, Nocc_b
         DO S = 1, Nocc_b
            If (r. EQ. s) Then
               C = Ocn_ob(r)
               T = T + W0_bb(p,s,q,r)*C*Half 
            Endif 
         Enddo
         Enddo 
         Fock_bb(p,q) = Fock_bb(p,q) + T
      Enddo
      Enddo 

      Do Q = 1, Nocc_b
      Do P = 1, Nocc_b
         T = 0.0D0
         DO R = 1, Nocc_b
         DO S = 1, Nocc_b
            If (r. EQ. s) Then
               C = Ocn_ob(r)
               T = T + W0_bb(s,p,r,q)*C*Half
            Endif 
         Enddo
         Enddo
         Fock_bb(p,q) = Fock_bb(p,q) + T
      Enddo
      Enddo

      Do Q = 1, Nocc_b
      Do P = 1, Nocc_b
         T = 0.0D0
         DO R = 1, Nocc_a
         DO S = 1, Nocc_a
            If (r. EQ. s) Then
               C = Ocn_oa(R)
               T = T + W0_ab(R,p,S,q)*C
            Endif
         Enddo
         Enddo
         Fock_bb(p,q) = Fock_bb(p,q) + T
      Enddo
      Enddo

#ifdef _DEBUG_LVLM
      Write(6,"(a)") "@-built_frac_fock: MO basis fock matrices"
      call output(Fock_aa,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
      call output(Fock_bb,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
#endif 
      Return
      End 


