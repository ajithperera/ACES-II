      Subroutine Form_canonical_dens(Densp_aa,Densp_bb,Densm_aa,
     +                               Densm_bb,T1aa,T1bb,Nocc_a,Nocc_b,
     +                               Nvrt_a,Nvrt_b,Nbasis)

      Implicit Double Precision(A-H,O-Z)
 
      Dimension Densp_aa(Nbasis,Nbasis)
      Dimension Densp_bb(Nbasis,Nbasis)
      Dimension Densm_aa(Nbasis,Nbasis)
      Dimension Densm_bb(Nbasis,Nbasis)
 
      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)
 
      Integer P,Q,A,I,J,E,M
      Data Ione /1/

#include "cse_vars.com"

C Form D_aa and (1.0-D_aa) for canonical occupations 

      Do I = 1, Nocc_a
      Do J = 1, Nocc_a
         T = 0.0D0
      Do E = 1, Nvrt_a 
         T = T + T1aa(E,I)*T1aa(E,J)
      Enddo
      If (I .EQ. J) Then
          Densp_aa(I,J) = Ocn_oa(I) - T
          Densm_aa(I,J) = (1.0D0-Ocn_oa(I)) + T
      Else
         Densp_aa(I,J) = -T
         Densm_aa(I,J) = +T
      Endif
      Enddo
      Enddo

      Ibegin = Nocc_a + 1
      Do P = Ibegin,Nbasis 
      Do Q = Ibegin,Nbasis
         T = 0.0D0
      Do M = 1, Nocc_a
         T = T1aa(P,M)*T1aa(Q,M)
      Enddo
      If (P .EQ. Q) Then
         Densp_aa(Q,P) =  Ocn_va(I) + T
         Densm_aa(Q,P) =  (1.0D0-Ocn_va(I))  - T
      Else
         Densp_aa(Q,P) =  T
         Densm_aa(Q,P) = -T
      Endif
      Enddo
      Enddo

      Do P = 1, Nalpha 
      Do Q = Ibegin , Nbasis 
         Densp_aa(Q,P) =  T1aa(Q,P)
         Densp_aa(P,Q) =  T1aa(Q,P)
         Densm_aa(Q,P) = -T1aa(Q,P)
         Densm_aa(P,Q) = -T1aa(Q,P)
      Enddo
      Enddo

C Form D_bb and (1.0-D_bb) for canonical occupations

      Do I = 1, NOcc_b
      Do J = 1, Nocc_b
         T = 0.0D0
      Do E = 1, Nvrt_b
         T = T + T1bb(E,I)*T1bb(E,J)
      Enddo
      If (I .EQ. J) Then
          Densp_bb(I,J) = Ocn_ob(I) - T
          Densm_bb(I,J) = (1.0D0-Ocn_oa(I)) + T
      Else
         Densp_bb(I,J) = -T
         Densm_bb(I,J) = +T
      Endif
      Enddo
      Enddo

      Ibegin = Nocc_b+1
      Do P = Ibegin,Nbasis
      Do Q = Ibegin,Nbasis
         T = 0.0D0
      Do M = 1, Nocc_b
         T = T1bb(P,M)*T1bb(Q,M)
      Enddo
      If (P .EQ. Q) Then
         Densp_bb(Q,P) =  Ocn_vb(I)  + T
         Densm_bb(Q,P) =  (1.0D0-Ocn_vb(I)) - T
      Else
         Densp_bb(Q,P) =  T
         Densm_bb(Q,P) = -T
      Endif
      Enddo
      Enddo

      Do P = 1, Nalpha
      Do Q = Ibegin , Nbasis
         Densp_bb(Q,P) =  T1bb(Q,P)
         Densp_bb(P,Q) =  T1bb(Q,P)
         Densm_bb(Q,P) = -T1bb(Q,P)
         Densm_bb(P,Q) = -T1bb(Q,P)
      Enddo
      Enddo
 
      Return
      End
