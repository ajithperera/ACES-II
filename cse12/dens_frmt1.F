      Subroutine Dens_frmt1(Densp_aa,Densp_bb,Densm_aa,Densm_bb,
     +                      T1aa,T1bb,Oneh_aa,Oneh_bb,Ref_energy,
     +                      Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,Nbasis,
     +                      Frac_occ)

      Implicit Double Precision(A-H,O-Z)
      Double Precision Oneh_aa,Onehh_bb

      Dimension Densp_aa(Nbasis,Nbasis)
      Dimension Densp_bb(Nbasis,Nbasis)
      Dimension Densm_aa(Nbasis,Nbasis)
      Dimension Densm_bb(Nbasis,Nbasis)

      Dimension Oneh_aa(Nbasis,Nbasis)
      Dimension Oneh_bb(Nbasis,Nbasis)
 
      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)
 
      Integer P,Q,A,I,J,E,M
      Logical Frac_occ

      Data Ione /1/

#include "cse12_vars.com"

      If (Frac_occ) Then

      Do P = 1, Nbasis
      Do Q = 1, Nbasis
         If (P .EQ. Q) Then
            Densp_aa(Q,P) =  T1aa(Q,P) + Ocn_oa(P)
            Densm_aa(Q,P) = -T1aa(Q,P) + (1.0D0-Ocn_oa(P))
         Else 
            Densp_aa(Q,P) =  T1aa(Q,P)
            Densm_aa(Q,P) = -T1aa(Q,P) 
         Endif 
      Enddo
      Enddo 
      
      Do P = 1, Nbasis
      Do Q = 1, Nbasis
         If (p .EQ. q) Then
            Densp_bb(q,p) =  T1bb(q,p) + Ocn_ob(p)
            Densm_bb(Q,P) = -T1bb(Q,P) + (1.0D0-Ocn_ob(p))
         Else
            Densp_bb(q,p) =  T1bb(q,p)
            Densm_bb(q,p) = -T1bb(q,p)
         Endif
      Enddo
      Enddo
 
      Else

      Call form_canonical_dens(Densp_aa,Densp_bb,Densm_aa,Densm_bb,
     +                         T1aa,T1bb,Nocc_a,Nocc_b,Nvrt_a,
     +                         Nvrt_b,Nbasis)
      Endif 
      E1_aa = Ddot(Nbasis*Nbasis,Oneh_aa,1,Densp_aa,1)
      E1_bb = Ddot(Nbasis*Nbasis,Oneh_bb,1,Densp_bb,1)
      E1    = E1_aa+E1_bb+Ref_energy
      Write(6,"(a,(1X,F15.10))") "D(p,q)*h(p,q)        :",E1

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Alpha/beta 1-density cumulants"
      call output(Densp_aa,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
      call output(Densp_bb,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
#endif 

      Return
      End
