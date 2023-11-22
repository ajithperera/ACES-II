      Subroutine Form_t1rai(T1aa,T1bb,Fmi_aa,Fmi_bb,Fae_aa,Fae_bb,
     +                      T1resid_aa,T1resid_bb,Nocc_a,Nocc_b,   
     +                      Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension T1resid_aa(Nvrt_a,Nocc_a)
      Dimension T1resid_bb(Nvrt_b,Nocc_b)

      Dimension Fmi_aa(Nocc_a,Nocc_a),Fmi_bb(Nocc_b,Nocc_b)
      Dimension Fae_aa(Nvrt_a,Nvrt_a),Fae_bb(Nvrt_b,Nvrt_b)

      Integer A,I,B,E,J,F,K

#include "cse12_vars.com"

C Htem_aa(A,I) <- Fmi_aa(M,I)*T(A,M) 
      
      Do I = 1, Nocc_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
         T = T  + Fmi_aa(M,I)*T1aa(A,M)
      Enddo
         T1resid_aa(A,I) = T1resid_aa(A,I) - T
      Enddo
      Enddo

C Htem_bb(a,i) <- Fmi_bb(m,i)*T(a,m) 

      Do I = 1, Nocc_b
      Do A = 1, Nvrt_b
         T = 0.0D0
      Do M = 1, Nocc_b
         T = T + Fmi_bb(m,i)*T1bb(a,m)
      Enddo
         T1resid_bb(a,i) = T1resid_bb(a,i) - T
      Enddo
      Enddo

C Htem_aa(A,I) <- Fae_aa(A,E)*T(E,I)

      Do I = 1, Nocc_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do E = 1, Nvrt_a
         T = T + Fae_aa(A,E)*T1aa(E,I)
      Enddo
         T1resid_aa(A,I) = T1resid_aa(A,I) + T
      Enddo
      Enddo

C Htem_bb(a,i) <- Fae_bb(a,e)*T(e,i)

      Do I = 1, Nocc_b
      Do A = 1, Nvrt_b
         T = 0.0D0
      Do E = 1, Nvrt_b
         T = T + Fae_bb(a,e)*T1bb(e,i)
      Enddo
         T1resid_bb(a,i) = T1resid_bb(a,i) + T
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      L_aa = Nocc_a*Nvrt_a
      L_bb = Nocc_b*Nvrt_b
      call checksum("T1resid_aa:",T1resid_aa,L_aa)
      call checksum("T1resid_bb:",T1resid_bb,L_bb)
      Write(6,*)
#endif

      Return
      End 
   


