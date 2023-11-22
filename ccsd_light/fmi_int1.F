      Subroutine Fmi_int1(T1aa,T1bb,T1resid_aa,T1resid_bb,Fmi_a,
     +                    Fmi_b,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b) 

      Implicit Double Precision(A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)
 
      Dimension T1resid_aa(Nvrt_a,Nocc_a)
      Dimension T1resid_bb(Nvrt_b,Nocc_b)

      Dimension Fmi_a(Nocc_a,Nocc_a)
      Dimension Fmi_b(Nocc_b,Nocc_b)

      Integer A,B,I,J,M

#include "ccsdlight_vars.com"

C AA block 

      Do I = 1, Nocc_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a 
CSSS         C = Ocn_oa(M) 
         C = 1.0D0
         T = T  + Fmi_a(M,I)*T1aa(A,M)*C
      Enddo 
         T1resid_aa(A,I) = T1resid_aa(A,I) - T 
      Enddo 
      Enddo 

C BB block 

      Do I = 1, Nocc_b
      Do A = 1, Nvrt_b
         T = 0.0D0
      Do M = 1, Nocc_b
CSSS         C = Ocn_ob(m) 
         C = 1.0D0
         T = T + Fmi_b(m,i)*T1bb(a,m)*C
      Enddo
         T1resid_bb(a,i) = T1resid_bb(a,i) - T
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      L_aa = Nocc_a*Nvrt_a
      L_bb = Nocc_b*Nvrt_b
      call checksum("T1resid_aa:",T1resid_aa,L_aa)
      call checksum("T1resid_bb:",T1resid_bb,L_bb)
#endif
      Return
      End 
