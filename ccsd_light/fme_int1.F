      Subroutine Fme_int1(T2aa,T2bb,T2ab,T1resid_aa,T1resid_bb, 
     +                    Fme_a,Fme_b,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b) 

      Implicit Double Precision(A-H,O-Z)
 
      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)
 
      Dimension T1resid_aa(Nvrt_a,Nocc_a)
      Dimension T1resid_bb(Nvrt_b,Nocc_b)

      Dimension Fme_a(Nocc_a,Nvrt_a)
      Dimension Fme_b(Nocc_b,Nvrt_b)

      Integer A,B,I,J,M

#include "ccsdlight_vars.com"

C AA block 
C T1resid_aa(A,I) = Fme_a(M,E)*T2aa(A,E,I,M) +Fme_b(m,e)*T(A,e,I,m)

      Do I = 1, Nocc_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do E = 1, Nvrt_a 
      Do M = 1, Nocc_a 
CSSS`         C = Ocn_oa(M)*(1.0D0-Ocn_va(E))
         C = 1.0D0
         T = T  + Fme_a(M,E)*T2aa(A,E,I,M)*C
      Enddo 
      Enddo 
         T1resid_aa(A,I) = T1resid_aa(A,I) + T 
      Enddo 
      Enddo 

      Do I = 1, Nocc_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do E = 1, Nvrt_b
      Do M = 1, Nocc_b
CSSS         C = Ocn_ob(m)*(1.0D0-Ocn_vb(e))
         C = 1.0D0
         T = T  + Fme_b(m,e)*T2ab(A,e,I,m)*C
      Enddo
      Enddo
         T1resid_aa(A,I) = T1resid_aa(A,I) + T
      Enddo
      Enddo

C BB block  
C T1resid_bb(a,i) = Fme_a(m,e)*T2aa(a,e,i,m) +Fme_a(M,E)*T2ab(E,a,M,i)

      Do I = 1, Nocc_b
      Do A = 1, Nvrt_b
         T = 0.0D0
      Do E = 1, Nvrt_b
      Do M = 1, Nocc_b
CSSS         C = Ocn_ob(m)*(1.0D0-Ocn_vb(e))
         C = 1.0D0
         T = T  + Fme_b(m,e)*T2bb(a,e,i,m)*C
      Enddo
      Enddo
         T1resid_bb(a,i) = T1resid_bb(a,i) + T
      Enddo
      Enddo

      Do I = 1, Nocc_b
      Do A = 1, Nvrt_b
         T = 0.0D0
      Do E = 1, Nvrt_a 
      Do M = 1, Nocc_a
CSSS         C = Ocn_oa(M)*(1.0D0-Ocn_va(E))
         C = 1.0D0
         T = T  + Fme_a(M,E)*T2ab(E,a,M,i)*C
      Enddo
      Enddo
         T1resid_bb(A,I) = T1resid_bb(A,I) + T
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
