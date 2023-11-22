      Subroutine T1_int1(T1aa,T1bb,T1resid_aa,T1resid_bb,W3_aaaa,
     +                   W3_bbbb,W3_abab,W3_baba,W3_abba,W3_baab,
     +                   Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)
 
      Dimension T1resid_aa(Nvrt_a,Nocc_a)
      Dimension T1resid_bb(Nvrt_b,Nocc_b)

      Dimension W3_aaaa(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W3_bbbb(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W3_abab(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W3_baba(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension W3_abba(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension W3_baab(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Integer M,N,A,F,E,I

#include "ccsdlight_vars.com"

C AA block 
C T1resid(A,I)= +T1(F,N)*W(N,A,F,I)
      
      Do i = 1, Nocc_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do n = 1, Nocc_a
      Do f = 1, Nvrt_a
CSSS         C = (1.0D0-Ocn_va(F))*Ocn_oa(N)
         C = (1.0D0-Ocn_va(A))*Ocn_oa(I)
         T = T + T1aa(F,N)*W3_aaaa(N,A,F,I)*C
      Enddo
      Enddo
      T1resid_aa(A,I) = T1resid_aa(A,I) + T
      Enddo 
      Enddo 

C T1resid(A,I)= +T1(f,n)*W(n,A,f,I)

      Do i = 1, Nocc_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do n = 1, Nocc_b
      Do f = 1, Nvrt_b
CSSS         C = (1.0D0-Ocn_vb(f))*Ocn_ob(n)
         C = (1.0D0-Ocn_va(A))*Ocn_oa(I)
         T = T + T1bb(f,n)*W3_baba(n,A,f,I)*C
      Enddo
      Enddo
      T1resid_aa(A,I) = T1resid_aa(A,I) + T
      Enddo
      Enddo

C BB block
C T1resid(a,i)= +T1(f,n)*W(n,a,f,i)

      Do i = 1, Nocc_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do n = 1, Nocc_b
      Do f = 1, Nvrt_b
CSSS         C = (1.0D0-Ocn_vb(f))*Ocn_ob(n)
         C = (1.0D0-Ocn_vb(a))*Ocn_ob(i)
         T = T + T1bb(f,n)*W3_bbbb(n,a,f,i)*C
      Enddo
      Enddo
      T1resid_bb(A,I) = T1resid_bb(A,I) + T
      Enddo
      Enddo

C T1resid(a,i)= +T1(F,N)*W(N,a,F,i)

      Do i = 1, Nocc_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do n = 1, Nocc_a
      Do f = 1, Nvrt_a
CSSS         C = (1.0D0-Ocn_va(F))*Ocn_oa(N)
         C = (1.0D0-Ocn_vb(a))*Ocn_ob(i)
         T = T + T1aa(F,N)*W3_abab(N,a,F,i)*C
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
