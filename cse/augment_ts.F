      Subroutine Augment_ts(T1aa,T1bb,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,
     +                      Daa,Dbb)

      Implicit Double Precision(A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension Daa(Nvrt_a,Nocc_a)
      Dimension Dbb(Nvrt_b,Nocc_b)

#include "cse_vars.com"

      Integer I,J,A,B

      Laa = Nvrt_a*Nocc_a
      LBb = Nvrt_b*Nocc_b
      Call Dzero(Daa, Laa)
      Call Dzero(Daa, Lbb)

      Index = 0
      Do I = 1, Nocc_a
         Index = Index + 1
      Do A = 1, Nvrt_a
         Daa(I,I) = Ocn_oa(index)
      Enddo
      Enddo

      Index = 0
      Do I = 1, Nocc_b
         Index = Index + 1
      Do A = 1, Nvrt_b
         Dbb(I,I) = Ocn_ob(index)
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The Alpha and Beta density matrices" 
      Call output(Daa,1,Nocc_a,1,Nvrt_a,Nocc_a,Nvrt_a,1)
      Write(6,*)
      Call output(Dbb,1,Nocc_a,1,Nvrt_a,Nocc_a,Nvrt_a,1)
#endif 

      Return
      End
