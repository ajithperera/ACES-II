      Subroutine Modf_htabci(T1aa,T1bb,W4_aa,W4_bb,W4_ab,W4_ba,
     +                       W_aa,W_bb,W_ab,Nocc_a,Nocc_b,Nvrt_a,
     +                       Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension W4_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W4_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W4_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W4_ba(Nvrt_a,Nvrt_b,Nocc_a,Nvrt_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Integer M,N,A,E,F 

#include "cse12_vars.com"

      Do N = 1, Nocc_a
      Do A = 1, Nvrt_a
      Do F = 1, Nvrt_a
      Do E = 1, Nvrt_a 
         T = 0.0D0
      DO M = 1, Nocc_a 
         If ((1.0D0-Ocn_oa(A)) .NE. 0.0D0) Then
            C = 1.0D0/(1.0D0-Ocn_oa(A))
         Else
            C = 1.0D0
         Endif
         T = T + T1aa(A,M)*W_aa(M,N,E,F)*C
      Enddo 
         W4_aa(E,F,A,N) = W4_aa(E,F,A,N) + T
      Enddo 
      Enddo
      Enddo 
      Enddo
      
      Do N = 1, Nocc_b
      Do A = 1, Nvrt_b
      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_B
         T = 0.0D0
      DO M = 1, Nocc_b
         If ((1.0D0-Ocn_ob(a)) .NE. 0.0D0)  Then
            C = 1.0D0/(1.0D0-Ocn_ob(a))
         Else
            C = 1.0D0
         Endif
         T = T + T1bb(a,m)*W_bb(m,n,e,f)*C
      Enddo
         W4_bb(e,f,a,n) = W4_bb(e,f,a,n) + T
      Enddo
      Enddo
      Enddo
      Enddo

      Do N = 1, Nocc_b
      Do A = 1, Nvrt_a
      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_a
         T = 0.0D0
      DO M = 1, Nocc_a
         If ((1.0D0-Ocn_oa(A)) .NE. 0.0D0) Then
            C = 1.0D0/(1.0D0-Ocn_oa(A))
         Else
            C = 1.0D0
         Endif
         T = T + T1aa(A,M)*W_ab(M,n,E,f)*C
      Enddo
         W4_ab(E,f,A,n) = W4_ab(E,f,A,n) + T
      Enddo
      Enddo
      Enddo
      Enddo

      Do N = 1, Nocc_a
      Do A = 1, Nvrt_b
      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_a
         T = 0.0D0
      DO M = 1, Nocc_b
         If ((1.0D0-Ocn_ob(a)) .NE. 0.0D0) Then
            C = 1.0D0/(1.0D0-Ocn_ob(a))
         Else
            C = 1.0D0
         Endif
         T = T + T1bb(a,m)*W_ab(N,m,E,f)*C 
      Enddo
         W4_ba(E,f,N,a) = W4_ba(E,f,N,a) + T
      Enddo
      Enddo
      Enddo
      Enddo

      Return
      End 

