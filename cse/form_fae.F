      Subroutine Form_fae(T1aa,T1bb,Fae_a,Fae_b,Fockvv_od_a,
     +                    Fockvv_od_b,Fockov_a,Fockov_b,Nocc_a,
     +                    Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension Fae_a(Nvrt_a,Nvrt_a)
      Dimension Fae_b(Nvrt_b,Nvrt_b)

      Dimension Fockov_a(Nocc_a,Nvrt_a)
      Dimension Fockov_b(Nocc_b,Nvrt_b)
 
      Dimension Fockvv_od_a(Nvrt_a,Nvrt_a)
      Dimension Fockvv_od_b(Nvrt_b,Nvrt_b)
   
      Integer M,N,A,F,E

#include "cse_vars.com"

C Fae_a(A,E)  = -1/2f(M,E)T(A,M)

      sum = 0.0D0
      Do e = 1, Nvrt_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do m = 1, Nocc_a
         T = T + 0.50D0*T1aa(A,M)*Fockov_a(M,E)
      Enddo
         C = (1.0D0-Ocn_va(A))
         Fae_a(A,E) =  Fockvv_od_a(A,E)*C - T
      Enddo
      Enddo

C Fae_a(a,e)  = -1/2f(m,e)T(a,m)

      Do e = 1, Nvrt_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do m = 1, Nocc_b
         T = T + 0.50D0*T1bb(a,m)*Fockov_b(m,e)
      Enddo
         C = (1.0D0-Ocn_vb(a))
         Fae_b(a,e) = Fockvv_od_b(a,e)*C - T
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      call checksum("Fae_a     :",Fae_a,Nvrt_a*Nvrt_a)
      call checksum("Fae_b     :",Fae_b,Nvrt_b*Nvrt_b)
      Write(6,"(a)") "From F(a,e)"
      do e=1,Nvrt_a
      do a=1,Nvrt_a
      Write(6,"(2I2,1x,F15.10)") a,e,Fae_a(a,e)
      enddo
      enddo
#endif  

      Return 
      End 
