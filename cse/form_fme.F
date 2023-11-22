      Subroutine Form_fme(Fme_a,Fme_b,Fockov_a,Fockov_b,Nocc_a,
     +                    Nocc_b,Nvrt_a,Nvrt_b,Fme_on)

      Implicit Double Precision(A-H,O-Z)

      Dimension Fme_a(Nocc_a,Nvrt_a),Fme_b(Nocc_b,Nvrt_b)
      Dimension Fockov_a(Nocc_a,Nvrt_a),Fockov_b(Nocc_b,Nvrt_b)

      Integer M,E,N,F
      Logical Fme_on

#include "cse_vars.com"
      
      Suma = 0.0D0
      DO E = 1, Nvrt_a
      Do M = 1, Nocc_a
          Fme_a(M,E) = Fockov_a(M,E)
          Suma = Suma + Fme_a(M,E)
      Enddo
      Enddo

      Sumb = 0.0D0
      DO E = 1, Nvrt_b
      Do M = 1, Nocc_b
          Fme_b(m,e) = Fockov_b(m,e)
          Sumb = Sumb + Fme_b(m,e)
      Enddo
      Enddo

      Fme_on = Suma .Gt. 0.0D0 .OR. Sumb .Gt. 0.0D0 

#ifdef _DEBUG_LVL0
      call checksum("Fme_a     :",Fme_a,Nocc_a*Nvrt_a)
      call checksum("Fme_b     :",Fme_b,Nocc_b*Nvrt_b)
#endif 
       
      Return 
      End
