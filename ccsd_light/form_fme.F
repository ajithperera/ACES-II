      Subroutine Form_fme(Fme_a,Fme_b,T1aa,T1bb,W_aa,W_bb,W_ab,
     +                    Fockov_a,Fockov_b,Nocc_a,Nocc_b,Nvrt_a,
     +                    Nvrt_b,Fme_on)

      Implicit Double Precision(A-H,O-Z)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension Fme_a(Nocc_a,Nvrt_a),Fme_b(Nocc_b,Nvrt_b)
      Dimension Fockov_a(Nocc_a,Nvrt_a),Fockov_b(Nocc_b,Nvrt_b)

      Integer M,E,N,F
      Logical Fme_on

#include "ccsdlight_vars.com"
      
#ifdef _DEBUG_LVLM
      call checksum("T1aa      :",T1aa,Nocc_a*Nvrt_a)
      call checksum("T1bb      :",T1bb,Nocc_b*Nvrt_b)
#endif 

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

C Fme(M,E) = T1(N,F)*W(MN,EF) + T1(n,f)*W(M,n,E,f)

      DO E = 1, Nvrt_a
      Do M = 1, Nocc_a
         T = 0.0D0
      Do N = 1, Nocc_a
      Do F = 1, Nvrt_a 
CSSS         C = Ocn_oa(N)*(1.0D0-Ocn_va(F))
         C = 1.0D0
         T = T + T1aa(F,N)*W_aa(M,N,E,F)*C
      Enddo
      Enddo
         Fme_a(M,E) = Fme_a(M,E) + T
      Enddo
      Enddo

      DO E = 1, Nvrt_a
      Do M = 1, Nocc_a
         T = 0.0D0
      Do N = 1, Nocc_b
      Do F = 1, Nvrt_b
CSSS         C = Ocn_ob(n)*(1.0D0-Ocn_vb(f))
         C = 1.0D0
         T = T + T1bb(f,n)*W_ab(M,n,E,f)*C 
      Enddo
      Enddo
         Fme_a(M,E) = Fme_a(M,E) + T
      Enddo
      Enddo

C Fme(e,m) = T1(n,f)*W(m,n,e,f) + T1(N,F)*W(m,N,e,F)

      DO E = 1, Nvrt_b
      Do M = 1, Nocc_b
         T = 0.0D0
      Do N = 1, Nocc_b
      Do F = 1, Nvrt_b
CSSS         C = Ocn_ob(n)*(1.0D0-Ocn_vb(f))
         C = 1.0D0
         T = T + T1bb(f,n)*W_bb(m,n,e,f)*C
      Enddo
      Enddo
         Fme_b(m,e) = Fme_b(m,e) + T
      Enddo
      Enddo

      DO E = 1, Nvrt_b
      Do M = 1, Nocc_b
         T = 0.0D0
      Do N = 1, Nocc_a
      Do F = 1, Nvrt_a
CSSS         C = Ocn_oa(N)*(1.0D0-Ocn_va(F))
         C = 1.0D0
         T =  T + T1aa(F,N)*W_ab(N,m,F,e)*C 
      Enddo
      Enddo
         Fme_b(m,e) = Fme_b(m,e) + T
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      call checksum("Fme_a     :",Fme_a,Nocc_a*Nvrt_a)
      call checksum("Fme_b     :",Fme_b,Nocc_b*Nvrt_b)
#endif 
       
      Return 
      End
