      subroutine Form_fmi(T1aa,T1bb,Fmi_a,Fmi_b,Fockoo_od_a,
     +                    Fockoo_od_b,Fockov_a,Fockov_b,Nocc_a,
     +                    Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)
     
      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension Fockoo_od_a(Nocc_a,Nocc_a)
      Dimension Fockoo_od_b(Nocc_b,Nocc_b)

      Dimension Fmi_a(Nocc_a,Nocc_a)
      Dimension Fmi_b(Nocc_b,Nocc_b)

      Dimension Fockov_a(Nocc_a,Nvrt_a)
      Dimension Fockov_b(Nocc_b,Nvrt_b)

      Integer M,N,E,F,I

#include "cse_vars.com"

C  Fmi_a(M,I) =  +1/2f(M,E)*T1(E,I)

      Do i = 1, Nocc_a
      Do m = 1, Nocc_a
         T = 0.0D0
      Do e = 1, Nvrt_a
         T = T + 0.50D0*T1aa(E,I)*Fockov_a(M,E)
      Enddo
         C = Ocn_oa(I) 
         Fmi_a(m,i) = Fockoo_od_a(M,I)*C + T
      Enddo
      Enddo

C  Fmi_b(m,i) =  +1/2f(m,e)*T1(e,i)

      Do i = 1, Nocc_b
      Do m = 1, Nocc_b
         T = 0.0D0
      Do e = 1, Nvrt_b
         T = T + 0.50D0*T1bb(e,i)*Fockov_b(m,e)
      Enddo
         C = Ocn_ob(i)
         Fmi_b(m,i) = Fockoo_od_b(m,i)*C  + T
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      call checksum("Fmi_a     :",Fmi_a,Nocc_a*Nocc_a)
      call checksum("Fmi_b     :",Fmi_b,Nocc_b*Nocc_b)
      do m=1,Nocc_a
      do i=1,Nocc_a
      Write(6,"(2I2,1x,F15.10)") m,i,Fmi_a(m,i)
      enddo
      enddo
#endif 

      Return 
      End 
