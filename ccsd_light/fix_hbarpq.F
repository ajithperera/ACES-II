      Subroutine Fix_Hbarpq(Hbar_mia,Hbar_mib,Hbar_aea,Hbar_aeb,
     +                      Hbar_mea,Hbar_meb,Fockoo_a,Fockoo_b,
     +                      Fockvv_a,Fockvv_b,Nocc_a,Nocc_b,
     +                      Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension Hbar_mia(Nocc_a,Nocc_a)
      Dimension Hbar_mib(Nocc_b,Nocc_b)
      Dimension Hbar_aea(Nvrt_a,Nvrt_a)
      Dimension Hbar_aeb(Nvrt_b,NvrT_b)
      Dimension Hbar_mea(Nocc_a,Nvrt_a)
      Dimension Hbar_meb(Nocc_b,Nvrt_b)

      Dimension Fockoo_a(Nocc_a,Nocc_a)
      Dimension Fockoo_b(Nocc_b,Nocc_b)
      Dimension Fockvv_a(Nvrt_a,Nvrt_a)
      Dimension Fockvv_b(Nvrt_b,Nvrt_b)

      Integer I,J,A,E

#include "ccsdlight_vars.com"

C Add the diagonal Fock elements to Hbar_mi and Hbar_ae intermediates.

      Do I = 1, Nocc_a
      Do M = 1, Nocc_a
         C = Ocn_oa(I)
         If (M .EQ. I) Then
            Hbar_mia(M,I) = Hbar_mia(M,I) +  Fockoo_a(M,I)*C 
         Endif 
      Enddo
      Enddo

      Do I = 1, Nocc_b
      Do M = 1, Nocc_b
         C = Ocn_ob(i)
         If (m .EQ. i) Then 
            Hbar_mib(m,i) = Hbar_mib(m,i) +  Fockoo_b(m,i)*C
         Endif 
      Enddo
      Enddo

      Do E = 1, Nvrt_a
      Do A = 1, Nvrt_a
         C = (1.0D0-Ocn_va(A))
         If (A .EQ. E) Then
            Hbar_aea(A,E) = Hbar_aea(A,E) + Fockvv_a(A,E)*C 
         Endif 
      Enddo
      Enddo

      Do E = 1, Nvrt_b
      Do A = 1, Nvrt_b
         C = (1.0D0-Ocn_vb(a))
         If (a .EQ. e) Then 
            Hbar_aeb(a,e) = Hbar_aeb(a,e) + Fockvv_b(a,e)*C 
         Endif 
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      call checksum("Hbar_mia  :",Hbar_mia,Nocc_a*Nocc_a)
      call checksum("Hbar_mib  :",Hbar_mib,Nocc_b*Nocc_b)
      call checksum("Hbar_aea  :",Hbar_aea,Nvrt_a*Nvrt_a)
      call checksum("Hbar_aeb  :",Hbar_aeb,Nvrt_b*Nvrt_b)
      call checksum("Hbar_mea  :",Hbar_mea,Nocc_a*Nvrt_a)
      call checksum("Hbar_meb  :",Hbar_meb,Nocc_b*Nvrt_b)
#endif
      
      Return
      End 


 
