      Subroutine Form_modf_f(Fae_a,Fae_b,Fmi_a,Fmi_b,Fme_a,Fme_b,
     +                       T1aa,T1bb,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,
     +                       Scale)

      Implicit Double Precision (A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension Fae_a(Nvrt_a,Nvrt_a)
      Dimension Fae_b(Nvrt_b,Nvrt_b)

      Dimension Fmi_a(Nocc_a,Nocc_a)
      Dimension Fmi_b(Nocc_b,Nocc_b)

      Dimension Fme_a(Nocc_a,Nvrt_a)
      Dimension Fme_b(Nocc_b,Nvrt_b)
   
      Integer M,N,A,F,E

#include "cse12_vars.com"

C Fae_aa(A,E)= -+1/2 T1(A,M)*F(m,e)
      
      Do e = 1, Nvrt_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do m = 1, Nocc_a
         T = T + 0.50D0*T1aa(A,M)*Fme_a(M,E)
      Enddo
      Fae_a(A,E) = Fae_a(A,E) - Scale*T
      Enddo 
      Enddo 

C Fae_aa(a,e)= +-1/2 T1(a,m)*F(m,e)

      Do e = 1, Nvrt_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do m = 1, Nocc_b
         T = T + 0.50D0*T1bb(a,m)*Fme_b(m,e)
      Enddo
         Fae_b(a,e) = Fae_b(a,e) - Scale*T
      Enddo 
      Enddo 

C Fmi_aa(M,I) = -+1/2T(E,I)*F(M,E)

      Do I = 1, Nocc_a
      Do M = 1, Nocc_a
         T = 0.0D0
      Do E = 1, Nvrt_a
         T = T + 0.50D0*T1aa(E,I)*Fme_a(M,E)
      Enddo
      Fmi_a(M,I) = Fmi_a(M,I) + Scale*T
      Enddo
      Enddo

C Fmi_bb(m,i) = -+1/2T(e,i)*F(m,e)

      Do i = 1, Nocc_b
      Do m = 1, Nocc_b
         T = 0.0D0
      Do e = 1, Nvrt_b
         T = T + 0.50D0*T1bb(e,i)*Fme_b(m,e)
      Enddo
      Fmi_b(m,i) = Fmi_b(m,i) + Scale*T
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " Form_mod_f"
      call checksum("Fae_a     :",Fae_a,Nvrt_a*Nvrt_a)
      call checksum("Fae_b     :",Fae_b,Nvrt_b*Nvrt_b)
      call checksum("Fmi_a     :",Fmi_a,Nocc_a*Nocc_a)
      call checksum("Fmi_b     :",Fmi_b,Nocc_b*Nocc_b)
      Write(6,*) 
#endif  

      Return 
      End 
