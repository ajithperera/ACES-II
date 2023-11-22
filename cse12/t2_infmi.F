      Subroutine T2_infmi(T2aa,T2bb,T2ab,Fmi_a,Fmi_b,W_aa,W_bb,
     +                    W_ab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_b) 
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b) 
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b) 

      Dimension Fmi_a(Nocc_a,Nocc_a)
      Dimension Fmi_b(Nocc_b,Nocc_b)

      Integer M,N,E,F,I

#include "cse12_vars.com"

C Fmi_aa(M,I) <- 1/2 T2(EF,IN)*W2(MN,EF) + T2(Ef,In)*W(Mn,Ef)
      
      Do I = 1, Nocc_a
      Do M = 1, Nocc_a
         T = 0.0D0 
      Do E = 1, Nvrt_a
      Do F = 1, Nvrt_a
      Do N = 1, Nocc_a
         T = T + 0.50D0*T2aa(E,F,I,N)*W_aa(M,N,E,F)
      Enddo
      Enddo
      Enddo 
      Fmi_a(M,I) = T 
      Enddo 
      Enddo 
 
      Do I = 1, Nocc_a
      Do M = 1, Nocc_a
         T = 0.0D0 
      Do E = 1, Nvrt_a
      Do f = 1, Nvrt_b
      Do n = 1, Nocc_b
         T = T + T2ab(E,f,I,n)*W_ab(M,n,E,f)
      Enddo
      Enddo
      Enddo 
      Fmi_a(M,I) = Fmi_a(M,I) + T
      Enddo 
      Enddo 

C Fmi_bb(m,i) <- 1/2 T2(ef,in)*W2(mn,ef) + T2(eF,iN)*W(mN,eF)

      Do i = 1, Nocc_b
      Do m = 1, Nocc_b
         T = 0.0D0
      Do e = 1, Nvrt_b
      Do f = 1, Nvrt_b
      Do n = 1, Nocc_b
         T = T + 0.50D0*T2bb(e,f,i,n)*W_bb(m,n,e,f)
      Enddo
      Enddo
      Enddo
      Fmi_b(m,i) = T 
      Enddo
      Enddo 

      Do i = 1, Nocc_b
      Do m = 1, Nocc_b
         T = 0.0D0
      Do e = 1, Nvrt_b
      Do F = 1, Nvrt_a
      Do N = 1, Nocc_a
         T = T + T2ab(F,e,N,i)*W_ab(N,m,F,e)
      Enddo
      Enddo
      Enddo
      Fmi_b(m,i) = Fmi_b(m,i) + T 
      Enddo
      Enddo 

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " @-T2_infmi" 
      call checksum("T2->Fmi_a :",Fmi_a,Nocc_a*Nocc_a)
      call checksum("T2->Fmi_b :",Fmi_b,Nocc_b*Nocc_b)
#endif 

      Return 
      End 
