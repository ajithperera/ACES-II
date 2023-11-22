      Subroutine T2_infae(T2aa,T2bb,T2ab,Fae_a,Fae_b,W_aa,W_bb, 
     +                    W_ab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a) 
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b) 
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b) 

      Dimension Fae_a(Nvrt_a,Nvrt_a)
      Dimension Fae_b(Nvrt_b,Nvrt_b)
 
      Integer M,N,A,F,E

#include "cse12_vars.com"

C Fae_aa(A,E)<- -1/2 T2(AF,MN)*W(MN,EF) - T2(Af,Mn)*W(Mn,Ef)
      
      Do e = 1, Nvrt_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do m = 1, Nocc_a
      Do n = 1, Nocc_a
      Do f = 1, Nvrt_a
         T = T - 0.50D0*T2aa(A,F,N,M)*W_aa(N,M,E,F)
      Enddo
      Enddo
      Enddo 
      Fae_a(A,E) = T 
      Enddo 
      Enddo 

      Do e = 1, Nvrt_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do m = 1, Nocc_b
      Do n = 1, Nocc_a
      Do f = 1, Nvrt_b
         T = T - T2ab(A,f,N,m)*W_ab(N,m,E,f)
      Enddo
      Enddo
      Enddo 
         Fae_a(A,E) = Fae_a(A,E) + T
      Enddo 
      Enddo 

C Fae_bb(a,e)<- -1/2 T2(af,mn)*W2(mn,ef) - T2(mN,aF)*W2(mN,eF)

      Do e = 1, Nvrt_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do m = 1, Nocc_b
      Do n = 1, Nocc_b
      Do f = 1, Nvrt_b
         T = T - 0.50D0*T2bb(a,f,n,m)*W_bb(n,m,e,f)
      Enddo
      Enddo
      Enddo
         Fae_b(a,e)  = T
      Enddo
      Enddo 

      Do e = 1, Nvrt_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do m = 1, Nocc_a
      Do n = 1, Nocc_b
      Do f = 1, Nvrt_a
         T = T - T2ab(F,a,M,n)*W_ab(M,n,F,e)
      Enddo
      Enddo
      Enddo
         Fae_b(a,e)  = Fae_b(a,e) + T
      Enddo
      Enddo 

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " @-T2_infae" 
      call checksum("T2->Fae_a :",Fae_a,Nvrt_a*Nvrt_a)
      call checksum("T2->Fae_b :",Fae_b,Nvrt_b*Nvrt_b)
#endif  

      Return 
      End 
