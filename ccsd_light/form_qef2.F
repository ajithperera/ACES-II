      Subroutine Form_qef2(Qae_aa,Qae_bb,W_aa,W_bb,W_ab,T2aa,T2bb,
     +                     T2ab,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)
     
      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a) 
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b) 
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b) 

      Dimension Qae_aa(Nvrt_a,Nvrt_a)
      Dimension Qae_bb(Nvrt_b,Nvrt_b)

      Integer M,N,A,F,E

#include "ccsdlight_vars.com"

C Qae_aa(A,E)=  T2(AF,MN)*W(MN,EF) + T2(Af,Mn)*W(Mn,Ef)
      
      Do e = 1, Nvrt_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do m = 1, Nocc_a
      Do n = 1, Nocc_a
      Do f = 1, Nvrt_a
         T = T + 0.50D0*T2aa(A,F,N,M)*W_aa(N,M,E,F)
      Enddo
      Enddo
      Enddo 
      Qae_aa(A,E) = -T 
      Enddo 
      Enddo 

      Do e = 1, Nvrt_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do m = 1, Nocc_b
      Do n = 1, Nocc_a
      Do f = 1, Nvrt_b
         C = 1.0D0
         T = T + T2ab(A,f,N,m)*W_ab(N,m,E,f)
      Enddo
      Enddo
      Enddo 
         Qae_aa(A,E) = Qae_aa(A,E) - T
      Enddo 
      Enddo 

C Fae_bb(a,e)= +1/2 T2(af,mn)*W(mn,ef) + T2(mN,aF)*W(mN,eF)

      Do e = 1, Nvrt_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do m = 1, Nocc_b
      Do n = 1, Nocc_b
      Do f = 1, Nvrt_b
         T = T + 0.50D0*T2bb(a,f,n,m)*W_bb(n,m,e,f)*C
      Enddo
      Enddo
      Enddo
      Qae_bb(a,e)  = -T 
      Enddo
      Enddo 

      Do e = 1, Nvrt_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do m = 1, Nocc_a
      Do n = 1, Nocc_b
      Do f = 1, Nvrt_a
         T = T + T2ab(F,a,M,n)*W_ab(M,n,F,e)
      Enddo
      Enddo
      Enddo
         Qae_bb(a,e)  = Qae_bb(a,e) - T
      Enddo
      Enddo 

#ifdef _DEBUG_LVL0
      call checksum("Qae_aa    :",Qae_aa,Nvrt_a*Nvrt_a)
      call checksum("Qae_bb    :",Qae_bb,Nvrt_b*Nvrt_b)
#endif  

      Return 
      End 
