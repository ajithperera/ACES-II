      Subroutine Form_qef1(Qef_aa,Qef_bb,Hbar_aibc1,Hbar_aibc2,
     +                     Hbar_aibc3,Hbar_aibc4,T1aa,T1bb,Nocc_a,
     +                     Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)
 
      Dimension Qef_aa(Nvrt_a,Nvrt_a)
      Dimension Qef_bb(Nvrt_b,Nvrt_b)

      Dimension Hbar_aibc1(Nvrt_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension Hbar_aibc2(Nvrt_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension Hbar_aibc3(Nvrt_a,Nocc_b,Nvrt_a,Nvrt_b)
      Dimension Hbar_aibc4(Nvrt_b,Nocc_a,Nvrt_b,Nvrt_a)

      Integer A,F,M,E

#include "ccsdlight_vars.com"

C AA block 
C Qef(E,F)= +T1(A,M)*W(E,M,F,A)
      
      Do E = 1, Nvrt_a
      Do F = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
      Do A = 1, Nvrt_a
         T = T + T1aa(A,M)*Hbar_aibc1(E,M,F,A)
      Enddo
      Enddo
      Qef_aa(E,F) =  T
      Enddo 
      Enddo 

C Qij(E,F)= -T1(a,m)*W(E,m,F,a)

      Do E = 1, Nvrt_a
      Do F = 1, Nvrt_a
         T = 0.0D0
      Do m = 1, Nocc_b
      Do a = 1, Nvrt_b
         T = T + T1bb(a,m)*Hbar_aibc3(E,m,F,a)
      Enddo
      Enddo
      Qef_aa(E,F) = Qef_aa(E,F) + T
      Enddo
      Enddo

C BB block
C Qef(e,f)= +T1(a,m)*W(e,m,f,a)

      Do E = 1, Nvrt_b
      Do F = 1, Nvrt_b
         T = 0.0D0
      Do m = 1, Nocc_b
      Do a = 1, Nvrt_b
         T = T + T1bb(a,m)*Hbar_aibc2(e,m,f,a)
      Enddo
      Enddo
      Qef_bb(e,f) =  T
      Enddo
      Enddo

C Qef(e,f)= +T1(A,M)*W(e,M,f,A)

      Do E = 1, Nvrt_b
      Do F = 1, Nvrt_b
         T = 0.0D0
      Do M = 1, Nocc_a
      Do A = 1, Nvrt_a
         T = T + T1aa(A,M)*Hbar_aibc4(e,M,f,A)
      Enddo
      Enddo
      Qef_bb(e,f) = Qef_bb(e,f) + T
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      L_aa = Nvrt_a*Nvrt_a
      L_bb = Nvrt_b*Nvrt_b
      call checksum("Qef_aa    :",Qef_aa,L_aa)
      call checksum("Qef_bb    :",Qef_bb,L_bb)
#endif  

      Return 
      End 
