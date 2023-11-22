      Subroutine Form_qij1(Qij_aa,Qij_bb,Hbar_ijka1,Hbar_ijka2,
     +                     Hbar_ijka3,Hbar_ijka4,T1aa,T1bb,Nocc_a,
     +                     Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension Qij_aa(Nocc_a,Nocc_a)
      Dimension Qij_bb(Nocc_b,Nocc_b)

      Dimension Hbar_ijka1(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension Hbar_ijka2(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension Hbar_ijka3(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension Hbar_ijka4(Nocc_b,Nocc_a,Nocc_b,Nvrt_a)
       
      Integer I,J,M,E

#include "ccsdlight_vars.com"

C AA block 
C Qij(I,J)= T1(E,M)*W(I,J,M,E)
      
      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
         T = 0.0D0
      Do m = 1, Nocc_a
      Do e = 1, Nvrt_a
         T = T + T1aa(E,M)*Hbar_ijka1(I,M,J,E)
      Enddo
      Enddo
      Qij_aa(I,J) = T
      Enddo 
      Enddo 

C Qij(I,J)= T1(e,m)*W(I,m,J,e)

      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
         T = 0.0D0
      Do m = 1, Nocc_b
      Do e = 1, Nvrt_b
         T = T + T1bb(e,m)*Hbar_ijka3(I,m,J,e)
      Enddo
      Enddo
      Qij_aa(I,J) = Qij_aa(I,J) + T
      Enddo
      Enddo

C BB block
C Qij(i,i)= T1(e,m)*W(i,m,j,e)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
         T = 0.0D0
      Do m = 1, Nocc_b
      Do e = 1, Nvrt_b
         T = T + T1bb(e,m)*Hbar_ijka2(i,m,j,e)
      Enddo
      Enddo
      Qij_bb(i,j) = T
      Enddo
      Enddo

C Qij(i,i)= T1(E,M)*W(i,M,j,E)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
         T = 0.0D0
      Do m = 1, Nocc_a
      Do e = 1, Nvrt_a
         T = T + T1aa(E,M)*Hbar_ijka4(i,M,j,E)
      Enddo
      Enddo
      Qij_bb(i,j) = Qij_bb(i,j) + T
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      L_aa = Nocc_a*Nocc_a
      L_bb = Nocc_b*Nocc_b
      call checksum("Qij_aa    :",Qij_aa,L_aa)
      call checksum("Qij_bb    :",Qij_bb,L_bb)
#endif  

      Return 
      End 
