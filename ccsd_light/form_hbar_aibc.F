      Subroutine Form_hbar_aibc(Hbar_aibc1,Hbar_aibc2,Hbar_aibc3,
     +                          Hbar_aibc4,T1aa,T1bb,W_aa,W_bb,
     +                          W_ab,W4_aa,W4_bb,W4_ab,W4_ba,Nocc_a,
     +                          Nocc_b,Nvrt_a,Nvrt_b,Scale)
     
      Implicit Double Precision(A-H,O-Z)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_b)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension W4_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W4_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W4_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W4_ba(Nvrt_a,Nvrt_b,Nocc_a,Nvrt_b)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension Hbar_aibc1(Nvrt_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension Hbar_aibc2(Nvrt_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension Hbar_aibc3(Nvrt_a,Nocc_b,Nvrt_a,Nvrt_b)
      Dimension Hbar_aibc4(Nvrt_b,Nocc_a,Nvrt_b,Nvrt_a)

#include "ccsdlight_vars.com"

      Integer A,B,C,I,M
      Data Dzero,Done /0.0D0,1.0D0/

C Hbar_aibc(A,I,B,C) = W(A,I,B,C) - T1(A,M)*<MI||BC>

      Do C = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do I = 1, Nocc_a 
      Do A = 1, Nvrt_a 
         If (Scale .Gt. Dzero) Then
            D = (1.0D0-Ocn_va(A))
         Else
            D = Done 
         Endif
         Hbar_aibc1(A,I,B,C) = W4_aa(B,C,A,I)*D
      Do M = 1, Nocc_a
         Hbar_aibc1(A,I,B,C) =  Hbar_aibc1(A,I,B,C) -
     +                          T1aa(A,M)*W_aa(M,I,B,C)*Scale
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C Hbar_aibc(a,i,b,c) = W(a,i,b,c) - T1(a,m)*<mi||bc>

      Do C = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do I = 1, Nocc_b
      Do A = 1, Nvrt_b
         If (Scale .Gt. Dzero) Then
            D = (1.0D0-Ocn_vb(a))
         Else
            D = Done 
         Endif
         Hbar_aibc2(a,i,b,c) = W4_bb(b,c,a,i)*D
      Do M = 1, Nocc_b
         Hbar_aibc2(a,i,b,c) =  Hbar_aibc2(a,i,b,c) -
     +                          T1bb(a,m)*W_bb(m,i,b,c)*Scale 
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C Hbar_aibc(A,i,B,c) = W(A,i,B,c) - T1(A,M)*<Mi||Bc>

      Do C = 1, Nvrt_b
      Do B = 1, Nvrt_a
      Do I = 1, Nocc_b
      Do A = 1, Nvrt_a
         If (Scale .Gt. Dzero) Then
            D = (1.0D0-Ocn_va(A))
         Else
            D = Done 
         Endif
         Hbar_aibc3(A,i,B,c) = W4_ab(A,c,B,i)*D
      Do M = 1, Nocc_a
         Hbar_aibc3(A,i,B,c) =  Hbar_aibc3(A,i,B,c) -
     +                          T1aa(A,M)*W_ab(M,i,B,c)*scale 
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C Hbar_aibc(a,I,b,C) = W(a,I,b,C) - T1(a,m)*<mI||bC>

      Do C = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do I = 1, Nocc_a
      Do A = 1, Nvrt_b
         If (Scale .Gt. Dzero) Then
            D = (1.0D0-Ocn_vb(a))
         Else
            D = Done 
         Endif
         Hbar_aibc4(a,I,b,C) = W4_ba(C,b,I,a)*D
      Do M = 1, Nocc_b
         Hbar_aibc4(a,I,b,C) = Hbar_aibc4(a,I,b,C) -
     +                         T1bb(a,m)*W_ab(I,m,C,b)*scale 
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      If (scale .eq. 1.0D0) Then
      L_aaaa =Nvrt_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb =Nvrt_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab =Nvrt_a*Nocc_b*Nvrt_a*Nvrt_b
      L_baba =Nvrt_b*Nocc_a*Nvrt_b*Nvrt_a
      call checksum("Hbar_aibc1:",Hbar_aibc1,L_aaaa)
      call checksum("Hbar_aibc2:",Hbar_aibc2,L_bbbb)
      call checksum("Hbar_aibc3:",Hbar_aibc3,L_abab)
      call checksum("Hbar_aibc4:",Hbar_aibc4,L_baba)
      Endif 
#endif

      Return
      End
