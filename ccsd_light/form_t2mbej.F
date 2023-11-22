      Subroutine Form_t2mbej(T2aa,T2bb,T2ab,W_aa,W_bb,W_ab,T2mbej_1,
     +                       T2mbej_2,T2mbej_3,T2mbej_4,T2mbej_5,
     +                       T2mbej_6,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_b)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension T2mbej_1(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension T2mbej_2(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension T2mbej_3(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension T2mbej_4(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension T2mbej_5(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension T2mbej_6(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Data Dzero /0.0D0/

      Integer M,N,A,F,E,B,J

#include "ccsdlight_vars.com"

C T2mbej_1(MB,EJ) = - 1/2T2(FB,JN)*W2(MN,EF) - 1/2T2(fB,Jn)*W2(Mn,Ef)
 
      Do J = 1, Nocc_a
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do M = 1, Nocc_a
         T2mbej_1(M,B,E,J) = Dzero 
      Do F = 1, Nvrt_a
      Do N = 1, Nocc_a
      T2mbej_1(M,B,E,J)= T2mbej_1(M,B,E,J) -
     +                   0.50D0*T2aa(F,B,J,N)*W_aa(M,N,E,F)
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_a
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do M = 1, Nocc_a
      Do F = 1, Nvrt_b
      Do N = 1, Nocc_b
      T2mbej_1(M,B,E,J)= T2mbej_1(M,B,E,J) + 
     +                   0.50D0*T2ab(B,f,J,n)*W_ab(M,n,E,f)
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C T2mbej_2(mb,ej) = - 1/2T2(fb,jn)*W2(mn,ef) + 1/2T2(Fb,jN)*W2(mN,eF)

      Do J = 1, Nocc_b
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do M = 1, Nocc_b
         T2mbej_2(m,b,e,j) = Dzero 
      Do F = 1, Nvrt_b
      Do N = 1, Nocc_b
      T2mbej_2(m,b,e,j)= T2mbej_2(m,b,e,j)  - 
     +                   0.50D0*T2bb(f,b,j,n)*W_bb(m,n,e,f)
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_b
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do M = 1, Nocc_b
      Do F = 1, Nvrt_a
      do N = 1, Nocc_a
      T2mbej_2(m,b,e,j)= T2mbej_2(m,b,e,j) + 
     +                   0.50D0*T2ab(F,b,N,j)*W_ab(N,m,F,e)
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C T2mbej_3(Mb,Ej) = - 1/2T2(fb,jn)*W(Mn,Ef) + 1/2T2(Fb,jN)*W(MN,EF)

      Do J = 1, Nocc_b 
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do M = 1, Nocc_a
         T2mbej_3(M,b,E,j) = Dzero 
      Do F = 1, Nvrt_b
      Do N = 1, Nocc_b
      T2mbej_3(M,b,E,j)= T2mbej_3(M,b,E,j)  + 
     +                   0.50D0*T2bb(b,f,j,n)*W_ab(M,n,E,f)
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
  
      Do J = 1, Nocc_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do M = 1, Nocc_a
      Do F = 1, Nvrt_a
      Do N = 1, Nocc_a
      T2mbej_3(M,b,E,j)= T2mbej_3(M,b,E,j) +  
     +                   0.50D0*T2ab(F,b,N,j)*W_aa(M,N,E,F)
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C T2mbej_4(mB,eJ) = - 1/2T2(BF,JN)*W2(mN,eF) + 1/2T2(Bf,Jn)*W2(mn,ef)

      Do J = 1, Nocc_a
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_a
      Do M = 1, Nocc_b
         T2mbej_4(m,B,e,J) = Dzero
      Do F = 1, Nvrt_a
      Do N = 1, Nocc_a
      T2mbej_4(m,B,e,J)= T2mbej_4(m,B,e,J) + 
     +                   0.50D0*T2aa(B,F,J,N)*W_ab(N,m,F,e)
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_a
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_a
      Do M = 1, Nocc_b
      Do F = 1, Nvrt_b
      Do N = 1, Nocc_b
      T2mbej_4(m,B,e,J)= T2mbej_4(m,B,e,J) + 
     +                   0.50D0*T2ab(B,f,J,n)*W_bb(m,n,e,f)
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C T2mbej_5(Mb,Je) = + 1/2T2(bF,JN)*W2(Mn,Fe) 

      Do J = 1, Nocc_a
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do M = 1, Nocc_a
         T2mbej_5(M,b,J,e) = Dzero 
      Do F = 1, Nvrt_a
      Do N = 1, Nocc_b
      T2mbej_5(M,b,J,e)=  T2mbej_5(M,b,J,e) -
     +                    0.50D0*T2ab(F,b,J,n)*W_ab(M,n,F,e)
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C Wmbej_6(mB,jE) = -W2(mB,jE) + 1/2T2(fB,jN)*W2(mN,fE) 

      Do J = 1, Nocc_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a 
      Do M = 1, Nocc_b
         T2mbej_6(m,B,j,E) = DZero
      Do F = 1, Nvrt_b
      Do N = 1, Nocc_a
      T2mbej_6(m,B,j,E)= T2mbej_6(m,B,j,E) -
     +                   0.50D0*T2ab(B,f,N,j)*W_ab(N,m,E,f)
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVLM
      L_aaaa = Nocc_a*Nvrt_a*Nvrt_a*Nocc_a
      L_bbbb = Nocc_b*Nvrt_b*Nvrt_b*Nocc_b
      L_abab = Nocc_a*Nvrt_b*Nvrt_a*Nocc_b
      L_baba = Nocc_b*Nvrt_a*Nvrt_b*Nocc_a
      L_abba = Nocc_a*Nvrt_b*Nocc_a*Nvrt_b
      L_baab = Nocc_b*Nvrt_a*Nocc_b*Nvrt_a
      call checksum("Form_wMBEJ:",T2mbej_1,L_aaaa)
      call checksum("Form_wmbej:",T2mbej_2,L_bbbb)
      call checksum("Form_wMbEj:",T2mbej_3,L_abab)
      call checksum("Form_wmBeJ:",T2mbej_4,L_baba)
      call checksum("Form_wmBEj:",T2mbej_5,L_abba)
      call checksum("Form_wMbeJ:",T2mbej_6,L_baab)
#endif

      Return
      End
