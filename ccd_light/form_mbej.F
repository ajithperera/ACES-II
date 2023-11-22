      Subroutine Form_mbej(T2aa,T2bb,T2ab,W3_aaaa,W3_bbbb,W3_abab,
     +                     W3_baba,W3_abba,W3_baab,W_aa,W_bb,W_ab,
     +                     Wmbej_1,Wmbej_2,Wmbej_3,Wmbej_4,Wmbej_5,
     +                     Wmbej_6,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_b)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension W3_aaaa(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W3_bbbb(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W3_abab(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W3_baba(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension W3_abba(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension W3_baab(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Dimension Wmbej_1(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Wmbej_2(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Wmbej_3(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Wmbej_4(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension Wmbej_5(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension Wmbej_6(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Integer M,N,A,F,E,B,J

#include "ccdlight_vars.com"

#ifdef _DEBUG_LVLM
      Write(6,*)
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("T2aa      :",T2aa,L_aaaa)
      call checksum("T2bb      :",T2bb,L_bbbb)
      call checksum("T2ab      :",T2ab,L_abab)
      L_aaaa = Nocc_a*Nvrt_a*Nvrt_a*Nocc_a
      L_bbbb = Nocc_b*Nvrt_b*Nvrt_b*Nocc_b
      L_abab = Nocc_a*Nvrt_b*Nvrt_a*Nocc_b
      L_baba = Nocc_b*Nvrt_a*Nvrt_b*Nocc_a
      L_abba = Nocc_a*Nvrt_b*Nocc_a*Nvrt_b
      L_baab = Nocc_b*Nvrt_a*Nocc_b*Nvrt_a
      call checksum("Form_wMBEJ:",W3_aaaa,L_aaaa)
      call checksum("Form_wmbej:",W3_bbbb,L_bbbb)
      call checksum("Form_wMbEj:",W3_abab,L_abab)
      call checksum("Form_wmBeJ:",W3_baba,L_baba)
      call checksum("Form_wmBEj:",W3_baab,L_baab)
      call checksum("Form_wMbeJ:",W3_abba,L_abba)
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("Form_wMNEF:",W_aa,L_aaaa)
      call checksum("Form_wmnef:",W_bb,L_bbbb)
      call checksum("Form_wMnEf:",W_ab,L_abab)
#endif

C Wmbej_1(MB,EJ) = W2(MB,EJ) - 1/2T2(FB,JN)*W2(MN,EF) - 1/2T2(fB,Jn)*W2(Mn,Ef)
 
      Do J = 1, Nocc_a
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do M = 1, Nocc_a
         Wmbej_1(M,B,E,J) = W3_aaaa(M,B,E,J)
      Do F = 1, Nvrt_a
      Do N = 1, Nocc_a
      C = Ocn_oa(N)*(1.0D0-Ocn_va(F))
      Wmbej_1(M,B,E,J)= Wmbej_1(M,B,E,J) -
     +                  0.50D0*T2aa(F,B,J,N)*W_aa(M,N,E,F)*C
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
      C = Ocn_ob(n)*(1.0D0-Ocn_vb(f))
      Wmbej_1(M,B,E,J)= Wmbej_1(M,B,E,J) + 
     +                  0.50D0*T2ab(B,f,J,n)*W_ab(M,n,E,f)*C
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C Wmbej_1(mb,ej) = W2(mb,ej) - 1/2T2(fb,jn)*W2(mn,ef) + 1/2T2(Fb,jN)*W2(mN,eF)

      Do J = 1, Nocc_b
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do M = 1, Nocc_b
         Wmbej_2(m,b,e,j) = W3_bbbb(m,b,e,j)
      Do F = 1, Nvrt_b
      Do N = 1, Nocc_b
      C = Ocn_ob(n)*(1.0D0-Ocn_vb(f))
      Wmbej_2(m,b,e,j)= Wmbej_2(m,b,e,j)  - 
     +                  0.50D0*T2bb(f,b,j,n)*W_bb(m,n,e,f)*C
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
      Do N = 1, Nocc_a
      C = Ocn_oa(N)*(1.0D0-Ocn_va(F))
      Wmbej_2(m,b,e,j)= Wmbej_2(m,b,e,j) + 
     +                  0.50D0*T2ab(F,b,N,j)*W_ab(N,m,F,e)*C 
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C Wmbej_3(Mb,Ej) = W(Mb,Ej) - 1/2T2(fb,jn)*W(Mn,Ef) + 1/2T2(Fb,jN)*W(MN,EF)

      Do J = 1, Nocc_b 
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do M = 1, Nocc_a
         Wmbej_3(M,b,E,j) = W3_abab(M,b,E,j) 
      Do F = 1, Nvrt_b
      Do N = 1, Nocc_b
      C = Ocn_ob(n)*(1.0D0-Ocn_vb(f))
      Wmbej_3(M,b,E,j)= Wmbej_3(M,b,E,j)  + 
     +                  0.50D0*T2bb(b,f,j,n)*W_ab(M,n,E,f)*C
 
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
      C = Ocn_oa(N)*(1.0D0-Ocn_va(F))
      Wmbej_3(M,b,E,j)= Wmbej_3(M,b,E,j) +  
     +                  0.50D0*T2ab(F,b,N,j)*W_aa(M,N,E,F)*C
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C Wmbej_4(mB,eJ) = W2(mB,eJ) - 1/2T2(BF,JN)*W2(mN,eF) + 1/2T2(Bf,Jn)*W2(mn,ef)

      Do J = 1, Nocc_a
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_a
      Do M = 1, Nocc_b
         Wmbej_4(m,B,e,J) = W3_baba(m,B,e,J) 
      Do F = 1, Nvrt_a
      Do N = 1, Nocc_a
      C = Ocn_oa(N)*(1.0D0-Ocn_va(F))
      Wmbej_4(m,B,e,J)= Wmbej_4(m,B,e,J) + 
     +                  0.50D0*T2aa(B,F,J,N)*W_ab(N,m,F,e)*C
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
      C = Ocn_ob(n)*(1.0D0-Ocn_vb(f))
      Wmbej_4(m,B,e,J)= Wmbej_4(m,B,e,J) + 
     +                  0.50D0*T2ab(B,f,J,n)*W_bb(m,n,e,f)*C
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C Wmbej_5(Mb,eJ) = W2(Mb,eJ) + 1/2T2(bF,JN)*W2(Mn,Fe) 

      Do J = 1, Nocc_a
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do M = 1, Nocc_a
         Wmbej_5(M,b,J,e) = W3_abba(M,b,J,e)
      Do F = 1, Nvrt_a
      Do N = 1, Nocc_b
      C = Ocn_ob(n)*(1.0D0-Ocn_va(F))
      Wmbej_5(M,b,J,e)=  Wmbej_5(M,b,J,e) -
     +                  0.50D0*T2ab(F,b,J,n)*W_ab(M,n,F,e)*C
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C Wmbej_6(mB,Ej) = -W2(mB,jE) + 1/2T2(fB,jN)*W2(mN,fE) 

      Do J = 1, Nocc_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a 
      Do M = 1, Nocc_b
         Wmbej_6(m,B,j,E) = W3_baab(m,B,j,E)
      Do F = 1, Nvrt_b
      Do N = 1, Nocc_a
      C = Ocn_oa(N)*(1.0D0-Ocn_vb(f))
      Wmbej_6(m,B,j,E)= Wmbej_6(m,B,j,E) -
     +                  0.50D0*T2ab(B,f,N,j)*W_ab(N,m,E,f)*C
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      L_aaaa = Nocc_a*Nvrt_a*Nvrt_a*Nocc_a
      L_bbbb = Nocc_b*Nvrt_b*Nvrt_b*Nocc_b
      L_abab = Nocc_a*Nvrt_b*Nvrt_a*Nocc_b
      L_baba = Nocc_b*Nvrt_a*Nvrt_b*Nocc_a
      L_abba = Nocc_a*Nvrt_b*Nocc_a*Nvrt_b
      L_baab = Nocc_b*Nvrt_a*Nocc_b*Nvrt_a
      call checksum("Form_wMBEJ:",Wmbej_1,L_aaaa)
      call checksum("Form_wmbej:",Wmbej_2,L_bbbb)
      call checksum("Form_wMbEj:",Wmbej_3,L_abab)
      call checksum("Form_wmBeJ:",Wmbej_4,L_baba)
      call checksum("Form_wmBEj:",Wmbej_5,L_abba)
      call checksum("Form_wMbeJ:",Wmbej_6,L_baab)
#endif

      Return
      End
