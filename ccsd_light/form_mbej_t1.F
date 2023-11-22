      Subroutine Form_mbej_t1(T1aa,T1bb,W_aa,W_bb,W_ab,W4_aa,W4_bb,
     +                        W4_ab,W4_ba,W5_aa,W5_bb,W5_ab,W5_ba,
     +                        Wmbej_1,Wmbej_2,Wmbej_3,Wmbej_4,Wmbej_5,
     +                        Wmbej_6,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_b)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension W4_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W4_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W4_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W4_ba(Nvrt_a,Nvrt_b,Nocc_a,Nvrt_b)

      Dimension W5_aa(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension W5_bb(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension W5_ab(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension W5_ba(Nocc_a,Nocc_b,Nvrt_a,Nocc_b)

      Dimension Wmbej_1(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Wmbej_2(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Wmbej_3(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Wmbej_4(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension Wmbej_5(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension Wmbej_6(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Integer M,N,A,F,E,B,J

#include "ccsdlight_vars.com"

#ifdef _DEBUG_LVLM
      Write(6,*)
      L_aa = Nocc_a*Nvrt_a
      L_bb = Nocc_b*Nvrt_b
      call checksum("T1aa      :",T1aa,L_aa)
      call checksum("T1bb      :",T1bb,L_bb)
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
      call checksum("Form_wmBEj:",Wmbej_5,L_baab)
      call checksum("Form_wMbeJ:",Wmbej_6,L_abba)
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("Form_wMNEF:",W_aa,L_aaaa)
      call checksum("Form_wmnef:",W_bb,L_bbbb)
      call checksum("Form_wMnEf:",W_ab,L_abab)
#endif

C Wmbej_1(MB,EJ) = -T1(F,J)*W(E,F,B,M) + T1(B,N)*W(M,N,J,E) -
C                   T1(F,J)*T1(B,N)*W(M,N,E,F) 
 
      Do J = 1, Nocc_a
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do M = 1, Nocc_a 
         T = 0.0D0 
      Do F = 1, Nvrt_a
CSSS         C = (1.0D0-Ocn_va(F))
         C = (1.0D0-Ocn_va(B))
         T = T + T1aa(F,J)*W4_aa(E,F,B,M)*C
      Enddo
         Wmbej_1(M,B,E,J) = Wmbej_1(M,B,E,J) - T
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_a
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do M = 1, Nocc_a
         T = 0.0D0
      Do N = 1, Nocc_a
CSSS         C = Ocn_oa(N)
         C = Ocn_oa(J)
         T = T + T1aa(B,N)*W5_aa(M,N,J,E)*C
      Enddo
         Wmbej_1(M,B,E,J) = Wmbej_1(M,B,E,J) + T
      Enddo
      Enddo
      Enddo
      Enddo
 
      Do J = 1, Nocc_a
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do M = 1, Nocc_a
         T = 0.0D0
      Do F = 1, Nvrt_a
      Do N = 1, Nocc_a
CSSS         C = Ocn_oa(N)*(1.0D0-Ocn_va(F))
         C = 1.0D0
         T = T + T1aa(F,J)*T1aa(B,N)*W_aa(M,N,E,F)*C
      Enddo
      Enddo
         Wmbej_1(M,B,E,J) = Wmbej_1(M,B,E,J) - T
      Enddo
      Enddo
      Enddo
      Enddo

C Wmbej_2(mb,ej) = -T1(f,j)*W(e,f,b,m) + T1(b,n)*W(m,n,j,e) -
C                   T1(f,j)*T1(b,n)*W(w,n,e,f)

      Do J = 1, Nocc_b
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do M = 1, Nocc_b
         T = 0.0D0 
      Do F = 1, Nvrt_b
CSSS        C = (1.0D0-Ocn_vb(f))
        C = (1.0D0-Ocn_vb(b))
        T = T + T1bb(f,j)*W4_bb(e,f,b,m)*C
      Enddo
        Wmbej_2(m,b,e,j) = Wmbej_2(m,b,e,j) - T
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_b
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do M = 1, Nocc_b
         T = 0.0D0
      Do N = 1, Nocc_b
CSSS        C = Ocn_ob(n)
        C = Ocn_ob(j)
        T = T + T1bb(b,n)*W5_bb(m,n,j,e)*C
      Enddo
        Wmbej_2(m,b,e,j) = Wmbej_2(m,b,e,j) + T
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_b
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do M = 1, Nocc_b
         T = 0.0D0
      Do F = 1, Nvrt_b
      Do N = 1, Nocc_b
CSSS         C = Ocn_ob(n)*(1.0D0-Ocn_vb(f))
         C = 1.0D0
         T = T + T1bb(f,j)*T1bb(b,n)*W_bb(m,n,e,f)*C
      Enddo
      Enddo
         Wmbej_2(m,b,e,j) = Wmbej_2(m,b,e,j) - T
      Enddo
      Enddo
      Enddo
      Enddo

C Wmbej_3(Mb,Ej) = T1(f,j)*W(E,f,M,b) - T1(b,n)*W(M,n,E,j) -
C                  T1(f,j)*T1(b,n)W(M,n,E,f)

      Do J = 1, Nocc_b 
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do M = 1, Nocc_a
         T = 0.0D0
      Do F = 1, Nvrt_b
CSSS         C = (1.0D0-Ocn_vb(f))
         C = (1.0D0-Ocn_vb(b))
         T = T + T1bb(f,j)*W4_ba(E,f,M,b)*C
      Enddo
         Wmbej_3(M,b,E,j) = Wmbej_3(M,b,E,j) + T
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do M = 1, Nocc_a
         T = 0.0D0
      Do N = 1, Nocc_b
CSSS         C = Ocn_ob(n)
         C = Ocn_ob(j)
         T = T + T1bb(b,n)*W5_ba(M,n,E,j)*C
      Enddo
         Wmbej_3(M,b,E,j) = Wmbej_3(M,b,E,j) - T
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do M = 1, Nocc_a
         T = 0.0D0
      Do N = 1, Nocc_b
      Do F = 1, Nvrt_b 
CSSS         C = Ocn_ob(n)*(1.0D0-Ocn_vb(f))
         C = 1.0D0
         T = T + T1bb(f,j)*T1bb(b,n)*W_ab(M,n,E,f)*C
      Enddo
      Enddo 
         Wmbej_3(M,b,E,j) = Wmbej_3(M,b,E,j) - T
      Enddo
      Enddo
      Enddo
      Enddo

C Wmbej_4(mB,eJ) = T1(J,F)*W(F,e,B,m) - T1(B,N)*W(N,m,J,e) -
C                  T1(J,F)*T1(B,N)W(N,m,F,e)

      Do J = 1, Nocc_a
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_a
      Do M = 1, Nocc_b
         T = 0.0D0
      Do F = 1, Nvrt_a
CSSS         C =(1.0D0-Ocn_va(F))
         C = (1.0D0-Ocn_va(B))
         T = T + T1aa(F,J)*w4_ab(F,e,B,m)*C
      Enddo
         Wmbej_4(m,B,e,J) = Wmbej_4(m,B,e,J) + T
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_a
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_a
      Do M = 1, Nocc_b
         T = 0.0D0
      Do N = 1, Nocc_a
CSSS         C = Ocn_oa(N)
         C = Ocn_oa(J)
         T = T + T1aa(B,N)*W5_ab(N,m,J,e)*C
      Enddo
         Wmbej_4(m,B,e,J) = Wmbej_4(m,B,e,J) - T
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_a
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_a
      Do M = 1, Nocc_b
         T = 0.0D0
      Do N = 1, Nocc_a
      Do F = 1, Nvrt_a 
CSSS         C = Ocn_oa(N)*(1.0D0-Ocn_va(F))
         C = 1.0D0
         T = T + T1aa(B,N)*T1aa(F,J)*W_ab(N,m,F,e)*c
      Enddo
      Enddo
         Wmbej_4(m,B,e,J) = Wmbej_4(m,B,e,J) - T
      Enddo
      Enddo
      Enddo
      Enddo

C Wmbej_5(Mb,Je) = -T1(F,J)*W(F,e,M,b) + T(b,n)*W(M,n,J,e) +
C                   T1(F,J)*T(b,n)*W(M,n,F,e)

      Do J = 1, Nocc_a
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do M = 1, Nocc_a
         T = 0.0D0
      Do F = 1, Nvrt_a
CSSS         C = (1.0D0-Ocn_va(F))
         C = (1.0D0-Ocn_vb(b))
         T = T + T1aa(F,J)*W4_ba(F,e,M,b)*C
      Enddo
         Wmbej_5(M,b,J,e) = Wmbej_5(M,b,J,e) + T
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_a
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do M = 1, Nocc_a
         T = 0.0D0
      Do N = 1, Nocc_b
CSSS         C = Ocn_ob(n)
         C = Ocn_oa(J)
         T = T + T1bb(b,n)*W5_ab(M,n,J,e)*C
      Enddo
         Wmbej_5(M,b,J,e) = Wmbej_5(M,b,J,e) - T
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_a
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do M = 1, Nocc_a
         T = 0.0D0
      Do N = 1, Nocc_b
      Do F = 1, Nvrt_a
CSSS         C = Ocn_ob(n)*(1.0D0-Ocn_va(F))
         C = 1.0D0
         T = T + T1aa(F,J)*T1bb(b,n)*W_ab(M,n,F,e)*C
      Enddo
      Enddo
         Wmbej_5(M,b,J,e) = Wmbej_5(M,b,J,e) - T
      Enddo
      Enddo
      Enddo
      Enddo

C Wmbej_6(mB,Ej) = -T1(f,j)*W(E,f,B,m) - T1(B,N)*W(N,m,J,e) + 
C                   T1(f,j)*T1(B,N)*W(N,m,E,f)

      Do J = 1, Nocc_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a 
      Do M = 1, Nocc_b
         T = 0.0D0
      Do F = 1, Nvrt_b
CSSS         C = (1.0D0-Ocn_vb(f))
         C = (1.0D0-Ocn_va(B))
         T = T + T1bb(f,j)*W4_ab(E,f,B,m)*C
      Enddo
         Wmbej_6(m,B,j,E) = Wmbej_6(m,B,j,E) + T
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do M = 1, Nocc_b
         T = 0.0D0
      Do N = 1, Nocc_a
CSSS         C = Ocn_oa(N)
         C = Ocn_ob(j)
         T = T + T1aa(B,N)*W5_ba(N,m,E,j)*C
      Enddo
         Wmbej_6(m,B,j,E) = Wmbej_6(m,B,j,E) - T
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do M = 1, Nocc_b
         T = 0.0D0
      Do N = 1, Nocc_a
      Do F = 1, Nvrt_b
CSSS         C = Ocn_oa(N)*(1.0D0-Ocn_vb(f))
         C = 1.0D0
         T = T + T1aa(B,N)*T1bb(f,j)*W_ab(N,m,E,f)*C 
      Enddo
      Enddo
         Wmbej_6(m,B,j,E) = Wmbej_6(m,B,j,E) - T
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
