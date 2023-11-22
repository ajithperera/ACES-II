      Subroutine Form_htmi(T1taup_aa,T1taup_bb,T1taup_ab,T1aa,T1bb,
     +                     Fmi_a,Fmi_b,W_aa,W_bb,W_ab,W5_aa,W5_bb,
     +                     W5_ab,W5_ba,Fockoo_od_a,Fockoo_od_b,
     +                     Fockov_a,Fockov_b,Nocc_a,Nocc_b,Nvrt_a,
     +                     Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension T1taup_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T1taup_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T1taup_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_b)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension W5_aa(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension W5_bb(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension W5_ab(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension W5_ba(Nocc_a,Nocc_b,Nvrt_a,Nocc_b)

      Dimension Fockoo_od_a(Nocc_a,Nocc_a)
      Dimension Fockoo_od_b(Nocc_b,Nocc_b)
      Dimension Fockov_a(Nocc_a,Nvrt_a)
      Dimension Fockov_b(Nocc_b,Nvrt_b)

      Dimension Fmi_a(Nocc_a,Nocc_a)
      Dimension Fmi_b(Nocc_b,Nocc_b)

      Integer M,N,E,F,I

#include "cse12_vars.com"

C  Fmi_a(M,I) <- T1(E,N)*W(M,N,I,E) + T1(e,n)*W(M,n,I,e) + (1-delta_mi)*fockoo(M,I)

      Do i = 1, Nocc_a
      Do m = 1, Nocc_a
         T = 0.0D0
         C = Ocn_oa(I)
         Fmi_a(M,I) = Fockoo_od_a(M,I)*C
      Do e = 1, Nvrt_a
      Do N = 1, Nocc_a
         C = Ocn_oa(I)
         T = T + T1aa(E,N)*W5_aa(M,N,I,E)*C 
      Enddo
      Enddo
      Fmi_a(M,I) = Fmi_a(M,I) + T
      Enddo
      Enddo

      Do i = 1, Nocc_a
      Do m = 1, Nocc_a
         T = 0.0D0
      Do e = 1, Nvrt_b
      Do N = 1, Nocc_b
         C = Ocn_oa(I)
         T = T + T1bb(e,n)*W5_ab(M,n,I,e)*C 
      Enddo
      Enddo
      Fmi_a(M,I) = Fmi_a(M,I) + T
      Enddo
      Enddo

C  Fmi_b(m,i) = T1(e,n)*W(m,n,i,e) + T1(E,N)*W(N,m,E,i) + (1-delta_mi)*Fockoo(m,i) 

      Do i = 1, Nocc_b
      Do m = 1, Nocc_b
         T = 0.0D0
         C = Ocn_ob(i)
         Fmi_b(m,i) = Fockoo_od_b(m,i)*C
      Do e = 1, Nvrt_b
      Do N = 1, Nocc_b
         C = Ocn_ob(i)
         T = T + T1bb(e,n)*W5_bb(m,n,i,e)*C
      Enddo
      Enddo
      Fmi_b(m,i) = Fmi_b(m,i) + T
      Enddo
      Enddo

      Do i = 1, Nocc_b
      Do m = 1, Nocc_b
         T = 0.0D0
      Do e = 1, Nvrt_a
      Do N = 1, Nocc_a
         C = Ocn_ob(i)
         T = T + T1aa(E,N)*W5_ba(N,m,E,i)*C
      Enddo
      Enddo
      Fmi_b(m,i) = Fmi_b(m,i) + T
      Enddo
      Enddo

C  Fmi_a(M,I) =  +1/2f(M,E)*T1(E,I)

      Do i = 1, Nocc_a
      Do m = 1, Nocc_a
         T = 0.0D0
      Do e = 1, Nvrt_a
         T = T + 0.50D0*T1aa(E,I)*Fockov_a(M,E)
      Enddo
      Fmi_a(m,i) = Fmi_a(m,i) + T
      Enddo
      Enddo

C  Fmi_b(m,i) =  +1/2f(m,e)*T1(e,i)

      Do i = 1, Nocc_b
      Do m = 1, Nocc_b
         T = 0.0D0
      Do e = 1, Nvrt_b
         T = T + 0.50D0*T1bb(e,i)*Fockov_b(m,e)
      Enddo
      Fmi_b(m,i) = Fmi_b(m,i) + T
      Enddo
      Enddo

C Fmi_a(M,I)<- 1/2 T1taup(EF,IN)*W2(MN,EF) + T1taup(Ef,In)*W(Mn,Ef)

      Do I = 1, Nocc_a
      Do M = 1, Nocc_a
         T = 0.0D0
      Do E = 1, Nvrt_a
      Do F = 1, Nvrt_a
      Do N = 1, Nocc_a
         T = T + 0.50D0*T1taup_aa(E,F,I,N)*W_aa(M,N,E,F)
      Enddo
      Enddo
      Enddo
      Fmi_a(m,i) = Fmi_a(m,i) + T
      Enddo
      Enddo

      Do I = 1, Nocc_a
      Do M = 1, Nocc_a
         T = 0.0D0
      Do E = 1, Nvrt_a
      Do f = 1, Nvrt_b
      Do n = 1, Nocc_b
         T = T + T1taup_ab(E,f,I,n)*W_ab(M,n,E,f)
      Enddo
      Enddo
      Enddo
      Fmi_a(M,I) =  Fmi_a(M,I) + T
      Enddo
      Enddo

C Fmi_b(m,i)= 1/2 T1taup(ef,in)*W2(mn,ef) + T1taup(eF,iN)*W(mN,eF)

      Do i = 1, Nocc_b
      Do m = 1, Nocc_b
         T = 0.0D0
      Do e = 1, Nvrt_b
      Do f = 1, Nvrt_b
      Do n = 1, Nocc_b
         T = T + 0.50D0*T1taup_bb(e,f,i,n)*W_bb(m,n,e,f)
      Enddo
      Enddo
      Enddo
      Fmi_b(m,i) = Fmi_b(m,i) + T
      Enddo
      Enddo

      Do i = 1, Nocc_b
      Do m = 1, Nocc_b
         T = 0.0D0
      Do e = 1, Nvrt_b
      Do F = 1, Nvrt_a
      Do N = 1, Nocc_a
         T = T + T1taup_ab(F,e,N,i)*W_ab(N,m,F,e)
      Enddo
      Enddo
      Enddo
      Fmi_b(m,i) = Fmi_b(m,i) + T
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      call checksum("Htmi_aa   :",Fmi_a,Nocc_a*Nocc_a)
      call checksum("Htmi_bb   :",Fmi_b,Nocc_b*Nocc_b)
#endif 

      Return 
      End 


