      Subroutine Form_htae(T1taup_aa,T1taup_bb,T1taup_ab,T1aa,T1bb,
     +                     Fae_a,Fae_b,W_aa,W_bb,W_ab,W4_aa,W4_bb,
     +                     W4_ab,W4_ba,Fockvv_od_a,Fockvv_od_b,
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
   
      Dimension W4_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W4_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W4_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W4_ba(Nvrt_a,Nvrt_b,Nocc_a,Nvrt_b)

      Dimension Fae_a(Nvrt_a,Nvrt_a)
      Dimension Fae_b(Nvrt_b,Nvrt_b)

      Dimension Fockov_a(Nocc_a,Nvrt_a)
      Dimension Fockov_b(Nocc_b,Nvrt_b)
 
      Dimension Fockvv_od_a(Nvrt_a,Nvrt_a)
      Dimension Fockvv_od_b(Nvrt_b,Nvrt_b)
   
      Integer M,N,A,F,E

#include "cse12_vars.com"

C Fae_a(A,E) <- +T1(F,M)*W(E,F,A,M) + T1(f,m)W(E,f,A,m) + (1=delta_ae)*fockvv(A,E)

      Do e = 1, Nvrt_a
      Do a = 1, Nvrt_a
         T = 0.0D0
         C = (1.0D0-Ocn_va(A))
         Fae_a(A,E) = Fockvv_od_a(A,E)*C
      Do m = 1, Nocc_a
      Do f = 1, Nvrt_a
         C = (1.0D0-Ocn_va(A))
         T = T + T1aa(F,M)*W4_aa(E,F,A,M)*C
      Enddo
      Enddo
         Fae_a(A,E) = Fae_a(A,E) + T
      Enddo
      Enddo

      Do e = 1, Nvrt_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do m = 1, Nocc_b
      Do f = 1, Nvrt_b
         C = (1.0D0-Ocn_va(A))
         T = T + T1bb(f,m)*W4_ab(E,f,A,m)*C
      Enddo
      Enddo
         Fae_a(A,E) = Fae_a(A,E) + T
      Enddo
      Enddo

C Fae_b(a,e) <- T1(f,m)*W(e,f,a,m) + T1(F,M)W(F,e,M,a) + (1=delta_ae)*fockvv(a,e))

      Do e = 1, Nvrt_b
      Do a = 1, Nvrt_b
         T = 0.0D0
         C = (1.0D0-Ocn_vb(a))
         Fae_b(a,e)  = Fockvv_od_b(a,e)*C
      Do f = 1, Nvrt_b
      Do m = 1, Nocc_b
         C = (1.0D0-Ocn_vb(a))
         T = T + T1bb(f,m)*W4_bb(e,f,a,m)*C
      Enddo
      Enddo
         Fae_b(a,e)  = Fae_b(a,e) + T
      Enddo
      Enddo

      Do e = 1, Nvrt_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do f = 1, Nvrt_a
      Do m = 1, Nocc_a
         C = (1.0D0-Ocn_vb(a))
         T = T + T1aa(F,M)*W4_ba(F,e,M,a)*C
      Enddo
      Enddo
         Fae_b(a,e)  = Fae_b(a,e) + T
      Enddo
      Enddo

C Fae_a(A,E) <- -1/2f(M,E)T(A,M)

      sum = 0.0D0
      Do e = 1, Nvrt_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do m = 1, Nocc_a
         T = T + 0.50D0*T1aa(A,M)*Fockov_a(M,E)
      Enddo
         Fae_a(A,E) = Fae_a(A,E) - T
      Enddo
      Enddo

C Fae_a(a,e) <- -1/2f(m,e)T(a,m)

      Do e = 1, Nvrt_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do m = 1, Nocc_b
         T = T + 0.50D0*T1bb(a,m)*Fockov_b(m,e)
      Enddo
         Fae_b(a,e) = Fae_b(a,e) - T
      Enddo
      Enddo

C Fae_a(A,E) <- -1/2 T1taup(AF,MN)*W(MN,EF) - T1taup(Af,Mn)*W(Mn,Ef)

      Do e = 1, Nvrt_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do m = 1, Nocc_a
      Do n = 1, Nocc_a
      Do f = 1, Nvrt_a
         T = T - 0.50D0*T1taup_aa(A,F,N,M)*W_aa(N,M,E,F)
      Enddo
      Enddo
      Enddo
         Fae_a(A,E) = Fae_a(A,E) + T
      Enddo
      Enddo

      Do e = 1, Nvrt_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do m = 1, Nocc_b
      Do n = 1, Nocc_a
      Do f = 1, Nvrt_b
         T = T - T1taup_ab(A,f,N,m)*W_ab(N,m,E,f)
      Enddo
      Enddo
      Enddo
         Fae_a(A,E) = Fae_a(A,E) + T
      Enddo
      Enddo

C Fae_b(a,e) <- -1/2 T1tau(af,mn)*W(mn,ef) - T1tau(mN,aF)*W(mN,eF)

      Do e = 1, Nvrt_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do m = 1, Nocc_b
      Do n = 1, Nocc_b
      Do f = 1, Nvrt_b
         T = T - 0.50D0*T1taup_bb(a,f,n,m)*W_bb(n,m,e,f)
      Enddo
      Enddo
      Enddo
         Fae_b(a,e)  = Fae_b(a,e) + T
      Enddo
      Enddo

      Do e = 1, Nvrt_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do m = 1, Nocc_a
      Do n = 1, Nocc_b
      Do f = 1, Nvrt_a
         T = T - T1taup_ab(F,a,M,n)*W_ab(M,n,F,e)
      Enddo
      Enddo
      Enddo
         Fae_b(a,e)  = Fae_b(a,e) + T
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      call checksum("Htae_aa   :",Fae_a,Nvrt_a*Nvrt_a)
      call checksum("Htae_bb   :",Fae_b,Nvrt_b*Nvrt_b)
#endif  

      Return 
      End 
