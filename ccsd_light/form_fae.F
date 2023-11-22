      Subroutine Form_fae(T2aa,T2bb,T2ab,T1aa,T1bb,Fae_a,Fae_b,W_aa,
     +                    W_bb,W_ab,W4_aa,W4_bb,W4_ab,W4_ba,Fockvv_od_a,
     +                    Fockvv_od_b,Fockov_a,Fockov_b,Nocc_a,
     +                    Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)
     
      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a) 
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

#include "ccsdlight_vars.com"

C Fae_aa(A,E)= -1/2 T2(AF,MN)*W(MN,EF) - T2(Af,Mn)*W(Mn,Ef)
      
      Do e = 1, Nvrt_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do m = 1, Nocc_a
      Do n = 1, Nocc_a
      Do f = 1, Nvrt_a
CSSS         C = (1.0D0-Ocn_va(F))*Ocn_oa(N)*Ocn_oa(M)
         C = 1.0D0
         T = T - 0.50D0*T2aa(A,F,N,M)*W_aa(N,M,E,F)*C
      Enddo
      Enddo
      Enddo 
      C1 = (1.0D0-Ocn_va(A))
      Fae_a(A,E) = T + Fockvv_od_a(A,E)*C1
      Enddo 
      Enddo 

      Do e = 1, Nvrt_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do m = 1, Nocc_b
      Do n = 1, Nocc_a
      Do f = 1, Nvrt_b
CSSS         C = (1.0D0-Ocn_vb(f))*Ocn_oa(N)*Ocn_ob(m)
         C = 1.0D0
         T = T - T2ab(A,f,N,m)*W_ab(N,m,E,f)*C
      Enddo
      Enddo
      Enddo 
         Fae_a(A,E) = Fae_a(A,E) + T
      Enddo 
      Enddo 

C Fae_bb(a,e)= -1/2 T2(af,mn)*W(mn,ef) - T2(mN,aF)*W(mN,eF)

      Do e = 1, Nvrt_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do m = 1, Nocc_b
      Do n = 1, Nocc_b
      Do f = 1, Nvrt_b
CSSS         C = (1.0D0-Ocn_vb(f))*Ocn_ob(n)*Ocn_ob(m)
         C = 1.0D0
         T = T - 0.50D0*T2bb(a,f,n,m)*W_bb(n,m,e,f)*C
      Enddo
      Enddo
      Enddo
      C1 = (1.0D0-Ocn_vb(a))
      Fae_b(a,e)  = T + Fockvv_od_b(a,e)*C1
      Enddo
      Enddo 

      Do e = 1, Nvrt_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do m = 1, Nocc_a
      Do n = 1, Nocc_b
      Do f = 1, Nvrt_a
CSSS         C = (1.0D0-Ocn_va(F))*Ocn_ob(n)*Ocn_oa(M)
         C = 1.0D0
         T = T - T2ab(F,a,M,n)*W_ab(M,n,F,e)*C
      Enddo
      Enddo
      Enddo
         Fae_b(a,e)  = Fae_b(a,e) + T
      Enddo
      Enddo 

C Fae_a(A,E)  = +T1(F,M)*W(E,F,A,M) + T1(f,m)W(E,f,A,m)

      Do e = 1, Nvrt_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do m = 1, Nocc_a
      Do f = 1, Nvrt_a
CSSS         C = (1.0D0-Ocn_va(F))*Ocn_oa(M)
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
CSSS         C = (1.0D0-Ocn_vb(f))*Ocn_ob(m)
         C = (1.0D0-Ocn_va(A))
         T = T + T1bb(f,m)*W4_ab(E,f,A,m)*C
      Enddo
      Enddo
         Fae_a(A,E) = Fae_a(A,E) + T
      Enddo
      Enddo

C Fae_b(a,e)  = T1(f,m)*W(e,f,a,m) + T1(F,M)W(F,e,M,a)

      Do e = 1, Nvrt_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do f = 1, Nvrt_b
      Do m = 1, Nocc_b
CSSS         C = (1.0D0-Ocn_vb(f))*Ocn_ob(m)
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
CSSS         C = (1.0D0-Ocn_va(F))*Ocn_oa(M)
         C = (1.0D0-Ocn_vb(a))
         T = T + T1aa(F,M)*W4_ba(F,e,M,a)*C
      Enddo
      Enddo
         Fae_b(a,e)  = Fae_b(a,e) + T
      Enddo
      Enddo

C Fae_a(A,E)  = -1/2f(M,E)T(A,M)

      sum = 0.0D0
      Do e = 1, Nvrt_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do m = 1, Nocc_a
CSSS         C = Ocn_oa(M)
         C = 1.0D0
         T = T + 0.50D0*T1aa(A,M)*Fockov_a(M,E)*C
      Enddo
         Fae_a(A,E) = Fae_a(A,E) - T
      Enddo
      Enddo

C Fae_a(a,e)  = -1/2f(m,e)T(a,m)

      Do e = 1, Nvrt_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do m = 1, Nocc_b
CSSS         C = Ocn_ob(m)
         C = 1.0D0
         T = T + 0.50D0*T1bb(a,m)*Fockov_b(m,e)*C
      Enddo
         Fae_b(a,e) = Fae_b(a,e) - T
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      call checksum("Fae_a     :",Fae_a,Nvrt_a*Nvrt_a)
      call checksum("Fae_b     :",Fae_b,Nvrt_b*Nvrt_b)
C      do e=1,Nvrt_a
C      do a=1,Nvrt_a
C      Write(6,"(2I2,1x,F15.10)") a,e,Fae_a(a,e)
C      enddo
C      enddo
#endif  

      Return 
      End 
