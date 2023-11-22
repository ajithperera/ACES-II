      subroutine Form_fmi(T2aa,T2bb,T2ab,T1aa,T1bb,Fmi_a,Fmi_b,W_aa,
     +                    W_bb,W_ab,W5_aa,W5_bb,W5_ab,W5_ba,
     +                    Fockoo_od_a,Fockoo_od_b,Fockov_a,Fockov_b,
     +                    Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b) 
     
      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_b) 
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b) 
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b) 

      Dimension W5_aa(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension W5_bb(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension W5_ab(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension W5_ba(Nocc_a,Nocc_b,Nvrt_a,Nocc_b)

      Dimension Fockoo_od_a(Nocc_a,Nocc_a)
      Dimension Fockoo_od_b(Nocc_b,Nocc_b)

      Dimension Fmi_a(Nocc_a,Nocc_a)
      Dimension Fmi_b(Nocc_b,Nocc_b)

      Dimension Fockov_a(Nocc_a,Nvrt_a)
      Dimension Fockov_b(Nocc_b,Nvrt_b)

      Integer M,N,E,F,I

#include "ccsdlight_vars.com"

#ifdef _DEBUG_LVLM
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("T2aa      :",T2aa,L_aaaa)
      call checksum("T2bb      :",T2bb,L_bbbb)
      call checksum("T2ab      :",T2ab,L_abab)
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("Form_wMNEF:",W_aa,L_aaaa)
      call checksum("Form_wmnef:",W_bb,L_bbbb)
      call checksum("Form_wMnEf:",W_ab,L_abab)
#endif 

C Fmi_aa(M,I)= 1/2 T2(EF,IN)*W2(MN,EF) + T2(Ef,In)*W(Mn,Ef)
      
      Do I = 1, Nocc_a
      Do M = 1, Nocc_a
         T = 0.0D0 
      Do E = 1, Nvrt_a
      Do F = 1, Nvrt_a
      Do N = 1, Nocc_a
CSSS         C = Ocn_oa(N)*(1.0D0-Ocn_va(F))*(1.0D0-Ocn_va(E))
         C = 1.0D0
         T = T + 0.50D0*T2aa(E,F,I,N)*W_aa(M,N,E,F)*C 
      Enddo
      Enddo
      Enddo 
      C1 = Ocn_oa(I)
      Fmi_a(M,I) = T + Fockoo_od_a(M,I)*C1
      Enddo 
      Enddo 
 
      Do I = 1, Nocc_a
      Do M = 1, Nocc_a
         T = 0.0D0 
      Do E = 1, Nvrt_a
      Do f = 1, Nvrt_b
      Do n = 1, Nocc_b
CSSS         C = Ocn_ob(n)*(1.0D0-Ocn_vb(f))*(1.0D0-Ocn_va(E))
         C = 1.0D0
         T = T + T2ab(E,f,I,n)*W_ab(M,n,E,f)*C
      Enddo
      Enddo
      Enddo 
      Fmi_a(M,I) =  Fmi_a(M,I) + T
      Enddo 
      Enddo 

C Fmi_bb(m,i)= 1/2 T2(ef,in)*W2(mn,ef) + T2(eF,iN)*W(mN,eF)

      Do i = 1, Nocc_b
      Do m = 1, Nocc_b
         T = 0.0D0
      Do e = 1, Nvrt_b
      Do f = 1, Nvrt_b
      Do n = 1, Nocc_b
CSSS         C = Ocn_ob(n)*(1.0D0-Ocn_vb(f))*(1.0D0-Ocn_vb(e))
         C = 1.0D0
         T = T + 0.50D0*T2bb(e,f,i,n)*W_bb(m,n,e,f)*C 
      Enddo
      Enddo
      Enddo
      C1 = Ocn_ob(i)
      Fmi_b(m,i) = T +  Fockoo_od_b(m,i)*C1
      Enddo
      Enddo 

      Do i = 1, Nocc_b
      Do m = 1, Nocc_b
         T = 0.0D0
      Do e = 1, Nvrt_b
      Do F = 1, Nvrt_a
      Do N = 1, Nocc_a
CSSS         C = Ocn_oa(N)*(1.0D0-Ocn_va(F))*(1.0D0-Ocn_vb(e))
         C = 1.0D0
         T = T + T2ab(F,e,N,i)*W_ab(N,m,F,e)*C
      Enddo
      Enddo
      Enddo
      Fmi_b(m,i) = Fmi_b(m,i) + T 
      Enddo
      Enddo 

C  Fmi_a(M,I) = T1(E,N)*W(M,N,I,E) + T1(e,n)*W(M,n,I,e)

      Do i = 1, Nocc_a
      Do m = 1, Nocc_a
         T = 0.0D0
      Do e = 1, Nvrt_a
      Do N = 1, Nocc_a
CSSS         C = Ocn_oa(N)*(1.0D0-Ocn_va(E))
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
CSSS         C = Ocn_ob(n)*(1.0D0-Ocn_vb(e))
         C = Ocn_oa(I)
         T = T + T1bb(e,n)*W5_ab(M,n,I,e)*C 
      Enddo
      Enddo
      Fmi_a(M,I) = Fmi_a(M,I) + T
      Enddo
      Enddo

C  Fmi_b(m,i) = T1(e,n)*W(m,n,i,e) + T1(E,N)*W(N,m,E,i)


      Do i = 1, Nocc_b
      Do m = 1, Nocc_b
         T = 0.0D0
      Do e = 1, Nvrt_b
      Do N = 1, Nocc_b
CSSS         C = Ocn_ob(n)*(1.0D0-Ocn_vb(e))
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
CSSS         C = Ocn_oa(N)*(1.0D0-Ocn_va(E))
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
CSSS         C = (1.0D0-Ocn_va(E))
         C = 1.0D0
         T = T + 0.50D0*T1aa(E,I)*Fockov_a(M,E)*C
      Enddo
      Fmi_a(m,i) = Fmi_a(m,i) + T
      Enddo
      Enddo

C  Fmi_b(m,i) =  +1/2f(m,e)*T1(e,i)

      Do i = 1, Nocc_b
      Do m = 1, Nocc_b
         T = 0.0D0
      Do e = 1, Nvrt_b
CSSS         C = (1.0D0-Ocn_vb(e))
         C = 1.0D0
         T = T + 0.50D0*T1bb(e,i)*Fockov_b(m,e)*C
      Enddo
      Fmi_b(m,i) = Fmi_b(m,i) + T
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      call checksum("Fmi_a     :",Fmi_a,Nocc_a*Nocc_a)
      call checksum("Fmi_b     :",Fmi_b,Nocc_b*Nocc_b)
#endif 

      Return 
      End 
