      Subroutine Form_mnij(T2aa,T2bb,T2ab,T1aa,T1bb,W1_aa,W1_bb,
     +                     W1_ab,W_aa,W_bb,W_ab,W5_aa,W5_bb,W5_ab,
     +                     W5_ba,Wmnij_1,Wmnij_2,Wmnij_3,Nocc_a,
     +                     Nocc_b,Nvrt_a,Nvrt_b,NR)

      Implicit Double Precision (A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension W1_aa(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension W1_bb(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension W1_ab(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Dimension W5_aa(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension W5_bb(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension W5_ab(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension W5_ba(Nocc_a,Nocc_b,Nvrt_a,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension Wmnij_1(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension Wmnij_2(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension Wmnij_3(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Integer M,N,I,J,E,F
      Logical NR 

#include "ccsdlight_vars.com" 

#ifdef _DEBUG_LVM
      Write(6,*)
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("T2aa      :",T2aa,L_aaaa)
      call checksum("T2bb      :",T2bb,L_bbbb)
      call checksum("T2ab      :",T2ab,L_abab)
      L_aaaa = Nocc_a*Nocc_a*Nocc_a*Nocc_a
      L_bbbb = Nocc_b*Nocc_b*Nocc_b*Nocc_b
      L_abab = Nocc_a*Nocc_b*Nocc_a*Nocc_b
      call checksum("Form_wMNIJ:",W1_aa,L_aaaa)
      call checksum("Form_wmnij:",W1_bb,L_bbbb)
      call checksum("Form_wMnIj:",W1_ab,L_abab)
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("Form_wMNEF:",W_aa,L_aaaa)
      call checksum("Form_wmnef:",W_bb,L_bbbb)
      call checksum("Form_wMnEf:",W_ab,L_abab)
#endif

      IF (NR) then 
         Fact_aa = 0.50D0 
         Fact_bb = 0.50D0
         Fact_ab = 1.00D0
      Else
         Fact_aa = 0.250D0 
         Fact_bb = 0.250D0
         Fact_ab = 0.500D0
      Endif 

C Wmnij_1(MN,IJ) = W1(MN,IJ) + 1/4T2(IJ,EF)*W(MN,EF) 
 
      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do N = 1, Nocc_a
      Do M = 1, Nocc_a
         C = Ocn_oa(I)*Ocn_oa(J)
         Wmnij_1(M,N,I,J) = W1_aa(M,N,I,J)*C
      Do E = 1, Nvrt_a
      Do F = 1, Nvrt_a
CSSS      C = (1.0D0-Ocn_va(E))*(1.0D0-Ocn_va(F))
      C = 1.0D0
      Wmnij_1(M,N,I,J) = Wmnij_1(M,N,I,J) + 
     +                   Fact_aa*T2aa(E,F,I,J)*W_aa(M,N,E,F)*C
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C Wmnij_1(MN,IJ) = P_(IJ) T1(E,J)*W(MN,IE) 

      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do N = 1, Nocc_a
      Do M = 1, Nocc_a
         T = 0.0D0
      Do E = 1, Nvrt_a
CSSS         C = (1.0D0-Ocn_va(E))
         C = Ocn_oa(I)
         T = T + T1aa(E,J)*W5_aa(M,N,I,E)*C 
      Enddo
         Wmnij_1(M,N,I,J) =  Wmnij_1(M,N,I,J) + T 
         Wmnij_1(M,N,J,I) =  Wmnij_1(M,N,J,I) - T
      Enddo
      Enddo
      Enddo
      Enddo

C Wmnij_2(mn,ij) = W1(mn,ij) + 1/4T2(ij,ef)*W2(mn,ef) 

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
      Do N = 1, Nocc_b
      Do M = 1, Nocc_b
         C = Ocn_ob(i)*Ocn_ob(j)
         Wmnij_2(m,n,i,j) = W1_bb(m,n,i,j)*C
      Do E = 1, Nvrt_b
      Do F = 1, Nvrt_b
CSSS      C = (1.0D0-Ocn_vb(e))*(1.0D0-Ocn_vb(f))
      C = 1.0D0
      Wmnij_2(m,n,i,j) = Wmnij_2(m,n,i,j) + 
     +                   Fact_bb*T2bb(e,f,i,j)*W_bb(m,n,e,f)*C 
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C Wmnij_2(mn,ij) = P_(ij) T1(e,j)*W(mn,ie) 

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
      Do N = 1, Nocc_b
      Do M = 1, Nocc_b
         T = 0.0D0
      Do E = 1, Nvrt_b
CSSS         C = (1.0D0-Ocn_vb(e))
         C = Ocn_ob(i)
         T = T + T1bb(e,j)*W5_bb(m,n,i,e)*C 
      Enddo 
         Wmnij_2(m,n,i,j) = Wmnij_2(m,n,i,j) + T
         Wmnij_2(m,n,j,i) = Wmnij_2(m,n,j,i) - T
      Enddo
      Enddo
      Enddo
      Enddo

C Wmnij_2(Mn,Ij) = W2(Mn,Ij) + T2(Ij,Ef)*W2(Mn,Ef)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do N = 1, Nocc_b
      Do M = 1, Nocc_a
         C = Ocn_oa(I)*Ocn_ob(j)
         Wmnij_3(M,n,I,j) =  W1_ab(M,n,I,j)*C
      Do E = 1, Nvrt_a
      Do F = 1, Nvrt_b
CSSS      C = (1.0D0-Ocn_va(E))*(1.0D0-Ocn_vb(f))
      C = 1.0D0
      Wmnij_3(M,n,I,j) = Wmnij_3(M,n,I,j) + 
     +                   Fact_ab*T2ab(E,f,I,j)*W_ab(M,n,E,f)*C
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C WmniJ_3(M,n,I,j) = T1(e,j)*W(M,n,I,e) 
     
      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do N = 1, Nocc_b
      Do M = 1, Nocc_a
         T = 0.0D0
      Do E = 1, Nvrt_b
CSSS         C = (1.0D0-Ocn_vb(e))
         C = Ocn_oa(I)
         T = T + T1bb(e,j)*W5_ab(M,n,I,e)*C 
      Enddo
         Wmnij_3(M,n,I,j) = Wmnij_3(M,n,I,j) + T
      Enddo
      Enddo
      Enddo
      Enddo

C WmniJ_3(M,n,I,j) = T1(E,I)*W(M,n,E,j)
    
      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do N = 1, Nocc_b
      Do M = 1, Nocc_a
         T = 0.0D0
      Do E = 1, Nvrt_a
CSSS         C = (1.0D0-Ocn_va(E))
         C = Ocn_ob(j)
         T = T + T1aa(E,I)*W5_ba(M,n,E,j)*C
      Enddo
         Wmnij_3(M,n,I,j) = Wmnij_3(M,n,I,j) + T
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      L_aaaa = Nocc_a*Nocc_a*Nocc_a*Nocc_a
      L_bbbb = Nocc_b*Nocc_b*Nocc_b*Nocc_b
      L_abab = Nocc_a*Nocc_b*Nocc_a*Nocc_b
      call checksum("Form_wMNIJ:",Wmnij_1,L_aaaa)
      call checksum("Form_wmnij:",Wmnij_2,L_bbbb)
      call checksum("Form_wMnIj:",Wmnij_3,L_abab)
#endif
      
      Return
      End
