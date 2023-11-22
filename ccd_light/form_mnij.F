      Subroutine Form_mnij(T2aa,T2bb,T2ab,W1_aa,W1_bb,W1_ab,
     +                     W_aa,W_bb,W_ab,Wmnij_1,Wmnij_2,
     +                     Wmnij_3,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension W1_aa(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension W1_bb(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension W1_ab(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension Wmnij_1(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension Wmnij_2(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension Wmnij_3(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Integer M,N,I,J,E,F

#include "ccdlight_vars.com" 

#ifdef _DEBUG_LVLM
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

C Wmnij_1(MN,IJ) = W1(MN,IJ) + 1/4T2(IJ,EF)*W(MN,EF) 
 
      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do N = 1, Nocc_a
      Do M = 1, Nocc_a
         Wmnij_1(M,N,I,J) = W1_aa(M,N,I,J)
      Do E = 1, Nvrt_a
      Do F = 1, Nvrt_a
      C = (1.0D0-Ocn_va(E))*(1.0D0-Ocn_va(F))
      Wmnij_1(M,N,I,J) = Wmnij_1(M,N,I,J) + 
     +                   0.25D0*T2aa(E,F,I,J)*W_aa(M,N,E,F)*C
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C Wmnij_2(mn,ij) = W1(mn,ij) + 1/4T2(ij,ef)*W2(mn,ef) 

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
      Do N = 1, Nocc_b
      Do M = 1, Nocc_b
         Wmnij_2(m,n,i,j) = W1_bb(m,n,i,j)
      Do E = 1, Nvrt_b
      Do F = 1, Nvrt_b
      C = (1.0D0-Ocn_vb(e))*(1.0D0-Ocn_vb(f))
      Wmnij_2(m,n,i,j) = Wmnij_2(m,n,i,j) + 
     +                   0.25D0*T2bb(e,f,i,j)*W_bb(m,n,e,f)*C 
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C Wmnij_2(Mn,Ij) = W2(Mn,Ij) + T2(Ij,Ef)*W2(Mn,Ef)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do N = 1, Nocc_b
      Do M = 1, Nocc_a
         Wmnij_3(M,n,I,j) =  W1_ab(M,n,I,j)
      Do E = 1, Nvrt_a
      Do F = 1, Nvrt_b
      C = (1.0D0-Ocn_va(E))*(1.0D0-Ocn_vb(f))
      Wmnij_3(M,n,I,j) = Wmnij_3(M,n,I,j) + 
     +                   T2ab(E,f,I,j)*W_ab(M,n,E,f)*C
      Enddo
      Enddo
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
