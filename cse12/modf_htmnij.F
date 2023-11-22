      Subroutine Modf_htmnij(T1aa,T1bb,W1_aa,W1_bb,W1_ab,W5_aa,W5_bb,
     +                       W5_ab,W5_ba,Wmnij_1,Wmnij_2,Wmnij_3,
     +                       Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension W1_aa(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension W1_bb(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension W1_ab(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Dimension W5_aa(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension W5_bb(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension W5_ab(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension W5_ba(Nocc_a,Nocc_b,Nvrt_a,Nocc_b)

      Dimension Wmnij_1(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension Wmnij_2(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension Wmnij_3(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Integer M,N,I,J,E,F

C These simply subtract these contributions since these terms have 
C have now twice what they suppose to be 

#include "cse12_vars.com" 

C Wmnij_1(MN,IJ) <- W1(MN,IJ)

      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do N = 1, Nocc_a
      Do M = 1, Nocc_a
         T = 0.0D0
         C = Ocn_oa(I)*Ocn_oa(J)
         Wmnij_1(M,N,I,J) = Wmnij_1(M,N,I,J) - W1_aa(M,N,I,J)*C
      Enddo
      Enddo
      Enddo
      Enddo

C Wmnij_1(MN,IJ) <- P_(IJ) T1(E,J)*W(MN,IE) 

      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do N = 1, Nocc_a
      Do M = 1, Nocc_a
         T = 0.0D0
      Do E = 1, Nvrt_a
         C = Ocn_oa(I)
         T = T + T1aa(E,J)*W5_aa(M,N,I,E)*C
      Enddo
         Wmnij_1(M,N,I,J) =  Wmnij_1(M,N,I,J) - T
         Wmnij_1(M,N,J,I) =  Wmnij_1(M,N,J,I) + T
      Enddo
      Enddo
      Enddo
      Enddo

C Wmnij_2(mn,ij) <- W1(mn,ij)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
      Do N = 1, Nocc_b
      Do M = 1, Nocc_b
         C = Ocn_ob(i)*Ocn_ob(j)
         Wmnij_2(m,n,i,j) = Wmnij_2(m,n,i,j) - W1_bb(m,n,i,j)*C
      Enddo
      Enddo
      Enddo
      Enddo

C Wmnij_2(mn,ij) <- P_(ij) T1(e,j)*W(mn,ie) 

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
      Do N = 1, Nocc_b
      Do M = 1, Nocc_b
         T = 0.0D0
      Do E = 1, Nvrt_b
         C = Ocn_ob(i)
         T = T + T1bb(e,j)*W5_bb(m,n,i,e)*C
      Enddo
         Wmnij_2(m,n,i,j) = Wmnij_2(m,n,i,j) - T
         Wmnij_2(m,n,j,i) = Wmnij_2(m,n,j,i) + T
      Enddo
      Enddo
      Enddo
      Enddo

C WmniJ_3(M,n,I,j) <- W1(Mn,Ij)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do N = 1, Nocc_b
      Do M = 1, Nocc_a
         C = Ocn_oa(I)*Ocn_ob(j)
         Wmnij_3(M,n,I,j) = Wmnij_3(M,n,I,j) - W1_ab(M,n,I,j)*C
      Enddo
      Enddo
      Enddo
      Enddo

C WmniJ_3(M,n,I,j) <- T1(E,I)*W(M,n,E,j)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do N = 1, Nocc_b
      Do M = 1, Nocc_a
         T = 0.0D0
      Do E = 1, Nvrt_b
         C = Ocn_oa(I)
         T = T + T1bb(e,j)*W5_ab(M,n,I,e)*C
      Enddo
         Wmnij_3(M,n,I,j) = Wmnij_3(M,n,I,j) - T
      Enddo
      Enddo
      Enddo
      Enddo

C WmniJ_3(M,n,I,j) <- T1(E,I)*W(M,n,E,j)
   
      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do N = 1, Nocc_b
      Do M = 1, Nocc_a
         T = 0.0D0
      Do E = 1, Nvrt_a
         C = Ocn_ob(j)
         T = T + T1aa(E,I)*W5_ba(M,n,E,j)*C
      Enddo
         Wmnij_3(M,n,I,j) = Wmnij_3(M,n,I,j) - T
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " @-modf_htmnij"
      L_aaaa = Nocc_a*Nocc_a*Nocc_a*Nocc_a
      L_bbbb = Nocc_b*Nocc_b*Nocc_b*Nocc_b
      L_abab = Nocc_a*Nocc_b*Nocc_a*Nocc_b
      call checksum("New_wMNIJ :",Wmnij_1,L_aaaa)
      call checksum("New_wmnij :",Wmnij_2,L_bbbb)
      call checksum("New_wMnIj :",Wmnij_3,L_abab)
#endif
      
      Return
      End
