      Subroutine Form_wmnij(W_pqrs_aaaa,W_pqrs_bbbb,W_pqrs_aabb,
     +                     W1_aa,W1_bb,W1_ab,Nbasis,Nocc_a,
     +                     Nocc_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension W_pqrs_aaaa(Nbasis,Nbasis,Nbasis,Nbasis)
      Dimension W_pqrs_bbbb(Nbasis,Nbasis,Nbasis,Nbasis)
      Dimension W_pqrs_aabb(Nbasis,Nbasis,Nbasis,Nbasis)

      Dimension W1_aa(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension W1_bb(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension W1_ab(Nocc_a,Nocc_b,Nocc_a,Nocc_b)
      Integer M,N,I,J

      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do N = 1, Nocc_a
      Do M = 1, Nocc_a
         WC = W_pqrs_aaaa(M,I,N,J)
         WE = W_pqrs_aaaa(M,J,N,I)
         W1_aa(M,N,I,J) = WC - WE
CSSS        Write(6,"(4(1x,i4),1x,F15.10)") m,n,i,j,W1_aa(M,N,I,J) 
      Enddo
      Enddo
      Enddo
      Enddo

      Do j = 1, Nocc_b
      Do i = 1, Nocc_b
      Do n = 1, Nocc_b
      Do m = 1, Nocc_b
         WC = W_pqrs_bbbb(m,i,n,j)
         WE = W_pqrs_bbbb(m,j,n,i)
         W1_bb(m,n,i,j) = WC - WE 
CSSS        Write(6,"(4(1x,i4),1x,F15.10)") m,n,i,j,W1_bb(M,N,I,J) 

      Enddo
      Enddo
      Enddo
      Enddo

      Do j = 1, Nocc_b
      Do i = 1, Nocc_a
      Do n = 1, Nocc_b
      Do m = 1, Nocc_a
         W1_ab(M,n,I,j) = W_pqrs_aabb(M,I,n,j)
CSSS         Write(6,"(4(1x,i4),1x,F15.10)") m,n,i,j,W1_ab(M,N,I,J) 
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      L_aaaa = Nocc_a*Nocc_a*Nocc_a*Nocc_a
      L_bbbb = Nocc_b*Nocc_b*Nocc_b*Nocc_b
      L_abab = Nocc_a*Nocc_b*Nocc_a*Nocc_b
      call checksum("Form_wMNIJ:",W1_aa,L_aaaa)
      call checksum("Form_wmnij:",W1_bb,L_bbbb)
      call checksum("Form_wMnIj:",W1_ab,L_abab)
#endif

      Return
      End
