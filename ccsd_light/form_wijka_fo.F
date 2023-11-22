      Subroutine Form_wijka_fo(W_pqrs_aaaa,W_pqrs_bbbb,W_pqrs_aabb,
     +                         W5_aa,W5_bb,W5_ab,W5_ba,Nbasis,Nocc_a,
     +                         Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension W_pqrs_aaaa(Nbasis,Nbasis,Nbasis,Nbasis)
      Dimension W_pqrs_bbbb(Nbasis,Nbasis,Nbasis,Nbasis)
      Dimension W_pqrs_aabb(Nbasis,Nbasis,Nbasis,Nbasis)

      Dimension W5_aa(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension W5_bb(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension W5_ab(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension W5_ba(Nocc_a,Nocc_b,Nvrt_a,Nocc_b)

      Integer I,J,K,A,V,P,R,Q,S

      P = 1
      R = Nvrt_a 
      A = 0
      Do V = P, R
         A = A + 1
      Do K = 1, Nocc_a
      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
         WC = W_pqrs_aaaa(I,K,J,V)
         WE = W_pqrs_aaaa(J,K,I,V)
         W5_aa(I,J,K,A) = WC-WE
      Enddo
      Enddo
      Enddo
      Enddo

      Q = 1
      S = Nvrt_b
      A = 0
      Do V = Q, S 
         A = A +  1
      Do K = 1, Nocc_b
      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
         WC = W_pqrs_bbbb(i,k,j,v)
         WE = W_pqrs_bbbb(j,k,i,v)
         W5_bb(i,j,k,a) = WC-WE 
      Enddo
      Enddo
      Enddo
      Enddo

      A = 0
      Do V = Q, S
         A = A + 1
      Do K = 1, Nocc_a
      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
         W5_ab(I,j,K,a) = W_pqrs_aabb(I,K,j,v)
      Enddo
      Enddo
      Enddo
      Enddo

      A = 0
      Do V = P, R
         A = A + 1
      Do K = 1, Nocc_b
      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
         W5_ba(I,j,A,k) = W_pqrs_aabb(I,V,j,k)
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      L_aaaa = Nocc_a*Nocc_a*Nocc_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nocc_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nocc_a*Nvrt_b
      L_baba = Nocc_a*Nocc_b*Nocc_b*Nvrt_a
      call checksum("Form_wIJKA:",W5_aa,L_aaaa)
      call checksum("Form_wijka:",W5_bb,L_bbbb)
      call checksum("Form_wIjKa:",W5_ab,L_abab)
      call checksum("Form_wIjAk:",W5_ba,L_baba)
#endif

      Return
      End
