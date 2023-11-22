      Subroutine Form_wabij(W_pqrs_aaaa,W_pqrs_bbbb,W_pqrs_aabb,
     +                      W0_aa,W0_bb,W0_ab,Nbasis,
     +                      Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension W_pqrs_aaaa(Nbasis,Nbasis,Nbasis,Nbasis)
      Dimension W_pqrs_bbbb(Nbasis,Nbasis,Nbasis,Nbasis)
      Dimension W_pqrs_aabb(Nbasis,Nbasis,Nbasis,Nbasis)

      Dimension W0_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension W0_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension W0_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Integer P,Q,A,B,I,J,R,S,C,D

      P    = Nocc_a + 1
      R    = Nvrt_a + Nocc_a
      S    = Nvrt_a + Nocc_a
      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
         B= 0
      Do D = P, R 
         A = 0 
         B = B + 1
      Do C = P, S
         A = A + 1
         WC = W_pqrs_aaaa(C,I,D,J)
         WE = W_pqrs_aaaa(C,J,D,I)
         W0_aa(A,B,I,J) = WC-WE
      Enddo
      Enddo
      Enddo
      Enddo

      Q    = Nocc_b + 1
      R    = Nvrt_b + Nocc_b
      S    = Nvrt_b + Nocc_b
      Do J = 1, Nocc_b
      Do I = 1, Nocc_b  
         B = 0
      Do D = Q, R
         A = 0 
         B = B + 1
      Do C = Q, S
         A = A + 1
         WC = W_pqrs_bbbb(c,i,d,j)
         WE = W_pqrs_bbbb(c,j,d,i)
         W0_bb(a,b,i,j) = WC-WE
      Enddo
      Enddo
      Enddo
      Enddo
      Write(6,*) 
 
      P    = Nocc_a + 1
      Q    = Nocc_b + 1
      R    = Nvrt_b + Nocc_b
      S    = Nvrt_a + Nocc_a
      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
         B = 0
      Do D = Q, R
         A = 0 
         B = B + 1
      Do C = P, S
         A = A + 1
         W0_ab(A,b,I,j) = W_pqrs_aabb(C,I,d,j)
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a 
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("Form_wABIJ:",W0_aa,L_aaaa)
      call checksum("Form_wabij:",W0_bb,L_bbbb)
      call checksum("Form_wAbIJ:",W0_ab,L_abab)
#endif 

      Return
      End
