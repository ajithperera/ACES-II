      Subroutine Form_wmnef(W_pqrs_aaaa,W_pqrs_bbbb,W_pqrs_aabb,
     +                     W_aa,W_bb,W_ab,Nbasis,Nocc_a,Nocc_b,
     +                     Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension W_pqrs_aaaa(Nbasis,Nbasis,Nbasis,Nbasis)
      Dimension W_pqrs_bbbb(Nbasis,Nbasis,Nbasis,Nbasis)
      Dimension W_pqrs_aabb(Nbasis,Nbasis,Nbasis,Nbasis)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Integer M,N,E,F,P,Q,R,S,T,U

      P = Nocc_a + 1
      R = Nocc_a + Nvrt_a 
      F = 0
      Do T = P, R
         E = 0
         F = F + 1
      Do U = P, R 
         E = E + 1
      Do N = 1, Nocc_a
      Do M = 1, Nocc_a
         WC = W_pqrs_aaaa(M,U,N,T)
         WE = W_pqrs_aaaa(M,T,N,U)
         W_aa(M,N,E,F) = WC-WE
      Enddo
      Enddo
      Enddo
      Enddo

      Q = Nocc_b + 1
      S = Nocc_b + Nvrt_b 
      F = 0
      Do T = Q, S
         E = 0
         F = F + 1
      Do U = Q, S
         E = E + 1
      Do N = 1, Nocc_b
      Do M = 1, Nocc_b
         WC = W_pqrs_bbbb(m,u,n,t)
         WE = W_pqrs_bbbb(m,t,n,u)
         W_bb(m,n,e,f) = WC-WE
      Enddo
      Enddo
      Enddo
      Enddo

      P = Nocc_a + 1
      Q = Nocc_b + 1
      F = 0
      Do T = Q, S
         E = 0
         F = F + 1
      Do U = P, R
         E = E + 1
      Do N = 1, Nocc_b
      Do M = 1, Nocc_a
         W_ab(M,n,E,f)= W_pqrs_aabb(M,U,n,t)
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("Form_wMNEF:",W_aa,L_aaaa)
      call checksum("Form_wmnef:",W_bb,L_bbbb)
      call checksum("Form_wMnEf:",W_ab,L_abab)
#endif

      Return
      End
