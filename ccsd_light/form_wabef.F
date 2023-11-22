      Subroutine Form_wabef(W_pqrs_aaaa,W_pqrs_bbbb,W_pqrs_aabb,
     +                      W2_aa,W2_bb,W2_ab,Nbasis,Nocc_a,
     +                      Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension W_pqrs_aaaa(Nbasis,Nbasis,Nbasis,Nbasis)
      Dimension W_pqrs_bbbb(Nbasis,Nbasis,Nbasis,Nbasis)
      Dimension W_pqrs_aabb(Nbasis,Nbasis,Nbasis,Nbasis)

      Dimension W2_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension W2_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension W2_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Integer P,Q,A,B,E,F,R,S,T,U,V,W

      P    = Nocc_a + 1
      R    = Nocc_a + Nvrt_a
      F    = 0
      Do W = P, R
         E = 0
         F = F + 1
      Do V = P, R
         B = 0 
         E = E + 1
      Do U = P, R
         A = 0
         B = B + 1
      Do T = P, R
         A = A + 1
         WC = W_pqrs_aaaa(T,V,U,W)
         WE = W_pqrs_aaaa(T,W,U,V)
         W2_aa(A,B,E,F) = WC-WE
      Enddo
      Enddo
      Enddo
      Enddo

      Q    = Nocc_b + 1
      S    = Nocc_b + Nvrt_b
      F    = 0
      Do W = Q, S
         E = 0
         F = F + 1
      Do V = Q, S
         B = 0 
         E = E + 1
      Do U = Q, S
         A = 0
         B = B + 1
      Do T = Q, S
         A = A + 1
         WC = W_pqrs_bbbb(t,v,u,w)
         WE = W_pqrs_bbbb(t,w,u,v)
         W2_bb(a,b,e,f) = WC-WE
      Enddo
      Enddo
      Enddo
      Enddo
      
      P    = Nocc_a + 1
      Q    = Nocc_b + 1
      R    = Nocc_a + Nvrt_a
      S    = Nocc_b + Nvrt_b
      F    = 0
      Do W = Q, S
         E = 0
         F = F + 1
      Do V = P, R
         B = 0 
         E = E + 1
      Do U = Q, S
         A = 0
         B = B + 1
      Do T = P, R 
         A = A + 1
         W2_ab(A,b,E,f) = W_pqrs_aabb(T,V,u,w)
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      L_aaaa = Nvrt_a*Nvrt_a*Nvrt_a*Nvrt_a
      L_bbbb = Nvrt_b*Nvrt_b*Nvrt_b*Nvrt_b
      L_abab = Nvrt_a*Nvrt_b*Nvrt_a*Nvrt_b
      call checksum("Form_wABEF:",W2_aa,L_aaaa)
      call checksum("Form_wabef:",W2_bb,L_bbbb)
      call checksum("Form_wAbEf:",W2_ab,L_abab)
#endif

      Return
      End
