      Subroutine Form_wmbej(W_pqrs_aaaa,W_pqrs_bbbb,W_pqrs_aabb,
     +                      W3_aaaa,W3_bbbb,W3_abab,W3_baba,
     +                      W3_abba,W3_baab,Nbasis,Nocc_a,Nocc_b,
     +                      Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension W_pqrs_aaaa(Nbasis,Nbasis,Nbasis,Nbasis)
      Dimension W_pqrs_bbbb(Nbasis,Nbasis,Nbasis,Nbasis)
      Dimension W_pqrs_aabb(Nbasis,Nbasis,Nbasis,Nbasis)

      Dimension W3_aaaa(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W3_bbbb(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W3_abab(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W3_baba(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension W3_abba(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension W3_baab(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Integer P,Q,M,B,E,J,R,S,U,V

      P    = Nocc_a + 1
      R    = Nocc_a + Nvrt_a
      Do J = 1, Nocc_a
         E = 0
      Do V = P, R 
         B = 0
         E = E + 1
      Do U = P, R
         B = B + 1
      Do M = 1, Nocc_a
         WC = W_pqrs_aaaa(M,V,U,J)
         WE = W_pqrs_aaaa(M,J,U,V )
         W3_aaaa(M,B,E,J) = WC-WE
      Enddo
      Enddo
      Enddo
      Enddo

      Q    = Nocc_b + 1
      S    = Nocc_b + Nvrt_b
      Do J = 1, Nocc_b
         E = 0
      Do V = Q, S
         B = 0
         E = E + 1
      Do U = Q, S
         B = B + 1
      Do M = 1, Nocc_b
         WC = W_pqrs_bbbb(m,v,u,j)
         WE = W_pqrs_bbbb(m,j,u,v)
         W3_bbbb(m,b,e,j) = WC-WE
      Enddo
      Enddo
      Enddo
      Enddo
 
      P    = Nocc_a + 1
      Q    = Nocc_b + 1
      R    = Nocc_a + Nvrt_a
      S    = Nocc_b + Nvrt_b
      Do J = 1, Nocc_b
         E = 0
      Do V = P, R
         B = 0
         E = E + 1
      Do U = Q, S
         B = B + 1
      Do M = 1, Nocc_a
         W3_abab(M,b,E,j) = W_pqrs_aabb(M,V,u,j)
      Enddo
      Enddo
      Enddo
      Enddo

      P    = Nocc_a + 1
      Q    = Nocc_b + 1
      Do J = 1, Nocc_a
         E = 0
      Do V = Q, S
         B = 0
         E = E + 1
      Do U = P, R
         B = B + 1
      Do M = 1, Nocc_b
         W3_baba(m,B,e,J) = W_pqrs_aabb(U,J,m,v)
      Enddo
      Enddo
      Enddo
      Enddo

      Q    = Nocc_b + 1
      Do J = 1, Nocc_a
         E = 0
      Do V = Q, S
         B = 0
         E = E + 1
      Do U = Q, S
         B = B + 1
      Do M = 1, Nocc_a
         W3_abba(M,b,J,e) = W_pqrs_aabb(M,J,u,v)
      Enddo
      Enddo
      Enddo
      Enddo

      P    = Nocc_a + 1
      Do J = 1, Nocc_b
         E = 0
      Do V = P, R 
         B = 0
         E = E + 1
      Do U = P, R 
         B = B + 1
      Do M = 1, Nocc_b
         W3_baab(m,B,j,E) =  W_pqrs_aabb(U,V,m,j)
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      L_aaaa = Nocc_a*Nvrt_a*Nvrt_a*Nocc_a
      L_bbbb = Nocc_b*Nvrt_b*Nvrt_b*Nocc_b
      L_abab = Nocc_a*Nvrt_b*Nvrt_a*Nocc_b
      L_baba = Nocc_b*Nvrt_a*Nvrt_b*Nocc_a
      L_abba = Nocc_a*Nvrt_b*Nocc_a*Nvrt_b
      L_baab = Nocc_b*Nvrt_a*Nocc_b*Nvrt_a 
      call checksum("Form_wMBEJ:",W3_aaaa,L_aaaa)
      call checksum("Form_wmbej:",W3_bbbb,L_bbbb)
      call checksum("Form_wMbEj:",W3_abab,L_abab)
      call checksum("Form_wmBeJ:",W3_baba,L_baba)
      call checksum("Form_wMbeJ:",W3_abba,L_abba)
      call checksum("Form_wmBEj:",W3_baab,L_baab)
#endif

      Return
      End
