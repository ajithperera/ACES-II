      Subroutine Form_wabci_fo(W_pqrs_aaaa,W_pqrs_bbbb,W_pqrs_aabb,
     +                         W4_aa,W4_bb,W4_ab,W4_ba,Nbasis,Nocc_a,
     +                         Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension W_pqrs_aaaa(Nbasis,Nbasis,Nbasis,Nbasis)
      Dimension W_pqrs_bbbb(Nbasis,Nbasis,Nbasis,Nbasis)
      Dimension W_pqrs_aabb(Nbasis,Nbasis,Nbasis,Nbasis)

      Dimension W4_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W4_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W4_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W4_ba(Nvrt_a,Nvrt_b,Nocc_a,Nvrt_b)

      Integer P,Q,A,B,C,F,R,S,T,U,V,W,I

      P    = 1
      R    = Nvrt_a
      Do I = 1, Nocc_a 
         C = 0
      Do V = P, R
         B = 0 
         C = C + 1
      Do U = P, R
         A = 0
         B = B + 1
      Do T = P, R
         A = A + 1
         WC = W_pqrs_aaaa(T,V,U,I)
         WE = W_pqrs_aaaa(U,V,T,I)
         W4_aa(A,B,C,I) = WC-WE 
      Enddo
      Enddo
      Enddo
      Enddo

      Q    = 1
      S    = Nvrt_b
      Do I = 1, Nocc_b
         C = 0
         F = F + 1
      Do V = Q, S
         B = 0 
         C = C + 1
      Do U = Q, S
         A = 0
         B = B + 1
      Do T = Q, S
         A = A + 1
         WC = W_pqrs_bbbb(t,v,u,i)
         WE = W_pqrs_bbbb(u,v,t,i)
         W4_bb(a,b,c,i) = WC-WE
      Enddo
      Enddo
      Enddo
      Enddo
      

      Do I = 1, Nocc_b
         C = 0
      Do V = P, R
         B = 0 
         C = C + 1
      Do U = Q, S
         A = 0
         B = B + 1
      Do T = P, R 
         A = A + 1
         W4_ab(A,b,C,i) = W_pqrs_aabb(T,V,u,i)
      Enddo
      Enddo
      Enddo
      Enddo

      Do I = 1, Nocc_a
         C = 0
      Do V = Q, S
         B = 0
         C = C + 1
      Do U = Q, S
         A = 0
         B = B + 1
      Do T = P, R
         A = A + 1
         W4_ba(A,b,I,c) = W_pqrs_aabb(T,I,u,v)
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      L_aaaa = Nvrt_a*Nvrt_a*Nvrt_a*Nocc_a 
      L_bbbb = Nvrt_b*Nvrt_b*Nvrt_b*Nocc_b
      L_abab = Nvrt_a*Nvrt_b*Nvrt_a*Nocc_b
      L_baba = Nvrt_a*Nvrt_b*Nvrt_b*Nocc_a
      call checksum("Form_wABCI:",W4_aa,L_aaaa)
      call checksum("Form_wabci:",W4_bb,L_bbbb)
      call checksum("Form_wAbCi:",W4_ab,L_abab)
      call checksum("Form_wAbic:",W4_ba,L_baba)
#endif

      Return
      End
