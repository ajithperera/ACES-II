      Subroutine debug_mod2p(T2aa,T2bb,T2ab,Wmnij_1,Wmnij_2,Wmnij_3,
     +                       Wabef_1,Wabef_2,Wabef_3,W_aa,W_bb,W_ab,
     +                       Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension Wmnij_1(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension Wmnij_2(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension Wmnij_3(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Dimension Wabef_1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
      Dimension Wabef_2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
      Dimension Wabef_3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do N = 1, Nocc_a
      Do M = 1, Nocc_a
         T = 0.0D0
      Do E = 1, Nvrt_a
      Do F = 1, Nvrt_a
         T = T + 0.25D0*T2aa(E,F,I,J)*W_aa(M,N,E,F)
      Enddo
      Enddo
          Wmnij_1(M,N,I,J) = Wmnij_1(M,N,I,J) -T 
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
      Do N = 1, Nocc_b
      Do M = 1, Nocc_b
         T = 0.0D0
      Do E = 1, Nvrt_b
      Do F = 1, Nvrt_b
         T  = T + 0.25D0*T2bb(e,f,i,j)*W_bb(m,n,e,f)
      Enddo
      Enddo
         Wmnij_2(m,n,i,j) = Wmnij_2(m,n,i,j) - T
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do N = 1, Nocc_b
      Do M = 1, Nocc_a
         T = 0.0D0
      Do E = 1, Nvrt_a
      Do F = 1, Nvrt_b
         T = T + 0.50D0*T2ab(E,f,I,j)*W_ab(M,n,E,f)
      Enddo
      Enddo
        Wmnij_3(M,n,I,j) = Wmnij_3(M,n,I,j) - T
      Enddo
      Enddo
      Enddo
      Enddo

      Do F = 1, Nvrt_a
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
      Do N = 1, Nocc_a
         T = T + 0.25D0*T2aa(A,B,M,N)*W_aa(M,N,E,F)
      Enddo
      Enddo
         Wabef_1(A,B,E,F) = Wabef_1(A,B,E,F) - T
      Enddo
      Enddo
      Enddo
      Enddo

      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         T = 0.0D0
      Do M = 1, Nocc_b
      Do N = 1, Nocc_b
         T = T + 0.25D0*T2bb(a,b,m,n)*W_bb(m,n,e,f)
      Enddo
      Enddo
         Wabef_2(a,b,e,f) = Wabef_2(a,b,e,f) - T
      Enddo
      Enddo
      Enddo
      Enddo

      Do F = 1, Nvrt_b
      Do E = 1, Nvrt_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
      Do N = 1, Nocc_b
         T = T + 0.50D0*T2ab(A,b,M,n)*W_ab(M,n,E,f)
      Enddo
      Enddo
         Wabef_3(A,b,E,f) = Wabef_3(A,b,E,f) - T
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
      L_aaaa = Nvrt_a*Nvrt_a*Nvrt_a*Nvrt_a
      L_bbbb = Nvrt_b*Nvrt_b*Nvrt_b*Nvrt_b
      L_abab = Nvrt_a*Nvrt_b*Nvrt_a*Nvrt_b
      call checksum("Form_wABEF:",Wabef_1,L_aaaa)
      call checksum("Form_wabef:",Wabef_2,L_bbbb)
      call checksum("Form_wAbEf:",Wabef_3,L_abab)
#endif

      Return
      End 


