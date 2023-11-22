      Subroutine debug_modf(T2aa,T2bb,T2ab,Fmi_a,Fmi_b,Fae_a,Fae_b,
     +                     W_aa,W_bb,W_ab,Nocc_a, Nocc_b,Nvrt_a,
     +                     Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension Fmi_a(Nocc_a,Nocc_a)
      Dimension Fmi_b(Nocc_b,Nocc_b)

      Dimension Fae_a(Nvrt_a,Nvrt_a)
      Dimension Fae_b(Nvrt_b,Nvrt_b)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_B,Nocc_a,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Do I = 1, Nocc_a
      Do M = 1, Nocc_a
         T = 0.0D0
      Do E = 1, Nvrt_a
      Do F = 1, Nvrt_a
      Do N = 1, Nocc_a
         T = T + 0.50D0*T2aa(E,F,I,N)*W_aa(M,N,E,F)
      Enddo
      Enddo
      Enddo
      Fmi_a(M,I) = Fmi_a(M,I) - T
      Enddo
      Enddo

      Do I = 1, Nocc_a
      Do M = 1, Nocc_a
         T = 0.0D0
      Do E = 1, Nvrt_a
      Do f = 1, Nvrt_b
      Do n = 1, Nocc_b
         T = T + T2ab(E,f,I,n)*W_ab(M,n,E,f)
      Enddo
      Enddo
      Enddo
      Fmi_a(M,I) =  Fmi_a(M,I) - T
      Enddo
      Enddo


      Do i = 1, Nocc_b
      Do m = 1, Nocc_b
         T = 0.0D0
      Do e = 1, Nvrt_b
      Do f = 1, Nvrt_b
      Do n = 1, Nocc_b
         T = T + 0.50D0*T2bb(e,f,i,n)*W_bb(m,n,e,f)
      Enddo
      Enddo
      Enddo
      Fmi_b(m,i) = Fmi_b(m,i) - T
      Enddo
      Enddo

      Do i = 1, Nocc_b
      Do m = 1, Nocc_b
         T = 0.0D0
      Do e = 1, Nvrt_b
      Do F = 1, Nvrt_a
      Do N = 1, Nocc_a
         T = T + T2ab(F,e,N,i)*W_ab(N,m,F,e)
      Enddo
      Enddo
      Enddo
      Fmi_b(m,i) = Fmi_b(m,i) - T
      Enddo
      Enddo

      Do e = 1, Nvrt_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do m = 1, Nocc_a
      Do n = 1, Nocc_a
      Do f = 1, Nvrt_a
         T = T - 0.50D0*T2aa(A,F,N,M)*W_aa(N,M,E,F)
      Enddo
      Enddo
      Enddo
      Fae_a(A,E) =  Fae_a(A,E) - T
      Enddo
      Enddo

      Do e = 1, Nvrt_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do m = 1, Nocc_b
      Do n = 1, Nocc_a
      Do f = 1, Nvrt_b
         T = T - T2ab(A,f,N,m)*W_ab(N,m,E,f)
      Enddo
      Enddo
      Enddo
         Fae_a(A,E) = Fae_a(A,E) - T
      Enddo
      Enddo

      Do e = 1, Nvrt_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do m = 1, Nocc_b
      Do n = 1, Nocc_b
      Do f = 1, Nvrt_b
         T = T - 0.50D0*T2bb(a,f,n,m)*W_bb(n,m,e,f)
      Enddo
      Enddo
      Enddo
      Fae_b(a,e)  =  Fae_b(a,e) - T 
      Enddo
      Enddo

      Do e = 1, Nvrt_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do m = 1, Nocc_a
      Do n = 1, Nocc_b
      Do f = 1, Nvrt_a
         T = T - T2ab(F,a,M,n)*W_ab(M,n,F,e)
      Enddo
      Enddo
      Enddo
         Fae_b(a,e)  = Fae_b(a,e) - T
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      call checksum("Fmi_a     :",Fmi_a,Nocc_a*Nocc_a)
      call checksum("Fmi_b     :",Fmi_b,Nocc_b*Nocc_b)
      call checksum("Fae_a     :",Fae_a,Nvrt_a*Nvrt_a)
      call checksum("Fae_b     :",Fae_b,Nvrt_b*Nvrt_b)
#endif

      Return
      End 


