      Subroutine make_vnm1_pot(Oneh,Focka,Fockb,Iuhf,Ldim1,Nocc)

      Implicit Double precision (A-H,O-z)

      Dimension Oneh(Ldim1),Focka(Ldim1),Fockb(Ldim1)
      Dimension Nocc(8,2)

      Nocca = 0 
      Noccb = 0 
      Do I= 1, 8
         Nocca = Nocca + Nocc(i,1)
         Noccb = Noccb + Nocc(i,2)
      Enddo
#ifdef _DEBUG_LVLM
      Write(6,"(a,2I4)") "The number of occupied electrons: ",
     &                    Nocca,Noccb
#endif 
      Factora = Float(Nocca - 1)
      Factorb = Float(Noccb - 1)
      Call Dscal(Ldim1, Factora, Focka, 1)
      If (Iuhf .ne. 0) Call Dscal(Ldim1, Factorb, Fockb, 1)
      Call Daxpy(Ldim1, 1.0D0, Oneh, 1, Focka, 1)
      If (Iuhf .ne. 0) Call Daxpy(Ldim1, 1.0D0, Oneh, 1, Fockb, 1)

      Factora = (1.0D0/Nocca)
      Factorb = (1.0D0/Noccb)
      Call Dscal(Ldim1, Factora, Focka, 1)
      If (Iuhf .ne. 0) Call Dscal(Ldim1, Factorb, Fockb, 1)

#ifdef _DEBUG_LVLM
      Write(6,"(a,2I4)") "H + (N-1)/N Fock matrices"
      Write(6,"(5(1x,F15.7))") (Focka(i),i=1,Ldim1)
      If (Iuhf.Ne.0) Write(6,"(5(1x,F15.7))") (Fockb(i),i=1,Ldim1)
#endif
      Return
      End


