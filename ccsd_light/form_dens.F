      Subroutine Form_dens(Dens,Evec,Scr,Occnum,Nbasis)

      Implicit Double Precision(A-H,O-Z)

      Dimension Dens(Nbasis,Nbasis),Ecve(Nbasis,Nbasis)
      Dimension Scr(Nbasis,Nbasis)
      Dimension Occnum(Nbasis,Nbasis)

      Call Dgemm("N","N",Nbasis,Nbasis,Nbasis,1.0D0,Evec,Nbasis,
     +            Occnum,Nbasis,0.0D0,Scr,Nbasis)

      Call Dgemm("N","T",Nbasis,Nbasis,Nbasis,1.0D0,Scr,Nbasis,
     +            Evec,Nbasis,0.0D0,Dens,Nbasis)

      Return
      End
