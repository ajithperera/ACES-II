      Subroutine Mkden(Evecs,Dens,Nbfns)

      Implicit Double Precision(A-H,O-Z)

      Dimension Evecs(Nbfns,1),Dens(Nbfns,Nbfns)

      Data Ione,Done,Dnull /1,1.0D0,0.0D0/

      Call Xgemm("N","T",Nbfns,Nbfns,Ione,Done,Evecs,Nbfns,
     &           Evecs,Nbfns,Dnull,Dens,Nbfns)
    
      Return
      End

