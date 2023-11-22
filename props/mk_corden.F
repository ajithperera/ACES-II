      Subroutine Mk_corden(Evecs,Dens,Nbfns,Naobfns,Docc)

      Implicit Double Precision(A-H,O-Z)

      Dimension Evecs(Naobfns,1),Dens(Naobfns,Naobfns)

      Data Ione,Done,Dnull /1,1.0D0,0.0D0/

      Call Xgemm("N","T",Naobfns,Naobfns,Ione,Docc,Evecs,Naobfns,
     &           Evecs,Naobfns,Dnull,Dens,Naobfns)
    
      Return
      End

