      Subroutine Debugout(Dens_fao, Nbfns)
      Implicit Double Precision (A-H, O-Z)
      
      Dimension Dens_fao(Nbfns, Nbfns)
     
          Write(6, "(4(2x,F12.10))") ((Dens_fao(i,j),
     &                                 i=1,nbfns), j=1,nbfns)

      Return
      End
  
