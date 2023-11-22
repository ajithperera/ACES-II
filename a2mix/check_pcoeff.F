      Subroutine  Check_pcoeff(PA,PB,Nrow,Ncol)

      Implicit Double Precision(A-H,O-Z)

      Dimension Pa(Nrow,Ncol), Pb(Nrow,Ncol)
      Dimension Pd(Nrow,Ncol)
 
      Call Dcopy(Nrow*Ncol,Pb,1,Pd,1)

      Call Daxpy(Nrow*Ncol,-1.0D0,Pa,1,Pd,1)
   
      Call checksum("PA-PB:",pd,Nrow*Ncol)
      Return
      End 
