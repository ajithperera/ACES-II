      Subroutine Print_hess(Hess,Nrow,Ncol)
     
      Implicit double precision(A-H,O-Z)

      Dimension Hess(Nrow, Ncol)

      Write(6,*) " Printing the Hessian"
      Call output(Hess,1,Nrow, 1, Ncol, Nrow, Ncol, 1)

      Return
      End
