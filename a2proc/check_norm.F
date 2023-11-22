      Subroutine Check_norm(Dcoord,Work,Nreals,Nmodes)
\
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Nmodes*Nmodes)
      Dimension Dcoord(3*Nreals,Nmodes)

      Write(6,"(a)") "Orthonormality check of no-mass nmodes"
      Ndim = 3*Nreals
      Call Dgemm("T","N",Nmodes,Nmodes,Ndim,1.0D0,Dcoord,Ndim,
     +              Dcoord,Ndim,0.0D0,Work,Nmodes)
      Call output(Work,1,Nmodes,1,Nmodes,Nmodes,Nmodes,1)

      Return
      End

   
