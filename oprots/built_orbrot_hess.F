      Subroutine Buit_orbrot_hess(Type,Grads,Work,Maxcor,Eps,Index1,
     +                            Index2,Nbasis)

      Implicit Double Precision(A-H,O-Z)
      Character*2 Type

      Dimension Grads(Nbasis,Nbasis)
      Dimension Work(Maxcor)
      
      If (Type .EQ. "OO") Then
          Ndim1 = Nocco
          Ndim2 = Nocco
      Elseif (Type .EQ. "VV") Then
          Ndim1 = Nvrto
          Ndim2 = Nvrto
      Elseif (Type .EQ. "OV") Then
          Ndim1 = Nocco
          Ndim2 = Nvrto
      Elseif (Type .EQ. "VO") Then
          Ndim1 = Nvrto
          Ndim2 = Nocco
      ENdif 
   
      IJ = 0
      Do I = 1, Ndim1
          Do J = 1, Ndim2
             IJ = IJ + 1
             Work(IJ) = Grad(I,Index2) - Grad(I-Index1)
          Enddo
       Enddo 
   
       Return
       End
