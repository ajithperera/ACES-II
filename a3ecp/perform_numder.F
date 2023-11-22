      Subroutine Perform_numder(X1,Y1,Z1,X2,Y2,Z2,Der2,Nbasis,Eps,
     &                          Ixyz,Jxyz)

      Implicit Double Precision(A-H,O-Z)

      Dimension X1(Nbasis,Nbasis) 
      Dimension Y1(Nbasis,Nbasis) 
      Dimension Z1(Nbasis,Nbasis) 
      Dimension X2(Nbasis,Nbasis) 
      Dimension Y2(Nbasis,Nbasis) 
      Dimension Z2(Nbasis,Nbasis) 
      Dimension Der2(Nbasis,Nbasis) 
      Logical XX,YY,ZZ,XY,XZ,YZ

      Eps2 = 2.0D0*Eps

      XX = .False.
      YY = .False.
      ZZ = .False.
      XY = .False.
      XZ = .False.
      YZ = .False.

      If (Ixyz .EQ. 1 .AND. Jxyz .EQ. 1) XX = .True.
      If (Ixyz .EQ. 2 .AND. Jxyz .EQ. 2) YY = .True.
      If (Ixyz .EQ. 3 .AND. Jxyz .EQ. 3) ZZ = .True.
      If (Ixyz .EQ. 1 .AND. Jxyz .EQ. 2) XY = .True.
      If (Ixyz .EQ. 1 .AND. Jxyz .EQ. 3) XZ = .True.
      If (Ixyz .EQ. 2 .AND. Jxyz .EQ. 3) YZ = .True.

      If (XX) Then
      DO I = 1, Nbasis 
      Do J = 1, Nbasis 
         Der2(I,J) = (X1(I,J) - X2(I,J))/Eps2
      Enddo 
      Enddo 

      Write(6,"(a)") " The XX derivative "
      Call output(Der2,1,Nbasis, 1,Nbasis,Nbasis,Nbasis,1)
      Endif

      If (YY) Then

      DO I = 1, Nbasis 
      Do J = 1, Nbasis 
         Der2(I,J) = (Y1(I,J) - Y2(I,J))/Eps2
      Enddo 
      Enddo 

      Write(6,"(a)") " The YY derivative "
      Call output(Der2,1,Nbasis, 1,Nbasis,Nbasis,Nbasis,1)
      Endif 

      If (ZZ) Then
      DO I = 1, Nbasis 
      Do J = 1, Nbasis 
         Der2(I,J) = (Z1(I,J) - Z2(I,J))/Eps2
      Enddo 
      Enddo 

      Write(6,"(a)") " The ZZ derivative "
      Call output(Der2,1,Nbasis, 1,Nbasis,Nbasis,Nbasis,1)
      Endif 

      If (XY) Then
      DO I = 1, Nbasis 
      Do J = 1, Nbasis 
         Der2(bI,J) = (X1(I,J) - Y2(I,J))/Eps2
      Enddo 
      Enddo 

      Write(6,"(a)") " The XY derivative "
      Call output(Der2,1,Nbasis, 1,Nbasis,Nbasis,Nbasis,1)
      Endif 

      If (XZ) Then
      DO I = 1, Nbasis 
      Do J = 1, Nbasis 
         Der2(I,J) = (X1(I,J) - Z2(I,J))/Eps2
      Enddo 
      Enddo 

      Write(6,"(a)") " The XZ derivative "
      Call output(Der2,1,Nbasis, 1,Nbasis,Nbasis,Nbasis,1)
      Endif 

      If (YZ) Then
      DO I = 1, Nbasis 
      Do J = 1, Nbasis 
         Der2(I,J) = (Y1(I,J) - Z2(I,J))/Eps2
      Enddo 
      Enddo 

      Write(6,"(a)") " The YZ derivative "
      Call output(Der2,1,Nbasis, 1,Nbasis,Nbasis,Nbasis,1)
      Endif 

      Return
      End

