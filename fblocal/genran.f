      Subroutine genran(Test, Vecin, Ranum, Nalpa, Nbasis)
C
C Description 
C   Generate a random number. The present
C   version is obtained from GAMESS program
C   system
C 
C Local variables
C   Vecin  : SCF eigenvectors
C   Nalpa  : Number of ALPHA electrons
C   Nbasis : Number of basis functions
C   Ranum  : Random number (output)
C   U, XY  : Arbitrary labels
C
C Declarations
C  
      Implicit Double Precision (A-H, O-Z)
      Dimension Vecin(Nbasis, Nbasis)
      Parameter(Zero = 0.0D0, One = 1.0D0)
C
      Save U
C
      PI = Acos(-One)
C
      If (Test .ne. Zero) Then
         N = Abs(Nalpa - Nbasis) + 1
         M = N + 5
         XY = Vecin(N, M) * Atan(one)
         U  = (PI + XY)**5
         XY = Real(Int(U))
         U =  U - XY
         Ranum = U
      else
         U = (PI + U)**5
         XY = Real(Int(U))
         U = U - XY
         Ranum = U
      Endif
C
      Return
      End


         
