      Subroutine Fbsum(Nbasis, Norbs, Sum, Dipole)
C
C Description
C    Calculate the FB localization sums : Sum{occupy}<i|ri|>**2
C 
C Variables
C    Nbasis    : Number of Basis functions
C    Norbs     : Number of occupied orbitals  
C    Sum       : Final sum returns in sum
C    Dipole    : MO dipole integrals
C
C Declarations
C 
      Implicit Double Precision (A-H, O-Z)
      Dimension Dipole(Nbasis, Nbasis)
C 
C Common block informations
C
      Common/Machsp/Iintln,Ifltln,Iintfp,Ialone,Ibitwd
      Common/Files/Luout,Moints
C 
C     
      Sum = 0.0D0
      Do 10 J = 1, Norbs
         Sum = Sum + Dipole(J, J)**2
 10   Continue
C
      Return
      End











