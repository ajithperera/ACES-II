      Subroutine MkFock(Iuhf, Label, Nbasis, Nbasmo,FockAO, Scr1,
     &                  Scr2, FockLO, FBLmo)
C 
C Description
C
C   Make the fock-matrix from FB localized molecular
C   orbitals. Note that resulting Fock matrix is 
C   non-diagonal. (Warnning: Make sure to use the 
C   proper flags when such Fock matrix is used in 
C   subsequent correlated calculation.
C
C Variables
C   Iuhf       : Flag for UHF and RHF runs
C   Label      : Distinguish ALPHA and BETA spin cases
C   Nbasis     : Number of symetry adapted basis functions 
C   Nbasmo     : Number of molecular orbitals
C   FockAO     : AO basis Fock matrix 
C   Scr1, Scr2 : Scratch Arrays
C   FockLO     : Fock matirx in FB localized molecular orbital basis
C   FBLmo      : FB localized orbitals
C
C Declarations
      Implicit Double Precision (A-H, O-Z)
      Dimension FockAO(Nbasis, Nbasis), Scr1(Nbasis, Nbasis),
     &          Scr2(Nbasis, Nbasis), FockLO(Nbasis, Nbasis),
     &          FBLmo(Nbasis, Nbasis)
      Character*5 Label
C
C Common Block informations
C
      Common/Machsp/Iintln,Ifltln,Iintfp,Ialone,Ibitwd
      Common/Files/Luout, Moints
C      
      Ione  = 1
      Zlich = 0.0D+00
      One   = 1.0D+00
C
C Write the FB localized orbitals into "JOBARC" file
C 
      IF (Label .eq. 'ALPHA') then
         Call Putrec(20, 'JOBARC', 'SCFEVECA', Nbasis*Nbasis*Iintfp,
     &               FBLmo)
      Else
         Call Putrec(20, 'JOBARC', 'SCFEVECB', Nbasis*Nbasis*Iintfp,
     &               FBLmo)
      Endif
C
C Initialize the all the arrays to zero. This is
C not really necessary but it is better to be safe, always.
C
      Call Zero (FockAO, Nbasis*Nbasis)
      Call Zero (Scr2, Nbasis*Nbasis)
C
      If (Label .eq. 'ALPHA') Then
         Call Getrec(20, 'JOBARC', 'FOCKA   ', Nbasis*Nbasis*Iintfp,
     &               FockAO)
      Else
         Call Getrec(20, 'JOBARC', 'FOCKB   ', Nbasis*Nbasis*Iintfp,
     &               FockAO)
      Endif
C
C Build the Fock matrix in FB-localized orbital basis
C The resulting Fock matrix is not diagonal

      Call Xgemm('N', 'N', Nbasis, Nbasis, Nbasis, One, 
     &     FockAO, Nbasis, FBLmo, Nbasis, Zlich, Scr1, Nbasis)
      
      Call Xgemm('T', 'N', Nbasis, Nbasis, Nbasis, One
     &     FBLmo, Nbasis, Scr1, Nbasis, Zlich, FockLO, Nbasis)
C 
C Copy the diagonal elements of transformed matrix to Scr2
C array
C 
      Call Dcopy(Nbasis, FockLO, Nbasis + 1, Scr2, 1)
C 
      If (Label .eq. 'ALPHA') Then
         Call Putrec(20, 'JOBARC', 'SCFEVALA', Nbasis*Iintfp, Scr2)
      Else
         Call Putrec(20, 'JOBARC', 'SCFEVALB', Nbasis*Iintfp, Scr2)
      Endif
C
C Compute the reference energy (should be identical to
C the HF energy)
C     
      Return
      End
    
      
