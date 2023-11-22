      Subroutine Rdvprp(Nbasis, Lun, Filnam, Label, Dipole, Ierr, Buf,
     &                  Ibuf, Libuf)
C
C Description
C
C   Read specified property integral matrix from the disk.
C   Because integrals may break the symmetry of the molecule,
C   they are returned in full square-matrix.
C
C Variables
C
C   Nbasis  : Number of basis function
C   Lun     : Unit number corresponds to vprop integral file (input)
C   Filnam  : File name for property integrals (input)
C   Dipole  : Area to store dipole integrals from vprop (output)
C   Label   : Vprop record label for properties (input)
C   Ierr    : Return code
C             = 0 Succesful completion
C             = 1 unable to open/read file
C   Nbuf    : Number of integral in Vprop integral buffer
C   Buf     : Buffer for reading integrals (scratch)
C   Ibuf    : Buffer for reading integrals (scartch)
C
C Declarations
C
      Implicit Double Precision (A-H, O-Z)
      Character*(*) Filnam, Label
      Dimension Ibuf(Libuf), Buf(Libuf), Dipole(Nbasis,Nbasis)
C
C Common Block informations
C
      Common/Machsp/Iintln,Ifltln,Iintfp,Ialone,Ibitwd
      Common/Files/Luout,Moints

C Open the integral File
 
      Call Opnfil(Lun, Filnam, 'UNFORMATTED', 0, Ierr)
      If (Mod(Ierr, 2) .ne. 0) then
         Write (Luout, *) ' @Rdvprp-I, Unable to open file ',Filnam,'.'
         Return
      endif
C 
C Find the pertinent integrals
C
      Rewind (Lun)
      Call Filabl(Lun, Label, Ierr)
      If (Ierr .eq. 1) then
         Write(Luout, 1000) Label, Filnam, Lun
         Return
      Endif
 1000 Format (1X, ' @Rdvprp-I, Label,',A,', not found on file ', A,
     &        ' (unit ', I2,').')
C
      Call zero(Dipole, Nbasis*Nbasis)

C     Read the vprop integral file

 200  Read (Lun) Buf, Ibuf, Nbuf
      If (Nbuf .ne. -1) then
C
C Map the indices from triangular to square
C
         Do 10 J = 1, Nbuf
C              
C Ibuf carry the index corresponding to elements of the square matrix
C So m*(m-1)/2 is greater than Ibuf(J). 

            m = 1 + (-1 + INT(Dsqrt(8.D0*Ibuf(J)+0.999D0)))/2
C
C Now fairly easily n 
C
            n = Ibuf(J) - (m*(m-1))/2
            
            Dipole(n, m) = Buf(J)
            Dipole(m, n) = Buf(J)
 10      continue
         Goto 200
      Endif
         
C
      Close(Lun)
      Return
      End


