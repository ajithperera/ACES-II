      SUBROUTINE DIAGMAT(EA,WORKEA,EIVELEA,EIVEREA,EIVAEA,WIEA,LWORKEA
     &,NEA,P)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
C
      DIMENSION EA(NEA,NEA)
      DIMENSION WORKEA(LWORKEA)
      DIMENSION EIVELEA(NEA,NEA)
      DIMENSION EIVEREA(NEA,NEA)
      DIMENSION EIVAEA(NEA)
      DIMENSION WIEA(NEA)
      LOGICAL P
C
      DONE = 1.0D+0
      DZERO = 0.D+0  
C
      CALL DGEEV('V','V',NEA,EA,NEA,EIVAEA,WIEA,EIVELEA,NEA,EIVEREA,NEA,
     &WORKEA,LWORKEA,INFO)
C
      IF (INFO.NE.0) THEN
         WRITE(6,*) 'ERROR IN DIAGONALIZATION'
         WRITE(6,*) 'INFO=',INFO
      ENDIF
C
      CALL DLASRT('I',NEA,EIVAEA,INFO)
C
      IF (P) THEN
         WRITE(6,*)
         WRITE(6,*) 'EIGENVALUES'
         CALL OUTPUT(EIVAEA,1,NEA,1,1,NEA,1,1)
      ENDIF
C
      RETURN
      END