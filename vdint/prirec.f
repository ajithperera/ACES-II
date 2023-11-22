      SUBROUTINE PRIREC(OUTMAT,NDIM1,NDIM2,HEAD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
      CHARACTER HEAD*(*)
      DIMENSION OUTMAT(NDIM1,NDIM2)
C
      WRITE (LUPRI, '(//,2A)') '  ', HEAD
      WRITE (LUPRI, '(80A)')   '  ', ('-', I = 1, LEN(HEAD))
      WRITE (LUPRI, '()')
      DO 100 I = 1, NDIM1
         WRITE (LUPRI, '(6F12.6)') (OUTMAT(I,J), J = 1, NDIM2)
         WRITE (LUPRI, '()')
  100 CONTINUE
      RETURN
      END
