      SUBROUTINE PRITRI(OUTMAT,NDIM,HEAD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
      CHARACTER HEAD*(*)
      DIMENSION OUTMAT(NDIM*(NDIM + 1)/2)
C
      WRITE (LUPRI, '(//,2A)') '  ', HEAD
      WRITE (LUPRI, '(80A)')   '  ', ('-', I = 1, LEN(HEAD))
      WRITE (LUPRI, '()')
      CALL OUTPAK(OUTMAT,NDIM,1,LUPRI)
C     IOFF = 0
C     DO 100 I = 1, NDIM
C        WRITE (LUPRI, '(6F12.6)') (OUTMAT(IOFF + J), J = 1, I)
C        WRITE (LUPRI, '()')
C        IOFF = IOFF + I
C 100 CONTINUE
      RETURN
      END
