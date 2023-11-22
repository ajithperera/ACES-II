      SUBROUTINE MOLLAB(A,LU,LUERR)
C
C  16-Jun-1986 hjaaj
C  (as SEARCH but CHARACTER*8 instead of REAL*8 labels)
C
C  Purpose:
C     Search for MOLECULE labels on file LU
C
      CHARACTER*8 A, B(4), C
      DATA C/'********'/
    1 READ (LU,END=3,ERR=6) B
      IF (B(1).NE.C) GO TO 1
      IF (B(4).NE.A) GO TO 1
      IF (LUERR.LT.0) LUERR = 0
      RETURN
C
    3 IF (LUERR.LT.0) THEN
         LUERR = -1
         RETURN
      ELSE
         WRITE(LUERR,4)A,LU
         CALL TRACE
         STOP
      END IF
    4 FORMAT(/' *** ERROR (MOLLAB), MOLECULE label ',A8,
     *        ' not found on unit',I4)
C
    6 IF (LUERR.LT.0) THEN
    8    LUERR = -2
         RETURN
      ELSE
         WRITE (LUERR,7) LU,A
         CALL TRACE
         STOP 
      END IF
    7 FORMAT(/' *** ERROR (MOLLAB), error reading unit',I4,
     *       /T22,'when searching for label ',A8)
      END
