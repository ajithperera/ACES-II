
      SUBROUTINE AMPT11(ICORE,MAXCOR,IUHF)
C
C THIS ROUTINE DRIVES A NUMBER OF ROUTINES WHICH DETERMINE AND PRINT
C  OUT THE LARGEST T2 AMPLITUDES.  ALL THREE SPIN CASES DONE FOR UHF,
C  AB ONLY FOR RHF.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ICORE(MAXCOR)
      COMMON /FLAGS/ IFLAGS(100)
      IF(IUHF.NE.0)THEN
       DO 10 I=3,4
        CALL SORTT1(ICORE,MAXCOR,I,90,IFLAGS(14),'T')
10     CONTINUE
      ELSE
       CALL SORTT1(ICORE,MAXCOR,3,90,IFLAGS(14),'T')
      ENDIF
      RETURN
      END
