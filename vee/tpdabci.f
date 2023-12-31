      SUBROUTINE TPDABCI(ICORE,MAXCOR,IUHF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
      DIMENSION ICORE(MAXCOR)
      COMMON/REFTYPE/MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
C
C DRIVER PROGRAM FOR THE G(IJ,KA) TERMS)
C
CEND
C
C CODED SEPTEMBER/93 JG
C
      CALL TPDABCI1(ICORE,MAXCOR,IUHF)
C
      IF (CC) CALL TPDABCI2(ICORE,MAXCOR,IUHF)
C
      CALL TPDABCI3(ICORE,MAXCOR,IUHF)
C
      IF (CC .AND. .NOT. CC2) CALL TPDABCI4(ICORE,MAXCOR,IUHF)
C
      IF (.NOT. CC2) CALL TPDABCI5(ICORE,MAXCOR,IUHF)
C
      CALL R2L2RNGD(ICORE,MAXCOR,IUHF,0.0D0,.TRUE.)
C
      IF (CC .AND. .NOT. CC2) CALL TPDABCI6(ICORE,MAXCOR,IUHF)
C
      CALL T2L2RNGD(ICORE,MAXCOR,IUHF,0.0D0)
C
      IF (.NOT. CC2) CALL TPDABCI7(ICORE,MAXCOR,IUHF)
C
      IF (CC) CALL TPDABCI8(ICORE,MAXCOR,IUHF)
C
C ALL DONE, RETURN
C
      RETURN
      END
