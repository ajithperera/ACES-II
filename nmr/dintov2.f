      SUBROUTINE DINTOV2(DIOV,ICORE,MAXCOR,IUHF,IP,ANTI)
C
CEND
C
C CODED SEPT/91 JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,ANTI 
C
      DIMENSION ICORE(MAXCOR),DIOV(1)
C
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW 
C
      CALL TIMER(1)
C
C CALCULATE FIRST CONTRIBUTIONS TO DI(I,A)
C
      CALL DIOV21(DIOV,ICORE,MAXCOR,IUHF,IP,ANTI)
C
C FOR HIGHER ORDER PERTURBATION THEORY AND CC METHODS, ADDITIONAL TERMS
C HAVE TO BE CONSIDERED
C
      IF(.NOT.MBPT2) THEN
C
c       CALL DIOV22(DIOV,ICORE,MAXCOR,IUHF,IP,ANTI)
c       CALL DIOV23(DIOV,ICORE,MAXCOR,IUHF,IP,ANTI)
C
C ADDITIONAL TERMS FOR FOURTH ORDER AND QCISD AND CCSD
C
       IF(.NOT.MBPT3.AND..NOT.CCD) THEN
C
c        CALL DIOV24(DIOV,ICORE,MAXCOR,IUHF,IP,ANTI)
c        CALL DIOV25(DIOV,ICORE,MAXCOR,IUHF,IP,ANTI)
c        CALL DIOV26(DIOV,ICORE,MAXCOR,IUHF,IP,ANTI)
C
       ENDIF
      ENDIF
C
C ALL DONE
C
      CALL TIMER(1)
      write(6,6001) TIMENEW
6001  FORMAT(' Calculation of the contributions of <pq||rs> to',
     &       ' d I(i,a)/ d chi',/,' required ',f5.1,' seconds.')
      RETURN
C
      END 