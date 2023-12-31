      SUBROUTINE DRE3EN(ICORE,MAXCOR,IUHF,IDOPPL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,QCISD,UCC
      LOGICAL ROHF4,ITRFLG
C
C THIRD-ORDER MBPT ENERGY DRIVER.
C
      DIMENSION ICORE(MAXCOR)
C
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /SWITCH/MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &               QCISD,UCC
      COMMON /ROHF/ ROHF4,ITRFLG
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW     
C
      EQUIVALENCE(METHOD,IFLAGS(2))
C
      IF(IUHF.EQ.1)IBOT=1
      IF(IUHF.EQ.0)IBOT=3
      ITOP=6
      IF(IDOPPL.NE.0)ITOP=1
      CALL TIMER(1)
C
      DO 10 ITYPE=1,ITOP,5
        IF(ITYPE.EQ.6.AND.IFLAGS(93).EQ.2)THEN
         CALL DRAOLAD(ICORE,MAXCOR,IUHF,.FALSE.,1,1,43,60,243,260)
        ELSE
         CALL DRLAD(ICORE,MAXCOR,IUHF,ITYPE)
        ENDIF
10    CONTINUE
C
      IF(METHOD.GT.9.AND.SING1.AND.(.NOT.QCISD).OR.
     &               (ROHF4.AND.ITRFLG)) THEN
        CALL T12INT2(ICORE,MAXCOR,IUHF)
      ENDIF
C
      DO 15 ISPIN=IBOT,3
       CALL DRRNG(ICORE,MAXCOR,ISPIN,ITYPE,IUHF)
15    CONTINUE
C
      CALL TIMER(1)
C
      IF(IFLAGS(1).GE.10)THEN
       WRITE(LUOUT,100)TIMENEW
100    FORMAT(T3,'@DRE3EN-I, T2-W contractions required ',F9.3,
     &           ' seconds.')
      ENDIF
      RETURN
      END
