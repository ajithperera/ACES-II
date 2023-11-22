      SUBROUTINE WOO(XOO,WMATOA,WMATOB,XSCR,IA,IB,IDIR)
C
C  NOTES:  I need to make sure I have both a forward and backward
C          transformation in this routine.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER DIRPRD,POP,VRT
      DIMENSION XOO(1),XSCR(1),WMATOA(1),WMATOB(1)
C
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON /SYM2/ POP(8,2),VRT(8,2),NJUNK(6)
C
      DATA ONE /1.0/
      DATA ZILCH /0.0/
C
C  Forward transformation (canonical MO --> semicanonical MO)
C
      IF(IDIR.EQ.0) THEN
        IF(IA.EQ.0) GOTO 101
C
C  First do the alpha spin case.
C
        IOFFX=1
        IOFFO=1
        DO 100 IRREP=1,NIRREP
          NOCC=POP(IRREP,1)
          IF(NOCC.NE.0) THEN
            CALL XGEMM('T','N',NOCC,NOCC,NOCC,ONE,WMATOA(IOFFO),NOCC,
     &                 XOO(IOFFX),NOCC,ZILCH,XSCR,NOCC)
            CALL XGEMM('N','N',NOCC,NOCC,NOCC,ONE,XSCR,NOCC,
     &                 WMATOA(IOFFO),NOCC,ZILCH,XOO(IOFFX),NOCC)
          ENDIF
C
          IOFFX=IOFFX+NOCC*NOCC
          IOFFO=IOFFO+NOCC*NOCC
  100   CONTINUE
C
C  Now do the beta spin case.
C
  101   IF(IB.EQ.0) GOTO 111
        IOFFO=1
        DO 110 IRREP=1,NIRREP
          NOCC=POP(IRREP,2)
          IF(NOCC.NE.0) THEN
            CALL XGEMM('T','N',NOCC,NOCC,NOCC,ONE,WMATOB(IOFFO),NOCC,
     &                 XOO(IOFFX),NOCC,ZILCH,XSCR,NOCC)
            CALL XGEMM('N','N',NOCC,NOCC,NOCC,ONE,XSCR,NOCC,
     &                 WMATOB(IOFFO),NOCC,ZILCH,XOO(IOFFX),NOCC)
          ENDIF
C
          IOFFX=IOFFX+NOCC*NOCC
          IOFFO=IOFFO+NOCC*NOCC
  110   CONTINUE
  111   CONTINUE
C
C  Reverse transformation (semicanonical MO --> canonical MO)
C
      ELSEIF(IDIR.EQ.1) THEN
        IF(IA.EQ.0) GOTO 201
C
C  First do the alpha spin case.
C
        IOFFX=1
        IOFFO=1
        DO 200 IRREP=1,NIRREP
          NOCC=POP(IRREP,1)
          IF(NOCC.NE.0) THEN
            CALL XGEMM('N','N',NOCC,NOCC,NOCC,ONE,WMATOA(IOFFO),NOCC,
     &                 XOO(IOFFX),NOCC,ZILCH,XSCR,NOCC)
            CALL XGEMM('N','T',NOCC,NOCC,NOCC,ONE,XSCR,NOCC,
     &                 WMATOA(IOFFO),NOCC,ZILCH,XOO(IOFFX),NOCC)
          ENDIF
C
          IOFFX=IOFFX+NOCC*NOCC
          IOFFO=IOFFO+NOCC*NOCC
  200   CONTINUE
C
C  Now do the beta spin case.
C
  201   IF(IB.EQ.0) GOTO 211
        IOFFO=1
        DO 210 IRREP=1,NIRREP
          NOCC=POP(IRREP,2)
          IF(NOCC.NE.0) THEN
            CALL XGEMM('N','N',NOCC,NOCC,NOCC,ONE,WMATOB(IOFFO),NOCC,
     &                 XOO(IOFFX),NOCC,ZILCH,XSCR,NOCC)
            CALL XGEMM('N','T',NOCC,NOCC,NOCC,ONE,XSCR,NOCC,
     &                 WMATOB(IOFFO),NOCC,ZILCH,XOO(IOFFX),NOCC)
          ENDIF
C
          IOFFX=IOFFX+NOCC*NOCC
          IOFFO=IOFFO+NOCC*NOCC
  210   CONTINUE
  211   CONTINUE
C
      ELSE
        WRITE(LUOUT,9000)
 9000   FORMAT(T3,'@WOO-F, Illegal selection for transformation!')
        CALL ERREX
      ENDIF
C
C  All done!
C
      RETURN
      END
