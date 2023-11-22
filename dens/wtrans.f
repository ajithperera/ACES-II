      SUBROUTINE WTRANS(XIA,WMATOA,WMATVA,WMATOB,WMATVB,XSCR,IDIR)
C
C  NOTES:  I need to make sure I have both a forward and backward
C          transformation in this routine.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER DIRPRD,POP,VRT
      DIMENSION XIA(1),XSCR(1),WMATOA(1),WMATVA(1),WMATOB(1),WMATVB(1)
C
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON /SYM2/ POP(8,2),VRT(8,2),NJUNK(6)
C
      DATA ONE /1.0/
      DATA ZILCH /0.0/
C
C  Forward transformation (canonical MO --> semicanonical MO)
C
      IF(IDIR.EQ.0) THEN
C
C  First do the alpha spin case.
C
        IOFFX=1
        IOFFO=1
        IOFFV=1
        DO 100 IRREP=1,NIRREP
          NOCC=POP(IRREP,1)
          NVRT=VRT(IRREP,1)
          IF(MIN(NOCC,NVRT).NE.0) THEN
            CALL XGEMM('T','N',NVRT,NOCC,NVRT,ONE,WMATVA(IOFFV),NVRT,
     &                 XIA(IOFFX),NVRT,ZILCH,XSCR,NVRT)
            CALL XGEMM('N','N',NVRT,NOCC,NOCC,ONE,XSCR,NVRT,
     &                 WMATOA(IOFFO),NOCC,ZILCH,XIA(IOFFX),NVRT)
          ENDIF
C
          IOFFX=IOFFX+NOCC*NVRT
          IOFFO=IOFFO+NOCC*NOCC
          IOFFV=IOFFV+NVRT*NVRT
  100   CONTINUE
C
C  Now do the beta spin case.
C
        IOFFO=1
        IOFFV=1
        DO 110 IRREP=1,NIRREP
          NOCC=POP(IRREP,2)
          NVRT=VRT(IRREP,2)
          IF(MIN(NOCC,NVRT).NE.0) THEN
            CALL XGEMM('T','N',NVRT,NOCC,NVRT,ONE,WMATVB(IOFFV),NVRT,
     &                 XIA(IOFFX),NVRT,ZILCH,XSCR,NVRT)
            CALL XGEMM('N','N',NVRT,NOCC,NOCC,ONE,XSCR,NVRT,
     &                 WMATOB(IOFFO),NOCC,ZILCH,XIA(IOFFX),NVRT)
          ENDIF
C
          IOFFX=IOFFX+NOCC*NVRT
          IOFFO=IOFFO+NOCC*NOCC
          IOFFV=IOFFV+NVRT*NVRT
  110   CONTINUE
C
C  Reverse transformation (semicanonical MO --> canonical MO)
C
      ELSEIF(IDIR.EQ.1) THEN
C
C  First do the alpha spin case.
C
        IOFFX=1
        IOFFO=1
        IOFFV=1
        DO 200 IRREP=1,NIRREP
          NOCC=POP(IRREP,1)
          NVRT=VRT(IRREP,1)
          IF(MIN(NOCC,NVRT).NE.0) THEN
            CALL XGEMM('N','N',NVRT,NOCC,NVRT,ONE,WMATVA(IOFFV),NVRT,
     &                 XIA(IOFFX),NVRT,ZILCH,XSCR,NVRT)
            CALL XGEMM('N','T',NVRT,NOCC,NOCC,ONE,XSCR,NVRT,
     &                 WMATOA(IOFFO),NOCC,ZILCH,XIA(IOFFX),NVRT)
          ENDIF
C
          IOFFX=IOFFX+NOCC*NVRT
          IOFFO=IOFFO+NOCC*NOCC
          IOFFV=IOFFV+NVRT*NVRT
  200   CONTINUE
C
C  Now do the beta spin case.
C
        IOFFO=1
        IOFFV=1
        DO 210 IRREP=1,NIRREP
          NOCC=POP(IRREP,2)
          NVRT=VRT(IRREP,2)
          IF(MIN(NOCC,NVRT).NE.0) THEN
            CALL XGEMM('N','N',NVRT,NOCC,NVRT,ONE,WMATVB(IOFFV),NVRT,
     &                 XIA(IOFFX),NVRT,ZILCH,XSCR,NVRT)
            CALL XGEMM('N','T',NVRT,NOCC,NOCC,ONE,XSCR,NVRT,
     &                 WMATOB(IOFFO),NOCC,ZILCH,XIA(IOFFX),NVRT)
          ENDIF
C
          IOFFX=IOFFX+NOCC*NVRT
          IOFFO=IOFFO+NOCC*NOCC
          IOFFV=IOFFV+NVRT*NVRT
  210   CONTINUE
C
      ELSE
        WRITE(6,9000)
 9000   FORMAT(T3,'@WTRANS-F, Illegal selection for transformation!')
        CALL ERREX
      ENDIF
C
C  All done!
C
      RETURN
      END
