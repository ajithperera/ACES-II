      SUBROUTINE SEMCPHF(XIA,WMATOA,WMATVA,WMATOB,WMATVB,
     &                   IRREPX,XSCR,IDIR)
C
C TRANSFORMS TWO-DIMENSIONAL MATRICES OF IRREP IRREPX FROM CANONICAL
C TO SEMI-CANONICAL REPRESENTATION AND VICE VERS
C
CEND
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*1 MOD1,MOD2
C
      INTEGER DIRPRD,POP,VRT
      DIMENSION XIA(1),XSCR(1),WMATOA(1),WMATVA(1),WMATOB(1),WMATVB(1)
C
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NJUNK(6)
C
      DATA ONE /1.0D0/
      DATA ZILCH /0.0D0/
C
      IF(IDIR.EQ.0) THEN
C
C  FORWARD TRANSFORMATION (CANONICAL MO --> SEMICANONICAL MO)
C
       MOD1='T'
       MOD2='N'
C
      ELSE IF(IDIR.EQ.1) THEN 
C
C  REVERSE TRANSFORMATION (CANONICAL MO --> SEMICANONICAL MO)
C
       MOD1='N'
       MOD2='T'
C
      ELSE
C
       WRITE(6,9000)
9000   FORMAT(T3,'@WTRANS-F, Illegal selection for transformation!')
       CALL ERREX
C
      ENDIF
C
C  FIRST DO ALPHA SPIN CASE
C
      IOFFX=1
      IOFFO=1
      DO 100 IRREP1=1,NIRREP
       IRREP2=DIRPRD(IRREPX,IRREP1)
       NOCC=POP(IRREP1,1)
       NVRT=VRT(IRREP2,1)
       IOFFV=0
       DO 109 IRR=1,IRREP2-1
        IOFFV=IOFFV+VRT(IRR,1)*VRT(IRR,1)
109     CONTINUE
        IF(MIN(NOCC,NVRT).NE.0) THEN
         CALL XGEMM(MOD1,'N',NVRT,NOCC,NVRT,ONE,WMATVA(IOFFV),NVRT,
     &              XIA(IOFFX),NVRT,ZILCH,XSCR,NVRT)
         CALL XGEMM('N',MOD2,NVRT,NOCC,NOCC,ONE,XSCR,NVRT,
     &              WMATOA(IOFFO),NOCC,ZILCH,XIA(IOFFX),NVRT)
        ENDIF
C
        IOFFX=IOFFX+NOCC*NVRT
        IOFFO=IOFFO+NOCC*NOCC
100    CONTINUE
C
C  NOW DO BETA SPIN CASE
C
       IOFFO=1
       DO 110 IRREP1=1,NIRREP
        IRREP2=DIRPRD(IRREPX,IRREP1)
        NOCC=POP(IRREP1,2)
        NVRT=VRT(IRREP2,2)
        IOFFV=0
        DO 119 IRR=1,IRREP2-1
         IOFFV=IOFFV+VRT(IRR,2)*VRT(IRR,2)
119     CONTINUE
        IF(MIN(NOCC,NVRT).NE.0) THEN
         CALL XGEMM(MOD1,'N',NVRT,NOCC,NVRT,ONE,WMATVB(IOFFV),NVRT,
     &              XIA(IOFFX),NVRT,ZILCH,XSCR,NVRT)
         CALL XGEMM('N',MOD2,NVRT,NOCC,NOCC,ONE,XSCR,NVRT,
     &              WMATOB(IOFFO),NOCC,ZILCH,XIA(IOFFX),NVRT)
        ENDIF
C
        IOFFX=IOFFX+NOCC*NVRT
        IOFFO=IOFFO+NOCC*NOCC
110    CONTINUE
C
C  ALL DONE, RETURN
C
      RETURN
      END
