      SUBROUTINE TSPABCI2(CORE,MAXCOR,LIST1,LIST2,ICODE,IRREPX)
      IMPLICIT NONE
C-----------------------------------------------------------------------
      DOUBLE PRECISION CORE
      INTEGER MAXCOR,LIST1,LIST2,ICODE,IRREPX
C-----------------------------------------------------------------------
      INTEGER IINTLN,IFLTLN,IINTFP,IALONE,IBITWD,
     &        NSTART,NIRREP,IRREPA,IRREPB,DIRPRD,IRPDPD,ISYTYP,ID,
     &        POP,VRT,NTAA,NTBB,NF1AA,NF2AA,NF1BB,NF2BB,IFLAGS,IPRINT
C-----------------------------------------------------------------------
      INTEGER IRPAB,IRPCI,I000,I010,I020,I030,I040,NEED,IDIS
      LOGICAL INCORE
C-----------------------------------------------------------------------
      DIMENSION CORE(1)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     1                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
      COMMON /FLAGS/  IFLAGS(100)
      EQUIVALENCE(IPRINT,IFLAGS( 1))

      Print*, "I am here"
C
C     Subroutine to switch AbCi/AbcI to bACi/bAcI or vice versa.
C
      IF(ICODE.EQ.1)THEN
      WRITE(6,1000)
 1000 FORMAT(' @TSPABCI2-I, (Alpha,Beta) ---> (Beta,Alpha). ')
      ENDIF
C
      IF(ICODE.EQ.2)THEN
      WRITE(6,1010)
 1010 FORMAT(' @TSPABCI2-I, (Beta,Alpha) ---> (Alpha,Beta). ')
      ENDIF
C
      IF(ICODE.NE.1.AND.ICODE.NE.2)THEN
      WRITE(6,1020) ICODE
 1020 FORMAT(' @TSPABCI2-I, Invalid value of ICODE. ICODE is ',I10)
      STOP
      ENDIF
C
      INCORE = .TRUE.
C
      DO   20 IRPCI = 1,NIRREP
      IRPAB = DIRPRD(IRREPX,IRPCI)
      I000 = 1
      I010 = I000 + MAX0(IRPDPD(IRPCI,12),IRPDPD(IRPAB,13))
      I020 = I010 + MAX0(IRPDPD(IRPCI,12),IRPDPD(IRPAB,13))
      I030 = I020 + MAX0(IRPDPD(IRPCI,12),IRPDPD(IRPAB,13))
      I040 = I030 +      IRPDPD(IRPCI,12)*IRPDPD(IRPAB,13)
      NEED = IINTFP * I040
C
      IF(NEED.LE.MAXCOR.AND.INCORE)THEN
C
        IF(IPRINT.GT.10) WRITE(6,1040)
C
        CALL GETLST(CORE(I030),1,IRPDPD(IRPCI,12),2,IRPCI,LIST1)
        IF(ICODE.EQ.1)THEN
          CALL SYMTR3(IRPAB,VRT(1,1),VRT(1,2),
     &                IRPDPD(IRPAB,13),IRPDPD(IRPCI,12),
     &                CORE(I030),CORE(I000),CORE(I010),CORE(I020))
        ELSE
          CALL SYMTR3(IRPAB,VRT(1,2),VRT(1,1),
     &                IRPDPD(IRPAB,13),IRPDPD(IRPCI,12),
     &                CORE(I030),CORE(I000),CORE(I010),CORE(I020))
        ENDIF
        CALL PUTLST(CORE(I030),1,IRPDPD(IRPCI,12),2,IRPCI,LIST1)
C
      ELSE
C
C     Handle one distribution at a time.
C
        WRITE(6,1050)
C
        DO 10 IDIS=1,IRPDPD(IRPCI,12)
C
        CALL GETLST(CORE(I030),IDIS,1,2,IRPCI,LIST1)
        IF(ICODE.EQ.1)THEN
          CALL SYMTR3(IRPAB,VRT(1,1),VRT(1,2),IRPDPD(IRPAB,13),1,
     &                CORE(I030),CORE(I000),CORE(I010),CORE(I020))
        ELSE
          CALL SYMTR3(IRPAB,VRT(1,2),VRT(1,1),IRPDPD(IRPAB,13),1,
     &                CORE(I030),CORE(I000),CORE(I010),CORE(I020))
        ENDIF
        CALL PUTLST(CORE(I030),IDIS,1,2,IRPCI,LIST1)
   10 CONTINUE
      ENDIF
   20 CONTINUE
C
      DO  120 IRPCI = 1,NIRREP
      IRPAB = DIRPRD(IRREPX,IRPCI)
      I000 = 1
      I010 = I000 + MAX0(IRPDPD(IRPCI,11),IRPDPD(IRPAB,13))
      I020 = I010 + MAX0(IRPDPD(IRPCI,11),IRPDPD(IRPAB,13))
      I030 = I020 + MAX0(IRPDPD(IRPCI,11),IRPDPD(IRPAB,13))
      I040 = I030 +      IRPDPD(IRPCI,11)*IRPDPD(IRPAB,13)
      NEED = IINTFP * I040
C
      IF(NEED.LE.MAXCOR.AND.INCORE)THEN
C
        IF(IPRINT.GT.10) WRITE(6,1040)
C
        CALL GETLST(CORE(I030),1,IRPDPD(IRPCI,11),2,IRPCI,LIST2)
        IF(ICODE.EQ.1)THEN
          CALL SYMTR3(IRPAB,VRT(1,1),VRT(1,2),
     &                IRPDPD(IRPAB,13),IRPDPD(IRPCI,11),
     &                CORE(I030),CORE(I000),CORE(I010),CORE(I020))
        ELSE
          CALL SYMTR3(IRPAB,VRT(1,2),VRT(1,1),
     &                IRPDPD(IRPAB,13),IRPDPD(IRPCI,11),
     &                CORE(I030),CORE(I000),CORE(I010),CORE(I020))
        ENDIF
        CALL PUTLST(CORE(I030),1,IRPDPD(IRPCI,11),2,IRPCI,LIST2)
C
      ELSE
C
C     Handle one distribution at a time.
C
        WRITE(6,1050)
C
        DO 110 IDIS=1,IRPDPD(IRPCI,11)
C
        CALL GETLST(CORE(I030),IDIS,1,2,IRPCI,LIST2)
        IF(ICODE.EQ.1)THEN
          CALL SYMTR3(IRPAB,VRT(1,1),VRT(1,2),IRPDPD(IRPAB,13),1,
     &                CORE(I030),CORE(I000),CORE(I010),CORE(I020))
        ELSE
          CALL SYMTR3(IRPAB,VRT(1,2),VRT(1,1),IRPDPD(IRPAB,13),1,
     &                CORE(I030),CORE(I000),CORE(I010),CORE(I020))
        ENDIF
        CALL PUTLST(CORE(I030),IDIS,1,2,IRPCI,LIST2)
  110   CONTINUE
      ENDIF
  120 CONTINUE
      RETURN
 1040 FORMAT(' @TSPABCI2-I, In-core algorithm is being used. ')
 1050 FORMAT(' @TSPABCI2-I, Out-of-core algorithm is being used. ')
      END
