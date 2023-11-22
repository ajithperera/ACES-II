      SUBROUTINE TSPABCI(CORE,MAXCOR,LIST,ICODE,IRREPX)
      IMPLICIT NONE
C-----------------------------------------------------------------------
      DOUBLE PRECISION CORE
      INTEGER MAXCOR,LIST,ICODE,IRREPX
C-----------------------------------------------------------------------
      INTEGER IINTLN,IFLTLN,IINTFP,IALONE,IBITWD,
     &        NSTART,NIRREP,IRREPA,IRREPB,DIRPRD,
     &        IRPDPD,ISYTYP,ID,
     &        POP,VRT,NTAA,NTBB,NF1AA,NF2AA,NF1BB,NF2BB,
     &        IFLAGS,IPRINT,IOFFVV,IOFFOO,IOFFVO
C-----------------------------------------------------------------------
      INTEGER IOFFOV,IRPAB,IRPCI,I000,I010,I020,I030,I040,NEED
      LOGICAL INCORE
C-----------------------------------------------------------------------
      DIMENSION CORE(1)
      DIMENSION IOFFOV(8,8,4)
C-----------------------------------------------------------------------
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     1                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
      COMMON /FLAGS/  IFLAGS(100)
      EQUIVALENCE(IPRINT,IFLAGS( 1))
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
C
      WRITE(6,1000)
 1000 FORMAT(' @TSPABCI-I, AbIc/AbcI record transposition. ')
C
      IF(ICODE.EQ.1)THEN
      WRITE(6,1010)
 1010 FORMAT(' @TSPABCI-I, (O,v) ---> (v,O) transposition. ')
      ENDIF
      IF(ICODE.EQ.2)THEN
      WRITE(6,1020)
 1020 FORMAT(' @TSPABCI-I, (v,O) ---> (O,v) transposition. ')
      ENDIF
      IF(ICODE.NE.1.AND.ICODE.NE.2)THEN
      WRITE(6,1030) ICODE
 1030 FORMAT(' @TSPABCI-I, Invalid value of ICODE. ICODE is ',I10)
      STOP
      ENDIF
C
      INCORE = .TRUE.
      CALL MKOFOV(IOFFOV)
C
      DO   10 IRPCI = 1,NIRREP
      IRPAB = DIRPRD(IRREPX,IRPCI)
      I000 = 1
      I010 = I000 + MAX0(IRPDPD(IRPAB,13),IRPDPD(IRPCI,12))
      I020 = I010 + MAX0(IRPDPD(IRPAB,13),IRPDPD(IRPCI,12))
      I030 = I020 + MAX0(IRPDPD(IRPAB,13),IRPDPD(IRPCI,12))
      I040 = I030 +      IRPDPD(IRPAB,13)*IRPDPD(IRPCI,12)
      NEED = I040 * IINTFP
C
      IF(NEED.LE.MAXCOR.AND.INCORE)THEN
C
        IF(IPRINT.GT.10) WRITE(6,1040)
C
        CALL GETLST(CORE(I030),1,IRPDPD(IRPCI,12),2,IRPCI,LIST)
C
        IF(ICODE.EQ.1)THEN
          CALL SYMTR1(IRPCI,POP(1,1),VRT(1,2),IRPDPD(IRPAB,13),
     &                CORE(I030),CORE(I000),CORE(I010),CORE(I020))
        ELSE
          CALL SYMTR1(IRPCI,VRT(1,2),POP(1,1),IRPDPD(IRPAB,13),
     &              CORE(I030),CORE(I000),CORE(I010),CORE(I020))
        ENDIF
C
        CALL PUTLST(CORE(I030),1,IRPDPD(IRPCI,12),2,IRPCI,LIST)
C
      ELSE
C
C     Otherwise we do the transposition out-of-core.
C
      WRITE(6,1050)
      if(irrepx.ne.1)then
        write(6,*) ' @tspabci-f, ocsymtr1 not coded for this case '
        call errex
      endif
C
        CALL OCSYMTR1(CORE,IRPDPD(IRPCI,13),IRPDPD(IRPCI,12),IRPCI,
     &                LIST,IOFFOV,ICODE)
      ENDIF
   10 CONTINUE
      RETURN
 1040 FORMAT(' @TSPABCI-I, In-core algorithm is being used. ')
 1050 FORMAT(' @TSPABCI-I, Out-of-core algorithm is being used. ')
      END
