      SUBROUTINE TSPABCI2(CORE,MAXCOR,ICODE,EXTRNL_CCSD)
      IMPLICIT INTEGER (A-Z)
      LOGICAL INCORE
      LOGICAL EXTRNL_CCSD
      DOUBLE PRECISION CORE(1)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     1                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
      COMMON /FLAGS/  IFLAGS(100)
      EQUIVALENCE(IPRINT,IFLAGS( 1))
      EQUIVALENCE(IDRLVL,IFLAGS( 3))
C
      COMMON /GAMLIS/ LISGO1,LISGO2,LISGO3,LISGO4,LISGV1,LISGV2,
     1                LISGV3,LISGV4
C
      COMMON /LISWI/  LWIC11,LWIC12,LWIC13,LWIC14,
     1                LWIC15,LWIC16,LWIC17,LWIC18,
     1                LWIC21,LWIC22,LWIC23,LWIC24,
     1                LWIC25,LWIC26,LWIC27,LWIC28,
     1                LWIC31,LWIC32,LWIC33,
     1                LWIC34,LWIC35,LWIC36,
     1                LWIC37,LWIC38,LWIC39,LWIC40,LWIC41,LWIC42

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
      I000 = 1
      I010 = I000 + MAX0(IRPDPD(IRPCI,12),IRPDPD(IRPCI,13))
      I020 = I010 + MAX0(IRPDPD(IRPCI,12),IRPDPD(IRPCI,13))
      I030 = I020 + MAX0(IRPDPD(IRPCI,12),IRPDPD(IRPCI,13))
      I040 = I030 +      IRPDPD(IRPCI,12)*IRPDPD(IRPCI,13)
      NEED = IINTFP * I040
C
      IF(NEED.LE.MAXCOR.AND.INCORE)THEN
C
      IF(IPRINT.GT.10) WRITE(6,1040)
 1040 FORMAT(' @TSPABCI2-I, In-core algorithm is being used. ')
C
      CALL GETLST(CORE(I030),1,IRPDPD(IRPCI,12),2,IRPCI,LWIC17)
      IF(ICODE.EQ.1)THEN
      CALL SYMTR3(IRPCI,VRT(1,1),VRT(1,2),
     1            IRPDPD(IRPCI,13),IRPDPD(IRPCI,12),
     1            CORE(I030),CORE(I000),CORE(I010),CORE(I020))
      ELSE
      CALL SYMTR3(IRPCI,VRT(1,2),VRT(1,1),
     1            IRPDPD(IRPCI,13),IRPDPD(IRPCI,12),
     1            CORE(I030),CORE(I000),CORE(I010),CORE(I020))
      ENDIF
      CALL PUTLST(CORE(I030),1,IRPDPD(IRPCI,12),2,IRPCI,LWIC17)
C
      IF(LWIC17.NE.LWIC27)THEN
      CALL GETLST(CORE(I030),1,IRPDPD(IRPCI,12),2,IRPCI,LWIC27)
      IF(ICODE.EQ.1)THEN
      CALL SYMTR3(IRPCI,VRT(1,1),VRT(1,2),
     1            IRPDPD(IRPCI,13),IRPDPD(IRPCI,12),
     1            CORE(I030),CORE(I000),CORE(I010),CORE(I020))
      ELSE
      CALL SYMTR3(IRPCI,VRT(1,2),VRT(1,1),
     1            IRPDPD(IRPCI,13),IRPDPD(IRPCI,12),
     1            CORE(I030),CORE(I000),CORE(I010),CORE(I020))
      ENDIF
      CALL PUTLST(CORE(I030),1,IRPDPD(IRPCI,12),2,IRPCI,LWIC27)
      ENDIF
C
      IF(IDRLVL.GT.0.OR.EXTRNL_CCSD)THEN
      CALL GETLST(CORE(I030),1,IRPDPD(IRPCI,12),2,IRPCI,LISGV3)
      IF(ICODE.EQ.1)THEN
      CALL SYMTR3(IRPCI,VRT(1,1),VRT(1,2),
     1            IRPDPD(IRPCI,13),IRPDPD(IRPCI,12),
     1            CORE(I030),CORE(I000),CORE(I010),CORE(I020))
      ELSE
      CALL SYMTR3(IRPCI,VRT(1,2),VRT(1,1),
     1            IRPDPD(IRPCI,13),IRPDPD(IRPCI,12),
     1            CORE(I030),CORE(I000),CORE(I010),CORE(I020))
      ENDIF
      CALL PUTLST(CORE(I030),1,IRPDPD(IRPCI,12),2,IRPCI,LISGV3)
      ENDIF
C
      ELSE
C
C     Handle one distribution at a time.
C
      WRITE(6,1050)
 1050 FORMAT(' @TSPABCI2-I, Out-of-core algorithm is being used. ')
C
      DO 10 IDIS=1,IRPDPD(IRPCI,12)
C
      CALL GETLST(CORE(I030),IDIS,1,2,IRPCI,LWIC17)
      IF(ICODE.EQ.1)THEN
      CALL SYMTR3(IRPCI,VRT(1,1),VRT(1,2),IRPDPD(IRPCI,13),1,
     1            CORE(I030),CORE(I000),CORE(I010),CORE(I020))
      ELSE
      CALL SYMTR3(IRPCI,VRT(1,2),VRT(1,1),IRPDPD(IRPCI,13),1,
     1            CORE(I030),CORE(I000),CORE(I010),CORE(I020))
      ENDIF
      CALL PUTLST(CORE(I030),IDIS,1,2,IRPCI,LWIC17)
C
      IF(LWIC17.NE.LWIC27)THEN
      CALL GETLST(CORE(I030),IDIS,1,2,IRPCI,LWIC27)
      IF(ICODE.EQ.1)THEN
      CALL SYMTR3(IRPCI,VRT(1,1),VRT(1,2),IRPDPD(IRPCI,13),1,
     1            CORE(I030),CORE(I000),CORE(I010),CORE(I020))
      ELSE
      CALL SYMTR3(IRPCI,VRT(1,2),VRT(1,1),IRPDPD(IRPCI,13),1,
     1            CORE(I030),CORE(I000),CORE(I010),CORE(I020))
      ENDIF
      CALL PUTLST(CORE(I030),IDIS,1,2,IRPCI,LWIC27)
      ENDIF
C
      IF(IDRLVL.GT.0.OR.EXTRNL_CCSD)THEN
      CALL GETLST(CORE(I030),IDIS,1,2,IRPCI,LISGV3)
      IF(ICODE.EQ.1)THEN
      CALL SYMTR3(IRPCI,VRT(1,1),VRT(1,2),IRPDPD(IRPCI,13),1,
     1            CORE(I030),CORE(I000),CORE(I010),CORE(I020))
      ELSE
      CALL SYMTR3(IRPCI,VRT(1,2),VRT(1,1),IRPDPD(IRPCI,13),1,
     1            CORE(I030),CORE(I000),CORE(I010),CORE(I020))
      ENDIF
      CALL PUTLST(CORE(I030),IDIS,1,2,IRPCI,LISGV3)
      ENDIF
C
   10 CONTINUE
C
      ENDIF
C
   20 CONTINUE
C
      DO  120 IRPCI = 1,NIRREP
      I000 = 1
      I010 = I000 + MAX0(IRPDPD(IRPCI,11),IRPDPD(IRPCI,13))
      I020 = I010 + MAX0(IRPDPD(IRPCI,11),IRPDPD(IRPCI,13))
      I030 = I020 + MAX0(IRPDPD(IRPCI,11),IRPDPD(IRPCI,13))
      I040 = I030 +      IRPDPD(IRPCI,11)*IRPDPD(IRPCI,13)
      NEED = IINTFP * I040
C
      IF(NEED.LE.MAXCOR.AND.INCORE)THEN
C
      IF(IPRINT.GT.10) WRITE(6,1040)
C
      CALL GETLST(CORE(I030),1,IRPDPD(IRPCI,11),2,IRPCI,LWIC18)
      IF(ICODE.EQ.1)THEN
      CALL SYMTR3(IRPCI,VRT(1,1),VRT(1,2),
     1            IRPDPD(IRPCI,13),IRPDPD(IRPCI,11),
     1            CORE(I030),CORE(I000),CORE(I010),CORE(I020))
      ELSE
      CALL SYMTR3(IRPCI,VRT(1,2),VRT(1,1),
     1            IRPDPD(IRPCI,13),IRPDPD(IRPCI,11),
     1            CORE(I030),CORE(I000),CORE(I010),CORE(I020))
      ENDIF
      CALL PUTLST(CORE(I030),1,IRPDPD(IRPCI,11),2,IRPCI,LWIC18)
C
      IF(LWIC18.NE.LWIC28)THEN
      CALL GETLST(CORE(I030),1,IRPDPD(IRPCI,11),2,IRPCI,LWIC28)
      IF(ICODE.EQ.1)THEN
      CALL SYMTR3(IRPCI,VRT(1,1),VRT(1,2),
     1            IRPDPD(IRPCI,13),IRPDPD(IRPCI,11),
     1            CORE(I030),CORE(I000),CORE(I010),CORE(I020))
      ELSE
      CALL SYMTR3(IRPCI,VRT(1,2),VRT(1,1),
     1            IRPDPD(IRPCI,13),IRPDPD(IRPCI,11),
     1            CORE(I030),CORE(I000),CORE(I010),CORE(I020))
      ENDIF
      CALL PUTLST(CORE(I030),1,IRPDPD(IRPCI,11),2,IRPCI,LWIC28)
      ENDIF
C
      IF(IDRLVL.GT.0.OR.EXTRNL_CCSD)THEN
      CALL GETLST(CORE(I030),1,IRPDPD(IRPCI,11),2,IRPCI,LISGV4)
      IF(ICODE.EQ.1)THEN
      CALL SYMTR3(IRPCI,VRT(1,1),VRT(1,2),
     1            IRPDPD(IRPCI,13),IRPDPD(IRPCI,11),
     1            CORE(I030),CORE(I000),CORE(I010),CORE(I020))
      ELSE
      CALL SYMTR3(IRPCI,VRT(1,2),VRT(1,1),
     1            IRPDPD(IRPCI,13),IRPDPD(IRPCI,11),
     1            CORE(I030),CORE(I000),CORE(I010),CORE(I020))
      ENDIF
      CALL PUTLST(CORE(I030),1,IRPDPD(IRPCI,11),2,IRPCI,LISGV4)
      ENDIF
C
      ELSE
C
C     Handle one distribution at a time.
C
      WRITE(6,1050)
C
      DO 110 IDIS=1,IRPDPD(IRPCI,11)
C
      CALL GETLST(CORE(I030),IDIS,1,2,IRPCI,LWIC18)
      IF(ICODE.EQ.1)THEN
      CALL SYMTR3(IRPCI,VRT(1,1),VRT(1,2),IRPDPD(IRPCI,13),1,
     1            CORE(I030),CORE(I000),CORE(I010),CORE(I020))
      ELSE
      CALL SYMTR3(IRPCI,VRT(1,2),VRT(1,1),IRPDPD(IRPCI,13),1,
     1            CORE(I030),CORE(I000),CORE(I010),CORE(I020))
      ENDIF
      CALL PUTLST(CORE(I030),IDIS,1,2,IRPCI,LWIC18)
C
      IF(LWIC18.NE.LWIC28)THEN
      CALL GETLST(CORE(I030),IDIS,1,2,IRPCI,LWIC28)
      IF(ICODE.EQ.1)THEN
      CALL SYMTR3(IRPCI,VRT(1,1),VRT(1,2),IRPDPD(IRPCI,13),1,
     1            CORE(I030),CORE(I000),CORE(I010),CORE(I020))
      ELSE
      CALL SYMTR3(IRPCI,VRT(1,2),VRT(1,1),IRPDPD(IRPCI,13),1,
     1            CORE(I030),CORE(I000),CORE(I010),CORE(I020))
      ENDIF
      CALL PUTLST(CORE(I030),IDIS,1,2,IRPCI,LWIC28)
      ENDIF
C
      IF(IDRLVL.GT.0.OR.EXTRNL_CCSD)THEN
      CALL GETLST(CORE(I030),IDIS,1,2,IRPCI,LISGV4)
      IF(ICODE.EQ.1)THEN
      CALL SYMTR3(IRPCI,VRT(1,1),VRT(1,2),IRPDPD(IRPCI,13),1,
     1            CORE(I030),CORE(I000),CORE(I010),CORE(I020))
      ELSE
      CALL SYMTR3(IRPCI,VRT(1,2),VRT(1,1),IRPDPD(IRPCI,13),1,
     1            CORE(I030),CORE(I000),CORE(I010),CORE(I020))
      ENDIF
      CALL PUTLST(CORE(I030),IDIS,1,2,IRPCI,LISGV4)
      ENDIF
C
  110 CONTINUE
C
      ENDIF
C
  120 CONTINUE
C
      RETURN
      END
