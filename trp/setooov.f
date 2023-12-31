      SUBROUTINE SETOOOV(CORE,MAXCOR,IUHF,ICODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CORE(1)
      INTEGER BUFSIZ,DIRPRD,POP,VRT,SYTYPL,SYTYPR,DISSIZ
      LOGICAL CHANGE
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
      COMMON /INFO/   NOCCO(2),NVRTO(2)
      COMMON /FLAGS/  IFLAGS(100)
      EQUIVALENCE(ICLLVL,IFLAGS( 2))
      EQUIVALENCE(IDRLVL,IFLAGS( 3))
      EQUIVALENCE(IREFNC,IFLAGS(11))
      EQUIVALENCE(IQRHFP,IFLAGS(32))
      EQUIVALENCE(IQRHFM,IFLAGS(33))
C
      COMMON /LISWI/  LWIC11,LWIC12,LWIC13,LWIC14,
     1                LWIC15,LWIC16,LWIC17,LWIC18,
     1                LWIC21,LWIC22,LWIC23,LWIC24,
     1                LWIC25,LWIC26,LWIC27,LWIC28,
     1                LWIC31,LWIC32,LWIC33,
     1                LWIC34,LWIC35,LWIC36,
     1                LWIC37,LWIC38,LWIC39,LWIC40,LWIC41,LWIC42
      COMMON /GAMLIS/ LISGO1,LISGO2,LISGO3,LISGO4,LISGV1,LISGV2,
     1                LISGV3,LISGV4
C
      CHANGE = .TRUE.
C
      WRITE(6,1000)
 1000 FORMAT(' @SETOOOV-I, Some OOOV lists are being transposed. ')
C
C     Routine to change certain ijka lists from (ij,ka) to (ka,ij) or
C     vice versa.
C
C     If IUHF=1 we deal with LWIC11, LWIC12, LWIC13, LWIC14
C     If IUHF=0 we deal with LWIC11 and LWIC14. If LWIC11 is 107, this
C     list does not exist for RHF, so we do nothing with it.
C
      BUFSIZ = 0
      DO 10 IRREP=1,NIRREP
      BUFSIZ = MAX0(BUFSIZ,
     1              IRPDPD(IRREP, 3),IRPDPD(IRREP, 4),IRPDPD(IRREP,14),
     1              IRPDPD(IRREP, 9),IRPDPD(IRREP,10),IRPDPD(IRREP,11),
     1              IRPDPD(IRREP,12))
   10 CONTINUE
C
C     Figure out if we have to change LWIC21, LWIC22, LWIC23, and LWIC24.
C
      NPASS = 1
      IF((ICLLVL.GE.14.AND.ICLLVL.LE.18).OR.ICLLVL.EQ.33.OR.
     &                                      ICLLVL.EQ.34)THEN
      NPASS = 2
      WRITE(6,1010)
 1010 FORMAT(' @SETOOOV-I, Two sets of OOOV lists are transposed. ')
      ENDIF
C
      DO 2000 IPASS=1,NPASS
C
      IF(ICODE.EQ.1)THEN
C
C     (ij,ka) ---> (ka,ij)
C
      IF(IUHF.EQ.1.OR.LWIC11.NE.107)THEN
      IF(IPASS.EQ.1)THEN
      LIST = LWIC11
      ELSE
      LIST = LWIC21
      ENDIF
      I000 = 1
      I010 = I000 + BUFSIZ
      DO 20 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP, 3)
      NDIS   = IRPDPD(IRREP,16)
      CALL GETTRN(CORE(I010),CORE(I000),DISSIZ,NDIS,2,IRREP,LIST)
C
      I010 = I010 + DISSIZ * NDIS
      NEED = I010 * IINTFP
C
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020)
      STOP
      ENDIF
   20 CONTINUE
C
C     Alter the list parameters.
C
      SYTYPL = 16
      SYTYPR =  3
      CALL NEWTYP(LIST,SYTYPL,SYTYPR,CHANGE)
C
      I000 = 1
      I010 = I000 + BUFSIZ
      DO 30 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP,16)
      NDIS   = IRPDPD(IRREP, 3)
      CALL PUTLST(CORE(I010),1,NDIS,2,IRREP,LIST)
C
      I010 = I010 + DISSIZ * NDIS
      NEED = I010 * IINTFP
C
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020)
      STOP
      ENDIF
   30 CONTINUE
C
      ENDIF
C
      IF(IUHF.EQ.1)THEN
      IF(IPASS.EQ.1)THEN
      LIST = LWIC12
      ELSE
      LIST = LWIC22
      ENDIF
      I000 = 1
      I010 = I000 + BUFSIZ
      DO 40 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP, 4)
      NDIS   = IRPDPD(IRREP,17)
      CALL GETTRN(CORE(I010),CORE(I000),DISSIZ,NDIS,2,IRREP,LIST)
C
      I010 = I010 + DISSIZ * NDIS
      NEED = I010 * IINTFP
C
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020)
      STOP
      ENDIF
   40 CONTINUE
C
C     Alter the list parameters.
C
      SYTYPL = 17
      SYTYPR =  4
      CALL NEWTYP(LIST,SYTYPL,SYTYPR,CHANGE)
C
      I000 = 1
      I010 = I000 + BUFSIZ
      DO 50 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP,17)
      NDIS   = IRPDPD(IRREP, 4)
      CALL PUTLST(CORE(I010),1,NDIS,2,IRREP,LIST)
C
      I010 = I010 + DISSIZ * NDIS
      NEED = I010 * IINTFP
C
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020)
      STOP
      ENDIF
   50 CONTINUE
C
      IF(IPASS.EQ.1)THEN
      LIST = LWIC13
      ELSE
      LIST = LWIC23
      ENDIF
      I000 = 1
      I010 = I000 + BUFSIZ
      DO 60 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP,14)
      NDIS   = IRPDPD(IRREP,11)
      CALL GETTRN(CORE(I010),CORE(I000),DISSIZ,NDIS,2,IRREP,LIST)
C
      I010 = I010 + DISSIZ * NDIS
      NEED = I010 * IINTFP
C
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020)
      STOP
      ENDIF
   60 CONTINUE
C
C     Alter the list parameters.
C
      SYTYPL = 11
      SYTYPR = 14
      CALL NEWTYP(LIST,SYTYPL,SYTYPR,CHANGE)
C
      I000 = 1
      I010 = I000 + BUFSIZ
      DO 70 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP,11)
      NDIS   = IRPDPD(IRREP,14)
      CALL PUTLST(CORE(I010),1,NDIS,2,IRREP,LIST)
C
      I010 = I010 + DISSIZ * NDIS
      NEED = I010 * IINTFP
C
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020)
      STOP
      ENDIF
   70 CONTINUE
C
      ENDIF
C
      IF(IPASS.EQ.1)THEN
      LIST = LWIC14
      ELSE
      LIST = LWIC24
      ENDIF
      I000 = 1
      I010 = I000 + BUFSIZ
      DO 80 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP,14)
      NDIS   = IRPDPD(IRREP,18)
      CALL GETTRN(CORE(I010),CORE(I000),DISSIZ,NDIS,2,IRREP,LIST)
C
      I010 = I010 + DISSIZ * NDIS
      NEED = I010 * IINTFP
C
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020)
      STOP
      ENDIF
   80 CONTINUE
C
C     Alter the list parameters.
C
      SYTYPL = 18
      SYTYPR = 14
      CALL NEWTYP(LIST,SYTYPL,SYTYPR,CHANGE)
C
      I000 = 1
      I010 = I000 + BUFSIZ
      DO 90 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP,18)
      NDIS   = IRPDPD(IRREP,14)
      CALL PUTLST(CORE(I010),1,NDIS,2,IRREP,LIST)
C
      I010 = I010 + DISSIZ * NDIS
      NEED = I010 * IINTFP
C
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020)
      STOP
      ENDIF
   90 CONTINUE
C
      ENDIF
C
      IF(ICODE.EQ.2)THEN
C
C     (ka,ij) ---> (ij,ka)
C
      IF(IUHF.EQ.1.OR.LWIC11.NE.107)THEN
      IF(IPASS.EQ.1)THEN
      LIST = LWIC11
      ELSE
      LIST = LWIC21
      ENDIF
      I000 = 1
      I010 = I000 + BUFSIZ
      DO 520 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP,16)
      NDIS   = IRPDPD(IRREP, 3)
      CALL GETTRN(CORE(I010),CORE(I000),DISSIZ,NDIS,2,IRREP,LIST)
C
      I010 = I010 + DISSIZ * NDIS
      NEED = I010 * IINTFP
C
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020)
      STOP
      ENDIF
  520 CONTINUE
C
C     Alter the list parameters.
C
      SYTYPL =  3
      SYTYPR = 16
      CALL NEWTYP(LIST,SYTYPL,SYTYPR,CHANGE)
C
      I000 = 1
      I010 = I000 + BUFSIZ
      DO 530 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP, 3)
      NDIS   = IRPDPD(IRREP,16)
      CALL PUTLST(CORE(I010),1,NDIS,2,IRREP,LIST)
C
      I010 = I010 + DISSIZ * NDIS
      NEED = I010 * IINTFP
C
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020)
      STOP
      ENDIF
  530 CONTINUE
C
      ENDIF
C
      IF(IUHF.EQ.1)THEN
      IF(IPASS.EQ.1)THEN
      LIST = LWIC12
      ELSE
      LIST = LWIC22
      ENDIF
      I000 = 1
      I010 = I000 + BUFSIZ
      DO 540 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP,17)
      NDIS   = IRPDPD(IRREP, 4)
      CALL GETTRN(CORE(I010),CORE(I000),DISSIZ,NDIS,2,IRREP,LIST)
C
      I010 = I010 + DISSIZ * NDIS
      NEED = I010 * IINTFP
C
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020)
      STOP
      ENDIF
  540 CONTINUE
C
C     Alter the list parameters.
C
      SYTYPL =  4
      SYTYPR = 17
      CALL NEWTYP(LIST,SYTYPL,SYTYPR,CHANGE)
C
      I000 = 1
      I010 = I000 + BUFSIZ
      DO 550 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP, 4)
      NDIS   = IRPDPD(IRREP,17)
      CALL PUTLST(CORE(I010),1,NDIS,2,IRREP,LIST)
C
      I010 = I010 + DISSIZ * NDIS
      NEED = I010 * IINTFP
C
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020)
      STOP
      ENDIF
  550 CONTINUE
C
      IF(IPASS.EQ.1)THEN
      LIST = LWIC13
      ELSE
      LIST = LWIC23
      ENDIF
      I000 = 1
      I010 = I000 + BUFSIZ
      DO 560 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP,11)
      NDIS   = IRPDPD(IRREP,14)
      CALL GETTRN(CORE(I010),CORE(I000),DISSIZ,NDIS,2,IRREP,LIST)
C
      I010 = I010 + DISSIZ * NDIS
      NEED = I010 * IINTFP
C
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020)
      STOP
      ENDIF
  560 CONTINUE
C
C     Alter the list parameters.
C
      SYTYPL = 14
      SYTYPR = 11
      CALL NEWTYP(LIST,SYTYPL,SYTYPR,CHANGE)
C
      I000 = 1
      I010 = I000 + BUFSIZ
      DO 570 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP,14)
      NDIS   = IRPDPD(IRREP,11)
      CALL PUTLST(CORE(I010),1,NDIS,2,IRREP,LIST)
C
      I010 = I010 + DISSIZ * NDIS
      NEED = I010 * IINTFP
C
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020)
      STOP
      ENDIF
  570 CONTINUE
C
      ENDIF
C
      IF(IPASS.EQ.1)THEN
      LIST = LWIC14
      ELSE
      LIST = LWIC24
      ENDIF
      I000 = 1
      I010 = I000 + BUFSIZ
      DO 580 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP,18)
      NDIS   = IRPDPD(IRREP,14)
      CALL GETTRN(CORE(I010),CORE(I000),DISSIZ,NDIS,2,IRREP,LIST)
C
      I010 = I010 + DISSIZ * NDIS
      NEED = I010 * IINTFP
C
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020)
      STOP
      ENDIF
  580 CONTINUE
C
C     Alter the list parameters.
C
      SYTYPL = 14
      SYTYPR = 18
      CALL NEWTYP(LIST,SYTYPL,SYTYPR,CHANGE)
C
      I000 = 1
      I010 = I000 + BUFSIZ
      DO 590 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP,14)
      NDIS   = IRPDPD(IRREP,18)
      CALL PUTLST(CORE(I010),1,NDIS,2,IRREP,LIST)
C
      I010 = I010 + DISSIZ * NDIS
      NEED = I010 * IINTFP
C
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020)
      STOP
      ENDIF
  590 CONTINUE
C
      ENDIF
 2000 CONTINUE
C
      IF(ICODE.NE.1.AND.ICODE.NE.2)THEN
      WRITE(6,1040) ICODE
      STOP
      ENDIF
C
 1020 FORMAT(' @SETOOOV-I, Insufficient memory to continue. ')
 1030 FORMAT(' @SETOOOV-I, Memory required ',I10,' double words. ')
 1040 FORMAT(' @SETOOOV-I, Invalid value of ICODE. ',I10,' was given\n',
     1       '             but must be 1 or 2. ')
      RETURN
      END
