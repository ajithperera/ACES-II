      SUBROUTINE SETVOVO(CORE,MAXCOR,IUHF,ICODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CORE(1)
      INTEGER BUFSIZ,DIRPRD,POP,VRT,SYTYPL,SYTYPR,DISSIZ
      LOGICAL CHANGE
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     1                DIRPRD(8,8)
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
      COMMON /FILES/  LUOUT,MOINTS
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
      WRITE(LUOUT,1000)
 1000 FORMAT(' @SETVOVO-I, Some VOVO lists are being reordered. ')
C
C     Routine to change the W(MB,EJ) lists (54-59) from E,M,B,J to E,B,M,J
C     or vice versa. Which resort is done depends on ICODE.
C
C     If IUHF=1 we deal with LWIC37,LWIC38,LWIC39,LWIC40,LWIC41,LWIC42
C     If IUHF=0 we deal with LWIC37,       LWIC39,       LWIC41
C     list does not exist for RHF, so we do nothing with it.
C
      LWIC37 = 54
      LWIC38 = 55
      LWIC39 = 56
      LWIC40 = 57
      LWIC41 = 58
      LWIC42 = 59
C
      BUFSIZ = 0
      DO 10 IRREP=1,NIRREP
      BUFSIZ = MAX0(BUFSIZ,
     1              IRPDPD(IRREP,19),IRPDPD(IRREP,20),IRPDPD(IRREP,13),
     1              IRPDPD(IRREP,21),IRPDPD(IRREP,22),IRPDPD(IRREP,14),
     1              IRPDPD(IRREP, 9),IRPDPD(IRREP,10),IRPDPD(IRREP,11),
     1              IRPDPD(IRREP,12))
   10 CONTINUE
C
      DO 2000 ISPIN=1,IUHF+1
C
      IF(ICODE.EQ.1)THEN
C
C     E,M,B,J ---> E,B,M,J
C     e,m,b,j ---> e,b,m,j
C
      IF(ISPIN.EQ.1)THEN
      LIST = LWIC37
      ELSE
      LIST = LWIC38
      ENDIF
      LEN = 0
      DO 20 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP, 8 + ISPIN)
      NDIS   = IRPDPD(IRREP, 8 + ISPIN)
      LEN    = LEN + DISSIZ * NDIS
   20 CONTINUE
C
      I000 = 1
      I010 = I000 + BUFSIZ
      I020 = I010 + LEN
      I030 = I020 + LEN
      NEED = I030 * IINTFP
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1020)
      STOP
      ENDIF
C
C     Read the whole list.
C
      CALL GETALL(CORE(I010),LEN,1,LIST)
C
C     Reorder the array.
C
      CALL SSTGEN(CORE(I010),CORE(I020),LEN,VRT(1,ISPIN),POP(1,ISPIN),
     1            VRT(1,ISPIN),POP(1,ISPIN),CORE(I000),1,'1324')
C
C     Alter the list parameters.
C
      SYTYPL = 18 + ISPIN
      SYTYPR = 20 + ISPIN
      CALL NEWTYP(LIST,SYTYPL,SYTYPR,CHANGE)
C
C     Write the reordered array back to the list.
C
      CALL PUTALL(CORE(I020),LEN,1,LIST)
C
C     E,M,b,j ---> E,b,M,j
C     e,m,B,J ---> e,B,m,J
C
      IF(ISPIN.EQ.1)THEN
      LIST = LWIC39
      ELSE
      LIST = LWIC40
      ENDIF
      LEN = 0
      DO 30 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP, 8 + ISPIN)
      NDIS   = IRPDPD(IRREP,11 - ISPIN)
      LEN    = LEN + DISSIZ * NDIS
   30 CONTINUE
C
      I000 = 1
      I010 = I000 + BUFSIZ
      I020 = I010 + LEN
      I030 = I020 + LEN
      NEED = I030 * IINTFP
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1020)
      STOP
      ENDIF
C
      CALL GETALL(CORE(I010),LEN,1,LIST)
      CALL SSTGEN(CORE(I010),CORE(I020),LEN,VRT(1,ISPIN),POP(1,ISPIN),
     1            VRT(1,3-ISPIN),POP(1,3-ISPIN),CORE(I000),1,'1324')
C
      SYTYPL = 13
      SYTYPR = 14
      CALL NEWTYP(LIST,SYTYPL,SYTYPR,CHANGE)
      CALL PUTALL(CORE(I020),LEN,1,LIST)
C
C     E,m,B,j ---> E,B,m,j
C     e,M,b,J ---> e,b,M,J
C
      IF(ISPIN.EQ.1)THEN
      LIST = LWIC41
      ELSE
      LIST = LWIC42
      ENDIF
      LEN = 0
      DO 40 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP,10 + ISPIN)
      NDIS   = IRPDPD(IRREP,10 + ISPIN)
      LEN    = LEN + DISSIZ * NDIS
   40 CONTINUE
C
      I000 = 1
      I010 = I000 + BUFSIZ
      I020 = I010 + LEN
      I030 = I020 + LEN
      NEED = I030 * IINTFP
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1020)
      STOP
      ENDIF
C
      CALL GETALL(CORE(I010),LEN,1,LIST)
      CALL SSTGEN(CORE(I010),CORE(I020),LEN,VRT(1,ISPIN),POP(1,3-ISPIN),
     1            VRT(1,ISPIN),POP(1,3-ISPIN),CORE(I000),1,'1324')
C
      SYTYPL = 18 + ISPIN
      SYTYPR = 23 - ISPIN
      CALL NEWTYP(LIST,SYTYPL,SYTYPR,CHANGE)
      CALL PUTALL(CORE(I020),LEN,1,LIST)
C
      ENDIF
C
      IF(ICODE.EQ.2)THEN
C
C     E,B,M,J ---> E,M,B,J
C     e,b,m,j ---> e,m,b,j
C
      IF(ISPIN.EQ.1)THEN
      LIST = LWIC37
      ELSE
      LIST = LWIC38
      ENDIF
      LEN = 0
      DO 120 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP, 8 + ISPIN)
      NDIS   = IRPDPD(IRREP, 8 + ISPIN)
      LEN    = LEN + DISSIZ * NDIS
  120 CONTINUE
C
      I000 = 1
      I010 = I000 + BUFSIZ
      I020 = I010 + LEN
      I030 = I020 + LEN
      NEED = I030 * IINTFP
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1020)
      STOP
      ENDIF
C
C     Read the whole list.
C
      CALL GETALL(CORE(I010),LEN,1,LIST)
C
C     Reorder the array.
C
      CALL SSTGEN(CORE(I010),CORE(I020),LEN,VRT(1,ISPIN),VRT(1,ISPIN),
     1            POP(1,ISPIN),POP(1,ISPIN),CORE(I000),1,'1324')
C
C     Alter the list parameters.
C
      SYTYPL =  8 + ISPIN
      SYTYPR =  8 + ISPIN
      CALL NEWTYP(LIST,SYTYPL,SYTYPR,CHANGE)
C
C     Write the reordered array back to the list.
C
      CALL PUTALL(CORE(I020),LEN,1,LIST)
C
C     E,b,M,j ---> E,M,b,j
C     e,B,m,J ---> e,m,B,J
C
      IF(ISPIN.EQ.1)THEN
      LIST = LWIC39
      ELSE
      LIST = LWIC40
      ENDIF
      LEN = 0
      DO 130 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP, 8 + ISPIN)
      NDIS   = IRPDPD(IRREP,11 - ISPIN)
      LEN    = LEN + DISSIZ * NDIS
  130 CONTINUE
C
      I000 = 1
      I010 = I000 + BUFSIZ
      I020 = I010 + LEN
      I030 = I020 + LEN
      NEED = I030 * IINTFP
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1020)
      STOP
      ENDIF
C
      CALL GETALL(CORE(I010),LEN,1,LIST)
      CALL SSTGEN(CORE(I010),CORE(I020),LEN,VRT(1,ISPIN),VRT(1,3-ISPIN),
     1            POP(1,ISPIN),POP(1,3-ISPIN),CORE(I000),1,'1324')
C
      SYTYPL =  8 + ISPIN
      SYTYPR = 11 - ISPIN
      CALL NEWTYP(LIST,SYTYPL,SYTYPR,CHANGE)
      CALL PUTALL(CORE(I020),LEN,1,LIST)
C
C     E,B,m,j ---> E,m,B,j
C     e,b,M,J ---> e,M,b,J
C
      IF(ISPIN.EQ.1)THEN
      LIST = LWIC41
      ELSE
      LIST = LWIC42
      ENDIF
      LEN = 0
      DO 140 IRREP=1,NIRREP
      DISSIZ = IRPDPD(IRREP,10 + ISPIN)
      NDIS   = IRPDPD(IRREP,10 + ISPIN)
      LEN    = LEN + DISSIZ * NDIS
  140 CONTINUE
C
      I000 = 1
      I010 = I000 + BUFSIZ
      I020 = I010 + LEN
      I030 = I020 + LEN
      NEED = I030 * IINTFP
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1020)
      STOP
      ENDIF
C
      CALL GETALL(CORE(I010),LEN,1,LIST)
      CALL SSTGEN(CORE(I010),CORE(I020),LEN,VRT(1,ISPIN),VRT(1,ISPIN),
     1            POP(1,3-ISPIN),POP(1,3-ISPIN),CORE(I000),1,'1324')
C
      SYTYPL = 10 + ISPIN
      SYTYPR = 10 + ISPIN
      CALL NEWTYP(LIST,SYTYPL,SYTYPR,CHANGE)
      CALL PUTALL(CORE(I020),LEN,1,LIST)
C
      ENDIF
 2000 CONTINUE
C
      IF(ICODE.NE.1.AND.ICODE.NE.2)THEN
      WRITE(LUOUT,1040) ICODE
      STOP
      ENDIF
C
 1020 FORMAT(' @SETVOVO-I, Insufficient memory to continue. ')
 1030 FORMAT(' @SETVOVO-I, Memory required ',I10,' double words. ')
 1040 FORMAT(' @SETVOVO-I, Invalid value of ICODE. ',I10,' was given\n',
     1       '             but must be 1 or 2.')
      RETURN
      END
