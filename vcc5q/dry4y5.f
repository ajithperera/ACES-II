      SUBROUTINE DRY4Y5(ICORE,MAXCOR,IUHF,ITYPE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD
      DIMENSION ICORE(MAXCOR)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NJUNK(18)
C
C     ITYPE = 0 : Do OOOO intermediate (Y4)
C     ITYPE = 1 : Do VVVV intermediate (Y5)
C
      IF(ITYPE.NE.0.AND.ITYPE.NE.1)THEN
      WRITE(6,1010) ITYPE
 1010 FORMAT(' @DRY4Y5-I, Invalid ITYPE. Given ',I10,' . Must be 0/1.')
      STOP 'DRY4Y5'
      ENDIF
C
      IF(ITYPE.EQ.0)THEN
C
      WRITE(6,1020)
 1020 FORMAT(' @DRY4Y5-I, Y4 intermediate being computed. ')
C
      DO   20 ISPIN=1,IUHF+1
C
      DO   10 IRREP=1,NIRREP
C
C     I000 Target
C     I010 
C     I020
C
      I000 = 1
      I010 = I000 + IINTFP * IRPDPD(IRREP,2 + ISPIN) ** 2
      I020 = I010 + IINTFP * IRPDPD(IRREP,    ISPIN) * 
     1                       IRPDPD(IRREP,2 + ISPIN)
      I030 = I020 + IINTFP * IRPDPD(IRREP,    ISPIN) * 
     1                       IRPDPD(IRREP,2 + ISPIN)
C
      IF(I030.GT.MAXCOR)THEN
      WRITE(6,1030) I030,MAXCOR
 1030 FORMAT(' @DRY4Y5-I, Insufficient memory. Got ',I10,
     1       ' Need ',I10)
      STOP
      ENDIF
C
      NROWS = IRPDPD(IRREP,2 + ISPIN)
      NCOLS = IRPDPD(IRREP,2 + ISPIN)
      NOVER = IRPDPD(IRREP,    ISPIN)
C
      IF(MIN(NROWS,NCOLS,NOVER).GT.0)THEN
C
      CALL GETLST(ICORE(I010),1,NCOLS,1,IRREP,43 + ISPIN)
      CALL GETLST(ICORE(I020),1,NCOLS,1,IRREP,43 + ISPIN)
C
      CALL ZERO(ICORE(I000),NROWS*NCOLS)
      CALL XGEMM('T','N',NROWS,NCOLS,NOVER,1.0D+00,
     1           ICORE(I010),NOVER,ICORE(I020),NOVER,
     1           0.0D+00,ICORE(I000),NROWS)
C
      CALL PUTLIST(ICORE(I000),1,NCOLS,99,IRREP,6 + ISPIN)
C
      ENDIF
C
   10 CONTINUE
   20 CONTINUE
C
      DO   30 IRREP=1,NIRREP
C
      I000 = 1
      I010 = I000 + IINTFP * IRPDPD(IRREP,14) ** 2
      I020 = I010 + IINTFP * IRPDPD(IRREP,13) * 
     1                       IRPDPD(IRREP,14)
      I030 = I020 + IINTFP * IRPDPD(IRREP,13) * 
     1                       IRPDPD(IRREP,14)
C
      IF(I030.GT.MAXCOR)THEN
      WRITE(6,1030) I030,MAXCOR
      STOP
      ENDIF
C
      NROWS = IRPDPD(IRREP,14)
      NCOLS = IRPDPD(IRREP,14)
      NOVER = IRPDPD(IRREP,13)
C
      IF(MIN(NROWS,NCOLS,NOVER).GT.0)THEN
C
      CALL GETLST(ICORE(I010),1,NCOLS,1,IRREP,46)
      CALL GETLST(ICORE(I020),1,NCOLS,1,IRREP,46)
C
      CALL ZERO(ICORE(I000),NROWS*NCOLS)
      CALL XGEMM('T','N',NROWS,NCOLS,NOVER,1.0D+00,
     1           ICORE(I010),NOVER,ICORE(I020),NOVER,
     1           0.0D+00,ICORE(I000),NROWS)
C
      CALL PUTLIST(ICORE(I000),1,NCOLS,99,IRREP,9)
C
      ENDIF
C
   30 CONTINUE
C
      ENDIF
C
      IF(ITYPE.EQ.1)THEN
C
      WRITE(6,1040)
 1040 FORMAT(' @DRY4Y5-I, Y5 intermediate being computed. ')
C
      DO  120 ISPIN=1,IUHF+1
C
      DO  110 IRREP=1,NIRREP
C
      I000 = 1
      I010 = I000 + IINTFP * IRPDPD(IRREP,    ISPIN) ** 2
      I020 = I010 + IINTFP * IRPDPD(IRREP,    ISPIN) * 
     1                       IRPDPD(IRREP,2 + ISPIN)
      I030 = I020 + IINTFP * IRPDPD(IRREP,    ISPIN) * 
     1                       IRPDPD(IRREP,2 + ISPIN)
C
      IF(I030.GT.MAXCOR)THEN
      WRITE(6,1030) I030,MAXCOR
      STOP
      ENDIF
C
      NROWS = IRPDPD(IRREP,    ISPIN)
      NCOLS = IRPDPD(IRREP,    ISPIN)
      NOVER = IRPDPD(IRREP,2 + ISPIN)
C
      IF(MIN(NROWS,NCOLS,NOVER).GT.0)THEN
C
      CALL GETLST(ICORE(I010),1,NOVER,1,IRREP,43 + ISPIN)
      CALL GETLST(ICORE(I020),1,NOVER,1,IRREP,43 + ISPIN)
C
      CALL ZERO(ICORE(I000),NROWS*NCOLS)
      CALL XGEMM('N','T',NROWS,NCOLS,NOVER,1.0D+00,
     1           ICORE(I010),NROWS,ICORE(I020),NCOLS,
     1           0.0D+00,ICORE(I000),NROWS)
C
C     CALL CHKSUM(ICORE(I000),NROWS*NCOLS)
C
      CALL PUTLIST(ICORE(I000),1,NCOLS,99,IRREP,    ISPIN)
C     CALL GETLIST(ICORE(I000),1,NCOLS,99,IRREP,    ISPIN)
C     CALL CHKSUM(ICORE(I000),NROWS*NCOLS)
C
      ENDIF
C
  110 CONTINUE
  120 CONTINUE
C
      DO  130 IRREP=1,NIRREP
C
      I000 = 1
      I010 = I000 + IINTFP * IRPDPD(IRREP,13) ** 2
      I020 = I010 + IINTFP * IRPDPD(IRREP,13) * 
     1                       IRPDPD(IRREP,14)
      I030 = I020 + IINTFP * IRPDPD(IRREP,13) * 
     1                       IRPDPD(IRREP,14)
C
      IF(I030.GT.MAXCOR)THEN
      WRITE(6,1030) I030,MAXCOR
      STOP
      ENDIF
C
      NROWS = IRPDPD(IRREP,13)
      NCOLS = IRPDPD(IRREP,13)
      NOVER = IRPDPD(IRREP,14)
C
      IF(MIN(NROWS,NCOLS,NOVER).GT.0)THEN
C
      CALL GETLST(ICORE(I010),1,NOVER,1,IRREP,46)
      CALL GETLST(ICORE(I020),1,NOVER,1,IRREP,46)
C
      CALL ZERO(ICORE(I000),NROWS*NCOLS)
      CALL XGEMM('N','T',NROWS,NCOLS,NOVER,1.0D+00,
     1           ICORE(I010),NROWS,ICORE(I020),NCOLS,
     1           0.0D+00,ICORE(I000),NROWS)
C
c     CALL CHKSUM(ICORE(I000),NROWS*NCOLS)
C
      CALL PUTLIST(ICORE(I000),1,NCOLS,99,IRREP,3)
C
      ENDIF
C
  130 CONTINUE
C
      ENDIF
C
      RETURN
      END