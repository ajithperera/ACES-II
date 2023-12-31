
C   INITIATES THE LIST FOR THE T1 AMPLITUDES AND INCREMENTS
C   AND ZEROS THEM

      SUBROUTINE INITSN(ICORE,MAXCOR,IUHF)
      IMPLICIT INTEGER(A-Z)
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,QCISD,UCC
      DIMENSION ICORE(MAXCOR)
      COMMON /SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPR(8,8)
      COMMON /SYM/ POP1(8),POP2(8),VRT1(8),VRT2(8),NT(2),
     &             NF1AA,NF1BB,NF2AA,NF2BB
      COMMON /FLAGS/ IFLAGS(100)

C     CREATE LISTS (FIRST CYCLE)
      IF (.NOT.SING1) THEN
         CALL UPDMOI(1,NT(1),1,90,0,0)
         CALL UPDMOI(1,NT(1),3,90,0,0)
         IF (IUHF.EQ.1) THEN
         CALL UPDMOI(1,NT(2),2,90,0,0)
         CALL UPDMOI(1,NT(2),4,90,0,0)
         END IF
      END IF
      IF ((IFLAGS(77)+IFLAGS(38)).NE.0) THEN
         CALL UPDMOI(1,NT(1),3,90,0,0)
         IF (IUHF.EQ.1) THEN
         CALL UPDMOI(1,NT(2),4,90,0,0)
         END IF
      END IF

C     ZERO THE T1 INCREMENTS (ALL CYCLES)
      CALL ZERLST(ICORE,NT(1),1,1,3,90)
      IF (IUHF.EQ.1) CALL ZERLST(ICORE,NT(2),1,1,4,90)

      RETURN
      END

