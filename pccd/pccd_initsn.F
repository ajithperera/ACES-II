      SUBROUTINE PCCD_INITSN(ICORE,MAXCOR,IUHF,ICYCLE,SIDE)

      IMPLICIT INTEGER(A-Z)
      DIMENSION ICORE(MAXCOR)
      LOGICAL START 
      CHARACTER*1 SIDE
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPR(8,8)
      COMMON /SYM/ POP1(8),POP2(8),VRT1(8),VRT2(8),NT(2),
     &             NF1AA,NF1BB,NF2AA,NF2BB

      START = ICYCLE .EQ. 1
      IF (START .AND. SIDE .EQ. "T") THEN
         CALL UPDMOI(1,NT(1),1,90,0,0)
         CALL UPDMOI(1,NT(1),3,90,0,0)
         IF (IUHF.EQ.1) THEN
            CALL UPDMOI(1,NT(2),2,90,0,0)
            CALL UPDMOI(1,NT(2),4,90,0,0)
         END IF
      END IF

      CALL ZERLST(ICORE,NT(1),1,1,3,90)
      IF (IUHF.EQ.1) CALL ZERLST(ICORE,NT(2),1,1,4,90)

      RETURN
      END

