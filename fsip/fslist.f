      
      SUBROUTINE FSLIST(ICORE,MAXCOR,IUHF,SECTOR)
C
C CREATES REQUIRED LISTS FOR FOCK-SPACE CALCULATIONS
C
      IMPLICIT INTEGER (A-Z)
      CHARACTER*2 SECTOR
      DIMENSION ICORE(MAXCOR)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /FSSYMPOP/ FSDPDAN(8,22),FSDPDNA(8,22),FSDPDAA(8,22),
     &                  FSDPDIN(8,22),FSDPDNI(8,22),FSDPDII(8,22),
     &                  FSDPDAI(8,22),FSDPDIA(8,22)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /FSSYM/ POPA(8,2),POPI(8,2),VRTA(8,2),VRTI(8,2),
     &               NTAN(2),NTNA(2),NTAA(2),
     &               NTIN(2),NTNI(2),NTII(2),
     &               NTAI(2),NTIA(2),
     &               NF1AN(2),NF1AA(2),NF1IN(2),NF1II(2),NF1AI(2),
     &               NF2AN(2),NF2AA(2),NF2IN(2),NF2II(2),NF2AI(2)
C
      IMODE=0
C
      IF(SECTOR.EQ.'01')THEN
       IF(IUHF.NE.0)THEN
        CALL INIPCK(1,21,16,1,IMODE,0,1)
        CALL INIPCK(1,22,17,2,IMODE,0,1)
        CALL INIPCK(1,16,22,3,IMODE,0,1)
        CALL INIPCK(1,21,16,81,IMODE,0,1)
        CALL INIPCK(1,22,17,82,IMODE,0,1)
        CALL INIPCK(1,16,22,83,IMODE,0,1)
        CALL INIPCK(1,21,17,84,IMODE,0,1)
        CALL UPDMOI(1,NFMI(2),2,94,0,0)
        CALL UPDMOI(1,NFMI(2),4,94,0,0)
        CALL UPDMOI(1,NFMI(2),4,91,0,0)
        CALL INIPCK(1,3,16,96,IMODE,0,1)
        CALL INIPCK(1,4,17,97,IMODE,0,1)
        CALL INIPCK(1,14,11,98,IMODE,0,1)
        CALL INIPCK(1,3,16,196,IMODE,0,1)
        CALL INIPCK(1,4,17,197,IMODE,0,1)
        CALL INIPCK(1,14,11,198,IMODE,0,1)
       ENDIF
       CALL INIPCK(1,21,17,4,IMODE,0,1)
       CALL INIPCK(1,14,18,99,IMODE,0,1)
       CALL INIPCK(1,14,18,199,IMODE,0,1)
       CALL UPDMOI(1,NFMI(1),1,94,0,0)
       CALL UPDMOI(1,NFMI(1),3,94,0,0)
       CALL UPDMOI(1,NFMI(1),3,91,0,0)
      ELSEIF(SECTOR.EQ.'10')THEN
       write(6,*)' sorry, not coded '
      ELSEIF(SECTOR.Eq.'11')THEN
       write(6,*)' sorry, not coded '
      ELSEIF(SECTOR.EQ.'20')THEN
       write(6,*)' sorry, not coded '
      ELSEIF(SECTOR.EQ.'02')THEN
       IF(IUHF.NE.0)THEN
        CALL INIPCK(1,3,3,111,IMODE,0,1)
        CALL INIPCK(1,4,4,112,IMODE,0,1)
       ENDIF
       CALL INIPCK (1,14,14,113,IMODE,0,1)
      ELSE
       write(6,*)' sorry, not known'
      ENDIF
      RETURN
      END
