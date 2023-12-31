      SUBROUTINE PCCD_INITIN(ICORE,MAXCOR,IUHF)
C
C INITIALIZES THE T2 INCREMENT LIST WITH THE <AB//IJ> INTEGRALS.
C
      IMPLICIT INTEGER (A-Z)
C
      DIMENSION ICORE(MAXCOR)
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),IJUNK(18)
C
      IBOT=3
      DO 15 ISPIN=IBOT,3
          NSIZE=0
          DO 10 IRREP=1,NIRREP
             NUMSYM=IRPDPD(IRREP,ISYTYP(2,13+ISPIN))
             DSZSYM=IRPDPD(IRREP,ISYTYP(1,13+ISPIN))
             NSIZE=NSIZE+NUMSYM*DSZSYM
10        CONTINUE
          CALL GETALL(ICORE,NSIZE,1,13+ISPIN)
          CALL PUTALL(ICORE,NSIZE,1,60+ISPIN)
15    CONTINUE

      RETURN
      END
