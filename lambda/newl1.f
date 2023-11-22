      SUBROUTINE NEWL1(L1,ICORE,MAXCOR,ISPIN)
C
C THIS ROUTINE PICKS UP THE FINAL T1 OR T1 INCREMENTS, DENOMINATOR
C WEIGHTS THEM AND THEN OVERWRITES THE T1 INCREMENT LIST WITH THE NEW VALUES.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,L1
C
      DIMENSION ICORE(MAXCOR)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP1(8),POP2(8),VRT1(8),VRT2(8),NT(2),NF1AA,
     &           NF1BB,NF2AA,NF2B
C
      DATA ONE /1.0D0/
C
      LSTDEN=63+ISPIN
      LSTINC=2+ISPIN
      NSIZE=NT(ISPIN)
      I010=1
      I020=I010+NSIZE*IINTFP
      I030=I020+NSIZE*IINTFP
      IF(I030.GT.MAXCOR)CALL INSMEM('NEWL1',I020,MAXCOR)
      CALL GETLST(ICORE(I010),1,1,1,9,LSTDEN)
      CALL VECDIV(L1,ICORE(I010),L1,NSIZE)
C
      RETURN
      END
