      SUBROUTINE PCCD_GAMMA3(ICORE,MAXCOR,IUHF)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYT,POP,VRT
      DIMENSION ICORE(MAXCOR)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      DATA ONE /1.0D+0/

      LISTG=113
      LISTT1=153

      DO 110 IRREP=1,NIRREP
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT1))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT1))
        I001=1
        I002=I001
        I003=I002

        IF(MIN(NUMSYT,DISSYT).NE.0)THEN
         I004=I003+IINTFP*NUMSYT*NUMSYT
         IF(I004.LT.MAXCOR) THEN
          CALL PCCD_G3ALL(ICORE(I001),ICORE(I002),ICORE(I003),
     &                    3,DISSYT,NUMSYT,LISTT1,LISTG,IRREP)
         ELSE
          CALL INSMEM('PCCD_GAMMA3',I004,MAXCOR) 
         ENDIF
        ENDIF
110    CONTINUE
C
#ifdef _DEBUG_LVL0
      call checkgam1(icore,13,113,ONE,IUHF,1,POP)
#endif 

      RETURN
      END
