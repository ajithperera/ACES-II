      
      SUBROUTINE PUTALLS2(Z,LENGTH,POP1, POP2, IRREPX,LSTNUM)
C
C THIS ROUTINE WRITES OUT ALL IRREPS OF A SYMMETRY PACKED DISTRIBUTION.
C
C      Z  - VECTOR HOLDING THE LIST.
C  LENGTH - THE TOTAL LENGTH.
C  POP1, POP2: DETERMINE THE NUMDIS LENGTH OF VECTOR FOR EACH IRREP
C  IRREPX - THE OVERALL SYMMETRY OF THE COEFFICIENTS
C  LSTNUM - THE LIST NUMBER WHICH Z IS WRITTEN TO.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DIMENSION POP1(8), POP2(8), NUMDSS(8)
      DOUBLE PRECISION Z(LENGTH)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),D(18)
C
C FIRST DETERMINE NUMDSS FROM POP1 AND POP2
C
      DO 10 XIRREP = 1, NIRREP
         NUMDSS(XIRREP) = 0
         DO 20 IRREP1 = 1, NIRREP
            IRREP2 = DIRPRD(IRREP1, XIRREP)
            NUMDSS(XIRREP) = NUMDSS(XIRREP)+POP1(IRREP1)*POP2(IRREP2)
 20      CONTINUE
 10   CONTINUE
C
C  LOOP OVER IRREPS.
C
      IOFF=1
      ILIST=LSTNUM
      DO 30 IRREPR=1,NIRREP
       IRREPL=DIRPRD(IRREPR,IRREPX)
C
C GET NUMBER OF DISTRIBUTIONS AND DISTRIBUTION SIZE.
C
       NUMDIS=NUMDSS(IRREPR)
       DISSIZ=IRPDPD(IRREPL,ISYTYP(1,ILIST))
C
C PLUNK VECTOR ON MOINTS.
C
       IF((NUMDIS*DISSIZ).EQ.0)GOTO 30
       CALL PUTLST(Z(IOFF),1,NUMDIS,1,IRREPR,LSTNUM)
       IOFF=IOFF+NUMDIS*DISSIZ
30    CONTINUE
      RETURN
      END