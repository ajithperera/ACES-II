C
C ***************************************************************
C  PROCEDURES THAT SKIP BETWEEN S-VECTORS AND S-LISTS
C ***************************************************************
C
      SUBROUTINE LOADS(S, NEA, ISPIN, IUHF, LISTS1, LISTS2, SPINAD,
     $   ICORE, MAXCOR)
C
C  THIS ROUTINE LOADS THE COMPLETE S-VECTOR INTO CORE, ORDERED AS 
C  FOLLOWS: 
C
C    SINGLES(ISPIN), DOUBLES (ISPIN, ISPIN), DOUBLES(ISPIN, 3-ISPIN)  [ UHF ]
C
C    SINGLES(ISPIN), DOUBLES (ISPIN, 3-ISPIN)                         [ RHF ]
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION S
      LOGICAL SPINAD
      DIMENSION S(NEA), LISTS2(2,2), ICORE(MAXCOR)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SINFO/ NS(8), SIRREP
C
      CALL GETLST(S(1), 1, 1, 1, ISPIN, LISTS1)
      IOFF = 1 + VRT(SIRREP, ISPIN) 
      IF (IUHF .NE. 0) THEN
         CALL GETLEN2(LENAA, IRPDPD(1, ISYTYP(1, LISTS2(ISPIN, ISPIN))),
     $      POP(1,ISPIN),NS)
         CALL GETALLS2(S(IOFF),LENAA, POP(1,ISPIN), NS, 1,
     $      LISTS2(ISPIN, ISPIN))
         IOFF = IOFF + LENAA
      ENDIF
      MSPIN = 3 - ISPIN
      CALL GETLEN2(LENAB, IRPDPD(1, ISYTYP(1, LISTS2(ISPIN, MSPIN))),
     $      POP(1,MSPIN+IUHF-1),NS)
      CALL GETALLS2(S(IOFF),LENAB, POP(1,MSPIN+IUHF-1), NS, 1,
     $      LISTS2(ISPIN, MSPIN))
C
      IF (IUHF .EQ. 0) THEN
         IF (SPINAD) THEN
            CALL SPNTSING(NEA, S, ICORE, MAXCOR, SIRREP, 2)
         ENDIF
      ENDIF
C
      RETURN
      END