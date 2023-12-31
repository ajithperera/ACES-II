      SUBROUTINE INIT101(ICORE,MAXCOR,IUHF)
C
C INITIALIZES SINGLE AMPLITUDE FOR (0,1) FS SECTOR WITH CORRESPONDING
C HBAR ELEMENTS
C
CEND
      IMPLICIT INTEGER (A-Z)
      LOGICAL RHF
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
C
      RHF=IUHF.EQ.0
C
C RHF REFERENCE CASE.  SPIN ADAPT AMPLITUDES.
C
      IF(.NOT.RHF)THEN
       CALL GETLST(ICORE,1,1,1,2,91)
       CALL PUTLST(ICORE,1,1,1,4,94)
      ENDIF
      CALL GETLST(ICORE,1,1,1,1,91)
      CALL PUTLST(ICORE,1,1,1,3,94)
C
      RETURN
      END
