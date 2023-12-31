      SUBROUTINE INIT201(ICORE,MAXCOR,IUHF)
C
C INITIALIZES DOUBLE AMPLITUDE FOR (0,1) FS SECTOR WITH CORRESPONDING
C DENOMINATOR-WEIGHTED HBAR ELEMENTS
C
CEND
      IMPLICIT INTEGER (A-Z)
      double precision snrm2
      LOGICAL RHF
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
C
      RHF=IUHF.EQ.0
C
C RHF REFERENCE CASE.  
C
      IF(RHF)THEN
       LISTW=110
       LISTT=199
       DO 10 IRREP=1,NIRREP
        NUMDST=IRPDPD (IRREP,ISYTYP(2,LISTT))
        DISSYT=IRPDPD (IRREP,ISYTYP(1,LISTT))
        I000=1
        I010=I000+IINTFP*NUMDST*DISSYT
        CALL GETLST(ICORE(I000),1,NUMDST,1,IRREP,LISTW)
        CALL PUTLST(ICORE(I000),1,NUMDST,1,IRREP,LISTT)
10     CONTINUE
C
      ELSE
C
C UHF REFERENCE CASE.  
C
       DO 20 ISPIN=1,4
        LISTW=106+ISPIN
        LISTT=195+ISPIN
        DO 21 IRREP=1,NIRREP
         NUMDST=IRPDPD(IRREP,ISYTYP(2,LISTT))
         DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
         I000=1
         I010=I000+IINTFP*NUMDST*DISSYT
         CALL GETLST(ICORE(I000),1,NUMDST,1,IRREP,LISTW)
         CALL PUTLST(ICORE(I000),1,NUMDST,1,IRREP,LISTT)
21      CONTINUE
20     CONTINUE
C
      ENDIF
C
      RETURN
      END
