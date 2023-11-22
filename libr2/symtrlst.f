      SUBROUTINE SYMTRLST(ICORE,MAXCOR,IUHF,POPP,POPQ,
     &                    POPR,POPS,IRREPX,
     &                    TYPBRA,TYPBRA2,TYPBRA2F,TYPKET,TYPKET2,
     &                    LISTXIN,LISTXOUT,ISWITCH,SPINAD,SYMEX,ANTI)
C
C THIS ROUTINE CONVERTS THE LIST X(pq,rs) ON DISK TO
C X(rq,ps) [ISWITCH=1] OR X(pr,qs) [ISWITCH=2]
C
C INPUT:
C
C     ICORE   : WORKING AREA.  MUST BE AT LEAST AS BIG AS 2*pqr
C     MAXCOR  : SIZE OF WORKING AREA
C     IUHF    : UHF
C     POPP    : POPULATION VECTOR FOR p
C     POPQ    : POPULATION VECTOR FOR q
C     POPR    : POPULATION VECTOR FOR r
C     POPS    : POPULATION VECTOR FOR s
C     IRREPX  : OVERALL SYMMETRY OF QUANTITY X
C     TYPBRA  : SYMMETRY NUMBER [FOR ISYMOFF] FOR pq [AS STORED ON
C               LISTXIN]
C     TYPBRA2 : SYMMETRY NUMBER FOR BRA INDICES IN TARGET [AS STORED ON
C               LISTXOUT]
C     TYPBRA2F: SYMMETRY NUMBER FOR FULL STORAGE MODE OF BRA
C               INDICES IN TARGET [SAME AS TYPBRA2 UNLESS ANTI=.TRUE.]
C     TYPKET  : SYMMETRY NUMBER [FOR ISYMOFF] FOR rs
C     TYPKET2 : SYMMETRY NUMBER FOR KET INDICES IN TARGET
C     LISTXIN : LIST X IS READ FROM
C     LISTXOUT: LIST TRANSPOSED X IS WRITTEN TO
C     ISWITCH : SEE ABOVE
C     SPINAD  : LHS INDICES OF X ARE TO BE SPIN-ADAPTED
C     SYMEX   : IF LHS INDICES OF X ARE TO BE EXPANDED
C               FROM A TRIANGULAR ARRAY WHEN READ.
C     ANTI    : IF LHS INDICES OF TARGET ARE TO BE ANTISYMMETRIZED
C               INTO A TRIANGULAR ARRAY WHEN WRITTEN.
C
CEND
      IMPLICIT INTEGER (A-Z)
      LOGICAL SPINAD,SYMEX,ANTI
      DIMENSION ICORE(MAXCOR),POPP(*),POPQ(*),POPR(*),POPS(*)
      DIMENSION POPDUM(8)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFEA(2),NFMI(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYMLOC/ISYMOFF(8,8,25)
      COMMON/SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
C
C LOOP OVER SYMMETRIES OF INDEX S
C
      DO 10 IRREPS=1,NIRREP
       NUMS=POPS(IRREPS)
       CALL IZERO(POPDUM,8)
       POPDUM(1)=1
C
C CALCULATE SIZE OF pqr DISTRIBUTIONS FOR THIS SYMMETRY
C
       NUMPQR=0
       MAXSCR=-1
       DO 11 IRREPPQ=1,NIRREP
        IRRTMP=DIRPRD(IRREPPQ,IRREPS)
        IRREPR=DIRPRD(IRREPX,IRRTMP)
        IF(SYMEX)THEN
         NUMPQ =IRPDPD(IRREPPQ,TYPBRA2)
        ELSE 
         NUMPQ =IRPDPD(IRREPPQ,ISYTYP(1,LISTXIN))
        ENDIF
        NUMR  =POPR(IRREPR)             
        NUMPQR=NUMR*NUMPQ+NUMPQR
        MAXSCR=MAX(MAXSCR,NUMPQ)
11     CONTINUE
       I000=1
       I010=I000+IINTFP*NUMPQR
       I020=I010+IINTFP*NUMPQR
c YAU - SSTGEN needs scratch space at I020. I hope this is enough.
       I030=I020+MAXSCR
       IF (I030.GT.MAXCOR) CALL INSMEM('SYMTRLST',I030,MAXCOR)
C
C LOOP OVER S VALUES
C
       DO 12 S=1,NUMS
C
C READ IN ALL pqr FOR THIS S
C
        CALL GET3(ICORE,LISTXIN,IRREPX,S,IRREPS,
     &            POPP,POPQ,POPR,POPS,TYPBRA2,TYPKET,'123',
     &            SYMEX,SPINAD,ICORE(I010))
C
C REORDER AND WRITE BACK OUT  
C
        IRRTMP=DIRPRD(IRREPS,IRREPX)
        IF(ISWITCH.EQ.1)THEN
         CALL SSTGEN(ICORE(I000),ICORE(I010),NUMPQR,POPP,POPQ,POPR,
     &               POPDUM,ICORE(I020),IRRTMP,'3214')
         CALL PUT3(ICORE(I010),LISTXOUT,IRREPX,S,IRREPS,POPP,POPS,
     &             TYPBRA2,TYPBRA2F,TYPKET2,'123',ANTI)
        ELSE
         CALL SSTGEN(ICORE(I000),ICORE(I010),NUMPQR,POPP,POPQ,POPR,
     &               POPDUM,ICORE(I020),IRRTMP,'1324')
         CALL PUT3(ICORE(I010),LISTXOUT,IRREPX,S,IRREPS,POPQ,POPS,
     &             TYPBRA2,TYPBRA2F,TYPKET2,'123',ANTI)
        ENDIF
12     CONTINUE
10    CONTINUE
      RETURN
      END
