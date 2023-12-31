      SUBROUTINE GET3(X,LISTX,IRREPX,GETIND,IRREPIND,
     &                POPP,POPQ,POPR,POPS,TYPEPQ,
     &                TYPERS,GETTYP,SYMEX,SPINAD,SCR)
C
C THIS ROUTINE READS THE LIST X(pq,rs) ON DISK AND RETURNS 
C X(pq;r--s) [ALL pqr FOR A GIVEN s].  THE INPUT IS THE IRREP 
C OF THE SPECIFIC s AND ITS SYMMETRY.
C
C THE ORDERING OF SYMMETRIES IN THE OUTPUT VECTOR IS AS FOLLOWS
C
C       X(p q r)
C
C         GX 1 1
C         GX 2 1
C         . . .
C         GX 1 2
C         GX 2 2
C         . . .
C
C GX IS DEFINED BY THE RESTRICTION 
C
C     GAMMA(P)*GAMMA(Q)*GAMMA(R)*GAMMA(S)=GAMMA(X)
C
CEND
      IMPLICIT INTEGER (A-Z)
      CHARACTER*3 GETTYP
      LOGICAL SPINAD,SYMEX
      DIMENSION X(*),POPR(*),POPS(*),SCR(*)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFEA(2),NFMI(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYMLOC/ISYMOFF(8,8,25)
      COMMON/SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
C
      IOFFB=1
C
      IF(GETTYP.EQ.'123')THEN
C
C COLLECT X(pqr) FOR s=GETIND, IRREP(s)=IRREPIND
C
       IRREPS=IRREPIND
       S=GETIND
       DO 10 IRREPR=1,NIRREP
        IRREPRS=DIRPRD(IRREPR,IRREPS)
        IRREPPQ=DIRPRD(IRREPX,IRREPRS)
        DISSYX=IRPDPD(IRREPPQ,ISYTYP(1,LISTX))
        IF(SYMEX)THEN
         DISSYX2=IRPDPD(IRREPPQ,TYPEPQ)
        ELSE
         DISSYX2=DISSYX
        ENDIF
        NUMDIS=IRPDPD(IRREPRS,ISYTYP(2,LISTX))
        NUMR=POPR(IRREPR)
C
C CALCULATE DISTRIBUTION OFFSET TO THE FIRST VALUE OF r FOR THIS s
C
        IDIS0=ISYMOFF(IRREPS,IRREPRS,TYPERS)+(S-1)*NUMR
        CALL GETLST(X(IOFFB),IDIS0,NUMR,1,IRREPRS,LISTX)
        IF(SPINAD)THEN
         MAXX=MAX(NUMR,DISSYX)
         ITMP1=1
         ITMP2=ITMP1+IINTFP*MAXX
         CALL SPINAD3(IRREPPQ,POPP,DISSYX,NUMR,X(IOFFB),
     &                SCR(ITMP1),SCR(ITMP2))
        ENDIF
        IF(SYMEX)THEN
         CALL SYMEXP2(IRREPPQ,POPP,DISSYX2,DISSYX,NUMR,X(IOFFB),
     &                X(IOFFB))
        ENDIF
        IOFFB=IOFFB+IINTFP*DISSYX2*NUMR
10     CONTINUE
C
      ELSEIF(GETTYP.EQ.'124')THEN
C
C COLLECT X(pqs) FOR r=GETIND, IRREP(r)=IRREPIND
C
       IRREPR=IRREPIND
       R=GETIND
       DO 20 IRREPS=1,NIRREP
        IRREPRS=DIRPRD(IRREPR,IRREPS)
        IRREPPQ=DIRPRD(IRREPX,IRREPRS)
        DISSYX=IRPDPD(IRREPPQ,ISYTYP(1,LISTX))
        IF(SYMEX)THEN
         DISSYX2=IRPDPD(IRREPPQ,TYPEPQ)
        ELSE
         DISSYX2=DISSYX
        ENDIF
        NUMDIS=IRPDPD(IRREPRS,ISYTYP(2,LISTX))
        NUMR=POPR(IRREPR)
        NUMS=POPS(IRREPS)
        DO 11 S=1,NUMS
C
C CALCULATE DISTRIBUTION OFFSET TO THIS rs DISTRIBUTION
C
         IDIS0=ISYMOFF(IRREPS,IRREPRS,TYPERS)+(S-1)*NUMR+R-1
         CALL GETLST(X(IOFFB),IDIS0,1,1,IRREPRS,LISTX)
         IF(SPINAD)THEN
          MAXX=MAX(DISSYX,1)
          ITMP1=1
          ITMP2=ITMP1+IINTFP*MAXX
          CALL SPINAD3(IRREPPQ,POPP,DISSYX,1,X(IOFFB),
     &                 SCR(ITMP1),SCR(ITMP2))
         ENDIF
         IF(SYMEX)THEN
          CALL SYMEXP2(IRREPPQ,POPP,DISSYX2,DISSYX,NUMR,X(IOFFB),
     &                 X(IOFFB))
         ENDIF
         IOFFB=IOFFB+IINTFP*DISSYX2
11      CONTINUE
20     CONTINUE
C
      ENDIF
C
      RETURN
      END
