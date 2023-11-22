      SUBROUTINE DISCO(ICORE,MAXCOR,IUHF,SPCASE,LIST)
C
C THIS ROUTINE EVALUATES THE DISCONNECTED CONTRIBUTION TO THE
C  LAMBDA(2) EQUATION [DIAGRAMS 2 AND 3 IN FIG. 11 JCP 90, 1752 (1989)]
C  FOR VARIOUS SPIN CASES.  THIS CODE IS CALLED BY THE RING ROUTINE AND
C  IS NEEDED FOR ALL CC GRADIENT METHODS WHICH INCLUDE T1 (QCISD AND CCSD).
C  
C
C       OP(L(I,A),F(J,B))   [SPCASE='AAAA']
C
C       OP[L(i,a),F(j,b)]   [SPCASE='BBBB']
C
C       OP[L(I,A),F(j,b)]   [SPCASE='AABB']
C
C       OP[L(i,a),F(J,B)]   [SPCASE='BBAA']
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,FACT
      CHARACTER*4 SPCASE
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      DATA ONE,ONEM /1.0,-1.0/
      IF(SPCASE.EQ.'AAAA')THEN
       LISTL1=1
       LISTF1=1
       LENL  =NT(1)
       LENF  =NT(1)
       LISTL2=190
       LISTF2=93
       FACT=ONEM
      ELSEIF(SPCASE.EQ.'BBBB')THEN
       LISTL1=2
       LISTF1=2
       LENL  =NT(2)
       LENF  =NT(2)
       LISTL2=190
       LISTF2=93
       FACT=ONEM
      ELSEIF(SPCASE.EQ.'AABB')THEN
       LISTL1=1
       LISTF1=2
       IF(IUHF.EQ.0) LISTF1=1
       LENL  =NT(1)
       LENF  =NT(2)
       LISTL2=190
       LISTF2=93
       FACT=ONE
      ELSEIF(SPCASE.EQ.'BBAA')THEN
C
C MODIFY SO THAT TRANSPOSE IS COMPUTED
C
       LISTL1=1
       LISTF1=2
       LENL  =NT(1)
       LENF  =NT(2)
       LISTL2=93
       LISTF2=190
       FACT=ONE
      ENDIF
      I000=1
      I010=I000+IINTFP*LENL*LENF
      I020=I010+IINTFP*LENL
      I030=I020+IINTFP*LENF
      CALL GETLST(ICORE(I010),1,1,1,LISTL1,LISTL2)
      CALL GETLST(ICORE(I020),1,1,1,LISTF1,LISTF2)
      CALL OUTPRD(ICORE(I010),ICORE(I020),ICORE(I000),LENL,
     &            LENF,FACT)
      I020=I010+IINTFP*LENL*LENF
      CALL GETLST(ICORE(I010),1,LENF,1,1,LIST)
      CALL VADD(ICORE(I000),ICORE(I000),ICORE(I010),LENL*LENF,ONE)
      CALL PUTLST(ICORE(I000),1,LENF,1,1,LIST)
      RETURN
      END