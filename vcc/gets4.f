      SUBROUTINE GETS4(ICORE,MAXCOR,E4S,IUHF)
C
C ROUTINE COMPUTES THE FOURTH-ORDER LINEAR SINGLES ENERGY FROM THE RELATION
C
C            E4S = |T(2)**2|*D(ia)
C
C AND IS USED IN ROHF-MBPT CALCULATIONS
C
C HERE, THE T1 INCREMENTS *BEFORE DENOMINATOR WEIGHTING* ARE PICKED UP
C  AND USED, SO THAT THE ACTUAL WORKING EQUATION IS 
C
C            E4S = |TINC(2)**2|/D(ia).
C
C  WHEN THINGS ARE DONE THIS WAY, THEN THE RECIPROCAL DENOMINATORS
C   STORED ON LISTS 48-50 CAN BE USED DIRECTLY WITH VECPRD.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(MAXCOR)
      DOUBLE PRECISION E4S,FACTOR,EAA,E4X,SDOT
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF(2,2)
      EAA=0.0
      E4S=0.0
      FACTOR=1.0
      IF(IUHF.EQ.0)FACTOR=2.0
C
C COMPUTE AA AND BB (BB FOR UHF ONLY) CONTRIBUTIONS TO E4S.
C
      DO 100 ISPIN=1,1+IUHF
       NSIZE=NT(ISPIN)
       I000=1
       I010=I000+IINTFP*NSIZE
       I020=I010+IINTFP*NSIZE
       CALL GETLST(ICORE(I000),1,1,1,2+ISPIN,90)
       CALL VECPRD(ICORE(I000),ICORE(I000),ICORE(I000),NSIZE)
       CALL GETLST(ICORE(I010),1,1,1,9,63+ISPIN)
       E4X=SDOT(NSIZE,ICORE(I000),1,ICORE(I010),1)
       EAA=EAA+E4X
100   CONTINUE
      E4S=E4S+FACTOR*EAA
      RETURN
      END
