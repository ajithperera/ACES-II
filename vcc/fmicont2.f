      SUBROUTINE FMICONT2(ICORE,MAXCOR,IUHF)
C
C     USED for PCCSD MODIFICATIONS
C
C     ARGUMENTS IN CALL ICORE,MAXCOR,IUHF
C
C THIS PROGRAM CALCULATES THE TERM
C
C  P(IJ) SUM M T(IM,EF) F(MJ)
C
C SYMMETRY ADAPTED 
C 
C
C IN RHF :
C
C - SUM m T(Im,Ab) F(mj) + SUM M T(Mj,Ab) F(MJ)
C
C IN UHF
C
C - SUM M T(IM,AB) F(MJ) + SUM M T(MJ,AB) F(MI)
C
C - SUM m T(Im,Ab) F(mj) + SUM M T(Mj,Ab) F(MI)
C
C - SUM m T(im,ab) F(mj) + SUM m T(mj,ab) F(mi)
C
C 
C FOR QCISD AND CCSD IN ADDITION THE TERMS
C
C - SUM M T(M,A) F(MI) 
C
C - SUM m T(m,a) F(mi)   (UHF only)
C
C ARE CALCULATED
C
C FOR CCSD METHODS THE FOLLOWING TERM
C
C 1/2 SUM E F(E,M) T(E,I)
C
C HAS TO BE ADDED TO THE F(MI) BEFORE THE
C CONTRACTION WITH THE T2 AMPLITUDES IS PERFORMED
C
CEND
C
C CODED JG JUNE/90
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYT,DISSYZ,POP,VRT
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,QCISD,UCC
      LOGICAL LINCC,CICALC,ROHF4,ITRFLG
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWO  
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,
     &             NF1AA,NF1BB,NF2AA,NF2BB
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /FLAGS/IFLAGS(100)
      COMMON /LINEAR/ LINCC,CICALC
      COMMON /SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC
      COMMON /ROHF/ ROHF4,ITRFLG
C
      EQUIVALENCE(METHOD,IFLAGS(2))
C
      DATA ONE,ONEM,HALF/1.0D0,-1.0D0,0.5D0/
C
C   CALCULATE SIZE OF F(M,I) ARRAY
C
      NFAA=NF1AA
      NFBB=NF1BB
      I0AA=MAXCOR+1-NFAA*IINTFP 
      MXCOR=MAXCOR-NFAA*IINTFP
      IF(IUHF.EQ.0) THEN
       I0BB=I0AA
      ELSE
       I0BB=I0AA-NFBB*IINTFP
       MXCOR=MXCOR-NFBB*IINTFP
      ENDIF
      CALL GETLST(ICORE(I0AA),1,1,1,1,91)
      IF(IUHF.EQ.1) THEN
       CALL GETLST(ICORE(I0BB),1,1,1,2,91)
      ENDIF
C
C     AA AND BB SPIN CASES
C
      IF(IUHF.EQ.1) THEN
C
C
       DO 100 ISPIN=1,2
C
        IF(ISPIN.EQ.1) THEN
         I000=I0AA
        ELSE
         I000=I0BB
        ENDIF
        LISTT=ISPIN+43
        LISTZ=ISPIN+60
C
        DO 50 IRREP=1,NIRREP 
C
        NOCCSQ=0
        DO 45 IRREPJ=1,NIRREP
         IRREPI=DIRPRD(IRREP,IRREPJ)
         NOCCSQ=NOCCSQ+POP(IRREPJ,ISPIN)*POP(IRREPI,ISPIN)
45      CONTINUE

C
C       RETRIEVE T2 AMPLITUDES AND CALCLUATE NEW ONES
C
         DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT)) 
         DISSYZ=IRPDPD(IRREP,ISYTYP(1,LISTZ))
         NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT))
         NUMSYZ=IRPDPD(IRREP,ISYTYP(2,LISTZ))
         I001=1
         xI002=I001+IINTFP*(1.*NOCCSQ)*DISSYT
         I002=I001+IINTFP*NOCCSQ*DISSYT
         xI003=xI002+IINTFP*(1.*NOCCSQ)*DISSYZ
         I003=I002+IINTFP*NOCCSQ*DISSYZ
         IF(MIN(NUMSYT,NUMSYZ,DISSYT,DISSYZ).NE.0) THEN
          I004=I003+IINTFP*MAX(DISSYT,DISSYZ)
          xI004=xI003+IINTFP*MAX(DISSYT,DISSYZ)
          IF(I004.LT.MXCOR.and. xI004.LT.MXCOR) THEN
C
C
C    IN CORE VERSION
C
          CALL FMIAA1(ICORE(I001),ICORE(I002),ICORE(I000),
     &                POP(1,ISPIN),NOCCSQ,DISSYT,DISSYZ,NUMSYT,
     &                NUMSYZ,NFAA,LISTT,LISTZ,IRREP,ICORE(I003))
          ELSE
           STOP 'FMIAA1'
          ENDIF
         ENDIF
50      CONTINUE
100    CONTINUE
      ENDIF
  
C
C      AB SPIN CASE
C
      LISTT=46
      LISTZ=63
C
C   LOOP OVER IRREPS
C
      DO 200 IRREP=1,NIRREP
C
C   RETRIEVE T2 AMPLITUDES AND CALCULATE Z-AMPLITUDES
C
       DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
       DISSYZ=IRPDPD(IRREP,ISYTYP(1,LISTZ))
       NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT))
       NUMSYZ=IRPDPD(IRREP,ISYTYP(2,LISTZ))
       I001=1
       I002=I001+IINTFP*NUMSYT*DISSYT
       xI002=I001+IINTFP*(1.*NUMSYT)*DISSYT
       I003=I002+IINTFP*NUMSYZ*DISSYZ
       xI003=xI002+IINTFP*(1.*NUMSYZ)*DISSYZ
       IF(MIN(NUMSYT,NUMSYZ,DISSYT,DISSYZ).NE.0) THEN
        I004=I003+IINTFP*MAX(DISSYT,DISSYZ,NUMSYT,NUMSYZ)*3
        xI004=xI003+IINTFP*MAX(DISSYT,DISSYZ,NUMSYT,NUMSYZ)*3
        IF(I004.LT.MXCOR.and.xI004.LT.MXCOR) THEN
C
C       IN CORE ALGORITHM
C
         CALL FMIAB1(ICORE(I001),ICORE(I002),ICORE(I0AA),ICORE(I0BB),
     &               POP(1,1),POP(1,2),VRT(1,1),VRT(1,2),DISSYT,
     &               DISSYZ,NUMSYT,NUMSYZ,NFAA,NFBB,LISTT,LISTZ,
     &               IRREP,IUHF,ICORE(I003))
        ELSE
C
C       OUT CORE ALGORITHM
C
         STOP 'FMIAB1'
        ENDIF
       ENDIF
200   CONTINUE
      RETURN
      END
