      SUBROUTINE FEACONT(ICORE,MAXCOR,IUHF)
C
C THIS ROUTINE CALCULATES THE TERM
C
C  P(AB) SUM E T(IJ,AE) F(E,B)
C
C USING SYMMETRY PACKED ARRAYS
C 
C
C IN RHF :
C
C SUM e T(Ij,Ae) F(eb) - SUM E  T(Ij,Eb) F(EA)
C
C IN UHF
C
C SUM E T(IJ,AE) F(EB) - SUM E T(IJ,EB) F(EA)
C
C SUM e T(Ij,Ae) F(eb) - SUM E T(Ij,Eb) F(EA)
C
C SUM e T(ij,ae) F(eb) - SUM e T(ij,eb) F(ea)
C
C
C FOR CCSD AND QCISD IN ADDITION THE TERMS
C
C SUM T(I,E) F(EA)
C
C SUM T(i,e) F(ea)
C
C ARE CALCULATED
C
C FURTHER FOR CCSD METHODS THE TERM
C
C - 1/2 SUM M F(E,M) T(A,M) 
C
C IS ADDED TO THE FEA CONTRIBUTION BEFORE THE
C CONTRACTION WITH THE T2 AMPLITUDES IS CARRIED OUT.
C
CEND
C
C
C CODED JG JUNE/90 JULY/90
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LINCC,CICALC,ROHF4,ITRFLG
      INTEGER DIRPRD,DISSYT,DISSYZ,POP,VRT
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,QCISD,UCC
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWO
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,
     &             NF1AA,NF1BB,NF2AA,NF2BB
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /LINEAR/ LINCC,CICALC
      COMMON /SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC
      COMMON /ROHF/ ROHF4,ITRFLG
C
      EQUIVALENCE(IFLAGS(2),METHOD)
C
      DATA ONE,HALFM /1.0D0,-0.5D0/
C
C    CALCULATE SIZE OF F(A,E) ARRAY
C    AND GET THESE ARRAYS
C
      NFAA=NF2AA
      NFBB=NF2BB
      I0AA=MAXCOR+1-NFAA*IINTFP
      MXCOR=MAXCOR-NFAA*IINTFP
      IF(IUHF.EQ.0) THEN
       I0BB=I0AA
      ELSE
       I0BB=I0AA-NFBB*IINTFP
       MXCOR=MXCOR-NFBB*IINTFP
      ENDIF
      CALL GETLST(ICORE(I0AA),1,1,1,1,92)
      IF(IUHF.EQ.1) THEN
       CALL GETLST(ICORE(I0BB),1,1,1,2,92)
      ENDIF
C
C     skip the crap and go directly to the f*t2 contraction if this is t4
C
      GOTO 1000
C
C     FOR QCISD AND CCSD METHODS GET THE T1-AMPLITUDES
C
C     CALCULATE ADDTIONAL TERMS FOR QCISD AND CCSD
C
      IF(METHOD.GT.9.AND.SING1.OR.METHOD.EQ.6.AND.SING1.OR.ROHF4)THEN
       I0TA=I0BB-NTAA*IINTFP
       I0ZA=I0TA-NTAA*IINTFP
       MXCOR=MXCOR-2*NTAA*IINTFP
       CALL GETLST(ICORE(I0TA),1,1,1,1,90)
       CALL GETLST(ICORE(I0ZA),1,1,1,3,90)
       IF(IUHF.EQ.0) THEN
        I0TB=I0TA
        I0ZB=I0ZA
       ELSE
        I0TB=I0ZA-NTBB*IINTFP
        I0ZB=I0TB-NTBB*IINTFP
        MXCOR=MXCOR-2*NTBB*IINTFP
        CALL GETLST(ICORE(I0TB),1,1,1,2,90)
        CALL GETLST(ICORE(I0ZB),1,1,1,4,90)
       ENDIF
C
       DO 300 ISPIN=1,IUHF+1
C
        IF(ISPIN.EQ.1) THEN
         IOFFF=I0AA
         IOFFT=I0TA
         IOFFZ=I0ZA
         I0Z=I0ZA
        ELSE
         IOFFF=I0BB
         IOFFT=I0TB
         IOFFZ=I0ZB
         I0Z=I0ZB
        ENDIF
C
        DO 250 IRREP=1,NIRREP
C
         NOCC=POP(IRREP,ISPIN)
         NVRT=VRT(IRREP,ISPIN)
         CALL XGEMM('T','N',NVRT,NOCC,NVRT,ONE,ICORE(IOFFF),
     &              NVRT,ICORE(IOFFT),NVRT,ONE,ICORE(IOFFZ),
     &              NVRT)
         IOFFF=IOFFF+NVRT*NVRT*IINTFP
         IOFFT=IOFFT+NVRT*NOCC*IINTFP
         IOFFZ=IOFFZ+NVRT*NOCC*IINTFP
250     CONTINUE
        CALL PUTLST(ICORE(I0Z),1,1,1,ISPIN+2,90)
300    CONTINUE
C
C  IF CCSD METHODS ADD TO FEA THE FOLLOWING CONTRIBUTION
C
C  - 1/2 SUM M F(E,M) T(A,M)
C
       IF(.NOT.QCISD.AND..NOT.LINCC.OR.ROHF4) THEN
C
C  USE ALLOCATION FOR Z TO STORE F(E,M)
C
        CALL GETLST(ICORE(I0ZA),1,1,1,1,93)
C
        IOFFEA=I0AA
        IOFFEM=I0ZA
        IOFT=I0TA
        DO 20 IRREP=1,NIRREP
         NOCC=POP(IRREP,1)
         NVRT=VRT(IRREP,1)
         CALL XGEMM('N','T',NVRT,NVRT,NOCC,HALFM,ICORE(IOFFEM),NVRT,
     &               ICORE(IOFT),NVRT,ONE,ICORE(IOFFEA),NVRT)
         IOFFEA=IOFFEA+NVRT*NVRT*IINTFP 
         IOFFEM=IOFFEM+NOCC*NVRT*IINTFP
         IOFT=IOFT+NOCC*NVRT*IINTFP
20      CONTINUE
        IF(IUHF.EQ.1) THEN
         CALL GETLST(ICORE(I0ZB),1,1,1,2,93)
C
         IOFFEA=I0BB
         IOFFEM=I0ZB
         IOFT=I0TB
         DO 21 IRREP=1,NIRREP
          NOCC=POP(IRREP,2)
          NVRT=VRT(IRREP,2)
          CALL XGEMM('N','T',NVRT,NVRT,NOCC,HALFM,ICORE(IOFFEM),NVRT,
     &                ICORE(IOFT),NVRT,ONE,ICORE(IOFFEA),NVRT)
          IOFFEA=IOFFEA+NVRT*NVRT*IINTFP 
          IOFFEM=IOFFEM+NOCC*NVRT*IINTFP
          IOFT=IOFT+NOCC*NVRT*IINTFP
21       CONTINUE
        ENDIF
       ENDIF
      ENDIF
C
C     for t4 stuff come here directly
C
 1000 CONTINUE
C
C     AA AND BB SPIN CASES
C
      IF(IUHF.EQ.1) THEN
C
C      THESE CASES ARE ONLY NECCESARY IN THE UHF CASE
C      IN RHF THE AAAA AMPLITUDES ARE CALCULATED FROM
C      THE ABAB AMPLITUDES
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
        NVRTSQ=0
        DO 45 IRREPJ=1,NIRREP
         IRREPI=DIRPRD(IRREP,IRREPJ)
         NVRTSQ=NVRTSQ+VRT(IRREPJ,ISPIN)*VRT(IRREPI,ISPIN)
45      CONTINUE
C
C     RETRIEVE T2 AMPLITUDES AND CALCULATE THE CORRESPONDING CONTRIBUTION
C
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
        DISSYZ=IRPDPD(IRREP,ISYTYP(1,LISTZ))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT))
        NUMSYZ=IRPDPD(IRREP,ISYTYP(2,LISTZ))
        I001=1
        I002=I001+IINTFP*NUMSYT*NVRTSQ
        I003=I002+IINTFP*NUMSYZ*NVRTSQ
        IF(MIN(NUMSYT,NUMSYZ,DISSYT,DISSYZ).NE.0) THEN
         I004=I003+IINTFP*MAX(NUMSYT,NUMSYZ)
         IF(I004.LT.MXCOR) THEN
C
C     IN CORE VERSION
C
         CALL FEAAA1(ICORE(I001),ICORE(I002),ICORE(I002),
     &               ICORE(I001),ICORE(I000),VRT(1,ISPIN),
     &               NVRTSQ,DISSYT,DISSYZ,NUMSYT,
     &               NUMSYZ,NFAA,LISTT,LISTZ,IRREP,ICORE(I003))
        ELSE
         STOP 'FEAAA1'
        ENDIF
       ENDIF
50    CONTINUE
100   CONTINUE
      ENDIF
  
C
C      AB SPIN CASE
C
       LISTT=46
       LISTZ=63
C
C    LOOP OVER IRREPS
C
       DO 200 IRREP=1,NIRREP
C
C    RETRIEVE AMPLITUDES AND CALCULATE CONTRIBUTION TO Z
C
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
        DISSYZ=IRPDPD(IRREP,ISYTYP(1,LISTZ))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT))
        NUMSYZ=IRPDPD(IRREP,ISYTYP(2,LISTZ))
        I001=1
        I002=I001+IINTFP*NUMSYT*DISSYT
        I003=I002+IINTFP*NUMSYZ*DISSYZ
        IF(MIN(NUMSYT,NUMSYZ,DISSYT,DISSYZ).NE.0) THEN
         I004=I003+IINTFP*MAX(NUMSYT,NUMSYZ,DISSYZ,DISSYT)*3
         IF(I004.LT.MXCOR) THEN
C
C       IN CORE ALGORITHM
C
         CALL FEAAB1(ICORE(I001),ICORE(I002),ICORE(I002),ICORE(I001),
     &               ICORE(I0AA),ICORE(I0BB),
     &               POP(1,1),POP(1,2),VRT(1,1),VRT(1,2),DISSYT,
     &               DISSYZ,NUMSYT,NUMSYZ,NFAA,NFBB,
     &               LISTT,LISTZ,IRREP,IUHF,ICORE(I003))
        ELSE
C
C       OUT CORE ALGORITHM
C
        STOP 'FEAAB1'
        ENDIF
       ENDIF
200   CONTINUE
      RETURN
      END