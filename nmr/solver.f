
      SUBROUTINE SOLVER(IRREP,NPERT,BAI,ICORE,MAXCOR,IUHF,ANTI)
C
C  THIS ROUTINE CONTROLS THE FORMATION OF THE A MATRIX
C  AND CALLS THE ACTUAL SOLVER FOR THE FIRST-ORDER
C  Z-VECTOR EQUATION
C
CEND
C
C  CODED SEPTEMBER/91 JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL INCORE,ROHF,SEMI,QRHF,ANTI
      INTEGER DIRPRD,POP,VRT
      DIMENSION BAI(1),ICORE(MAXCOR)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/FLAGS/ IFLAGS(100)
      COMMON/REF/ROHF,QRHF,SEMI
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
C
      DATA TOL/0.0D+0/
C
      CALL TIMER(1)
C
C  SET PARAMETERS FOR SOLVING THE FIRST-ORDER Z-VECTOR-EQUATION
C
      KMAX=IFLAGS(31)
      CONV=10.D0**(-IFLAGS(30))
C
C  ALLOCATE CORE MEMORY FOR THE ORBITAL ENERGIES
C
      MXCOR=MAXCOR
      NORBA=NOCCO(1)+NVRTO(1)
      IEVAA=MXCOR+1-IINTFP*NORBA
      MXCOR=MXCOR-IINTFP*NORBA
      CALL GETREC(20,'JOBARC','SCFEVALA',IINTFP*NORBA,ICORE(IEVAA))
      IF(IUHF.EQ.1) THEN
       NORBB=NOCCO(2)+NVRTO(2)
       IEVBB=IEVAA-IINTFP*NORBB
       MXCOR=MXCOR-IINTFP*NORBB
       CALL GETREC(20,'JOBARC','SCFEVALB',IINTFP*NORBB,ICORE(IEVBB))
      ENDIF
C
C   FORM NOW THE A MATRIX AND SOLVE FOR d DIA/ d chi *(E(I)-E(A))
C
      IF(IUHF.EQ.0) THEN
C
       NUMSYW=IRPDPD(IRREP,ISYTYP(1,19))
C I010 ... A MATRIX OF DIMENSION NUMSYW*NUMSYW
       I010=1
C I020 ... SCRATCH VECTOR FOR CONSTRUCTING THE A MATRIX
       I020=I010+IINTFP*NUMSYW*NUMSYW
C I030 ... TOTAL LENGTH OF CORE REQUIRED
       I030=I020+IINTFP*NUMSYW*NUMSYW
       IF(I030.LT.MXCOR) THEN
       IF(ANTI) THEN
C
C A-MATRIX FOR COMPLEX PERTURBATIONS
C
        CALL MKARHF3(IRREP,ICORE(I010),ICORE(I020))
C
       ELSE
C
C A-MATRIX FOR REAL PERTURBATIONS
C
        CALL MKARHF2(IRREP,ICORE(I010),ICORE(I020))
C
       ENDIF
       ELSE
        CALL INSMEM('MKARHF',I030,MXCOR)
       ENDIF
C
C  SCALE BAI WITH THE ORBITAL ENERGY DENOMINATOR
C
       CALL FORMU(IRREP,NPERT,NUMSYW,BAI,ICORE(IEVAA),POP(1,1),
     &            VRT(1,1),NOCCO(1))
C
       CALL TIMER(1)
       write(6,6001) TIMENEW
6001   FORMAT(' Construction of A-matrix required ',f5.1,' seconds.')
C
C       ALLOCATE MEMORY FOR SOLVING THE CPHF-EQUATION
C
C    I010  .... A-MATRIX (ALREADY ALLOCATED
C
C    I020  .... SCRATCH VECTOR FOR BAINEW, SOLVING THE LINEAR EQUATION
C
      I030=I020+IINTFP*MAX(NUMSYW,2*KMAX)*NPERT
C
C    I030  .... THE VECTOR WHICH SPAN THE ITERATIVE SUBSPACE
C
      I040=I030+IINTFP*(KMAX+1)*NUMSYW*NPERT
C
C    I040  .... THE SMALL A MATRIX
C
      I050=I040+IINTFP*KMAX*KMAX*NPERT
C
C    I050  .... THE NORM OF ALL EXPANSION VECTORS
C
      I060=I050+IINTFP*(KMAX+1)*NPERT
C
C    I060  .... A VECTOR HOLDING NPERT SCALING FACTORS
C
      I070=I060+IINTFP*NPERT
C
C    I070  .... A SCRATCH VECTOR OF LENGTH NPERT
C
      I080=I070+NPERT
C
      IF(I080.LT.MXCOR) THEN
       CALL LINEQ1(IRREP,NPERT,ICORE(I010),BAI,ICORE(I020),
     &             ICORE(I030),ICORE(I040),ICORE(I050),
     &             ICORE(I060),ICORE(I070),
     &             ICORE(IEVAA),CONV,NUMSYW,KMAX,NOCCO(1))
      ELSE
       CALL INSMEM('LINEQ1',I080,MXCOR)
      ENDIF
C
C   UHF CASE
C
       ELSE IF(IUHF.EQ.1) THEN
C
C CHECK IF ROHF
C
        IF(ROHF) THEN
C
C  ROHF RUN
C
         NUMSYWA=IRPDPD(IRREP,ISYTYP(2,19))
         NUMSYWB=IRPDPD(IRREP,ISYTYP(2,20))
C
C I010 .... 3 BLOCKS OF A MATRICES
C
         I010=1
C
C  I020 ... SCRATCH ARRAY FOR CONTRUCTING THE A-MATRIX
C
         I020=I010+IINTFP*(NUMSYWA*NUMSYWA+
     &               NUMSYWB*NUMSYWB+NUMSYWA*NUMSYWB)
C
C I031 ... FOCK MATRICES (OCCUPIED-OCCUPIED BLOCK)
C
         I031=I020+IINTFP*MAX(NUMSYWA*NUMSYWA,NUMSYWB*NUMSYWB)
C
C I041 ... FOCK MATRICES (VIRTUAL-VIRTUAL BLOCK)
C
         I041=I031+IINTFP*(NF1(1)+NF1(2))
C
C I051 ... FOCK MATRIX (OCCUPIED-VIRTUAL BETA SPIN CASE)
C
         I051=I041+IINTFP*(NF2(1)+NF2(2))
C
C IEND1 ... TOTAL LENGTH REQUIRED FOR A-MATRIX CONSTRUCTION
C
         IEND1=I051+IINTFP*NT(2)
C
C  ALLOCATE CORE FOR SOLUTION OF THE CPHF-EQUATION
C
C  I010 .... A MATRIX ALREADY ALLOCATED
C
C  I020 .... SCRATCH VECTOR FOR BAINEW, SOLVING THE LINEAR EQUATION
C
         I025=I020+IINTFP*NUMSYWA*NPERT
         I030=I020+IINTFP*MAX(NUMSYWA+NUMSYWB,2*KMAX)*NPERT
         I035=I030+IINTFP*(KMAX+1)*NUMSYWA*NPERT
C
C  I030, I035 ... THE VECTOR WHICH SPAN THE ITERATIVE SUBSPACE
C
         I040=I035+IINTFP*(KMAX+1)*NUMSYWB*NPERT
C
C     I040  ... THE SMALL A MATRIX
C
         I050=I040+IINTFP*KMAX*KMAX*NPERT
C
C     I050 .... THE NORM OF ALL EXPANSION VECTORS
C
         I060=I050+IINTFP*(KMAX+1)*NPERT
C
C    I060  .... A VECTOR HOLDING NPERT SCALING FACTORS
C
         I070=I060+IINTFP*NPERT
C
C    I070  .... A SCRATCH VECTOR OF LENGTH NPERT
C
         I080=I070+NPERT
C
      IEND2=I080+NPERT
C
      IF(MAX(IEND1,IEND2).LT.MAXCOR) THEN
       INCORE=.TRUE.
       CALL MKAROHF(IRREP,ICORE(I010),ICORE(I020),ICORE(I041),
     &              ICORE(I031),ICORE(I051),NF2(1),NF1(1),NT(2))
      ELSE
       INCORE=.FALSE.
       CALL ERREX
      ENDIF
C
      CALL ROHFU(IRREP,NPERT,NUMSYWA,NUMSYWB,BAI,
     &           BAI(1+NUMSYWA*NPERT),
     &           ICORE(IEVAA),ICORE(IEVBB),POP(1,1),POP(1,2),
     &           VRT(1,1),VRT(1,2),NOCCO(1),NOCCO(2))
C
      CALL LINEQ2(IRREP,NPERT,ICORE(I010),BAI,BAI(1+NUMSYWA*NPERT),
     &            ICORE(I020),ICORE(I025),ICORE(I030),
     &            ICORE(I035),ICORE(I040),ICORE(I050),ICORE(I060),
     &            ICORE(I070),ICORE(IEVAA),ICORE(IEVBB),CONV,NUMSYWA,
     &            NUMSYWB,KMAX,NOCCO(1),NOCCO(2),INCORE,.TRUE.)
C
      ELSE
C
C  UHF REQUIRES MUCH MORE CORE MEMORY THAN RHF, THEREFORE CONSIDER
C  IN CORE AND OUT CORE VERSION
C
C  IF INCORE POSSIBLE, LOAD THE INTEGRALS IN THE APPROBIATE ALLOCATION
C  AND RUN THE CALCULATION
C  OTHERWISE READ THEM IN WHEN THE A*Y PRODUCT IS FORMED
C
C  TRY NOW TO ALLOCATE FOR IN-CORE
C
      NUMSYWA=IRPDPD(IRREP,ISYTYP(2,19))
      NUMSYWB=IRPDPD(IRREP,ISYTYP(2,20))
C
C  I010 ... 3 BLOCKS OF THE A-MATRIX
C
      I010=1
C
C  I020 ... SCRATCH ARRAY FOR CONTRUCTING THE A-MATRIX
C
      I020=I010+IINTFP*(NUMSYWA*NUMSYWA+NUMSYWB*NUMSYWB+NUMSYWA*NUMSYWB)
C
C  IEND1... TOTAL LENGTH OF CORE REQUIRED FOR A-MATRIX CONSTRUCTION
C
      IEND1=I020+IINTFP*MAX(NUMSYWA*NUMSYWA,NUMSYWB*NUMSYWB)
C
C  ALLOCATE CORE FOR SOLUTION OF THE CPHF-EQUATION
C
C  I010 .... A MATRIX ALREADY ALLOCATED
C
C  I020 .... SCRATCH VECTOR FOR BAINEW, SOLVING THE LINEAR EQUATION
C
      I025=I020+IINTFP*NUMSYWA*NPERT
      I030=I020+IINTFP*MAX(NUMSYWA+NUMSYWB,2*KMAX)*NPERT
      I035=I030+IINTFP*(KMAX+1)*NUMSYWA*NPERT
C
C  I030, I035 ... THE VECTOR WHICH SPAN THE ITERATIVE SUBSPACE
C
      I040=I035+IINTFP*(KMAX+1)*NUMSYWB*NPERT
C
C     I040  ... THE SMALL A MATRIX
C
      I050=I040+IINTFP*KMAX*KMAX*NPERT
C
C     I050 .... THE NORM OF ALL EXPANSION VECTORS
C
      I060=I050+IINTFP*(KMAX+1)*NPERT
C
C    I060  .... A VECTOR HOLDING NPERT SCALING FACTORS
C
      I070=I060+IINTFP*NPERT
C
C    I070  .... A SCRATCH VECTOR OF LENGTH NPERT
C
      I080=I070+NPERT
C
      IEND2=I080+NPERT
C
      IF(MAX(IEND1,IEND2).LT.MXCOR) THEN
       INCORE=.TRUE.
       CALL MKAUHF2(IRREP,ICORE(I010),ICORE(I020))
      ELSE
       INCORE=.FALSE.
C
C   ALLOCATE MEMORY FOR OUT-CORE STEP
C
C   A-MATRIX BUFFER
C
       I020=I010+IINTFP*2*MAX(NUMSYWA*NUMSYWA,NUMSYWB*NUMSYWB)
C
       I025=I020+IINTFP*NUMSYWA*NPERT
       I030=I020+IINTFP*MAX(NUMSYWA+NUMSYWB,2*KMAX)*NPERT
       I035=I030+IINTFP*(KMAX+1)*NUMSYWA*NPERT
C
       I040=I035+IINTFP*(KMAX+1)*NUMSYWB*NPERT
C
       I050=I040+IINTFP*KMAX*KMAX*NPERT
C
      I060=I050+IINTFP*(KMAX+1)*NPERT
C
C    I060  .... A VECTOR HOLDING NPERT SCALING FACTORS
C
      I070=I060+IINTFP*NPERT
C
C    I070  .... A SCRATCH VECTOR OF LENGTH NPERT
C
      I080=I070+NPERT
C
       IEND=I080+NPERT
C
       IF(IEND.GE.MXCOR) CALL INSMEM('LINEQ2',IEND,MXCOR)
C
      ENDIF
C
C  SCALE XIA WITH THE ORBITAL ENERGY DENOMINATORS
C
      CALL FORMU(IRREP,NPERT,NUMSYWA,BAI,ICORE(IEVAA),
     &           POP(1,1),VRT(1,1),NOCCO(1))
      CALL FORMU(IRREP,NPERT,NUMSYWB,BAI(1+NUMSYWA*NPERT),
     &           ICORE(IEVBB),POP(1,2),VRT(1,2),NOCCO(2))
C
      CALL LINEQ2(IRREP,NPERT,ICORE(I010),BAI,BAI(1+NUMSYWA*NPERT),
     &            ICORE(I020),ICORE(I025),ICORE(I030),
     &            ICORE(I035),ICORE(I040),ICORE(I050),ICORE(I060),
     &            ICORE(I070),ICORE(IEVAA),ICORE(IEVBB),CONV,NUMSYWA,
     &            NUMSYWB,KMAX,NOCCO(1),NOCCO(2),INCORE,.FALSE.)
      ENDIF
      ENDIF
C
      CALL TIMER(1)
C
      write(6,6002) TIMENEW
6002  FORMAT(' Iterative solution of the linear equations',
     &        ' required ',f5.1,' seconds.')
      RETURN
      END