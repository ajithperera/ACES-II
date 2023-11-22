      SUBROUTINE FANCYPUT1(Z,BUF,PUTTYP,PERM,ANTI,AUGTYP,FACTOR,
     &                     NUMP,NUMQ,NUMR,NUMS,
     &                     IFIRST,NUMREC,ICACHE,IRREPGET,LIST,RHF)
C
C WITH THIS SUBROUTINE, CERTAIN GAMES CAN BE PLAYED WITH I/O
C  WHICH ALLOW SOME CONTRACTIONS TO BE CARRIED OUT EASILY WITH
C  MATRIX MULTIPLICATION.  THE PURPOSE OF THE ROUTINE IS TO WRITE
C  A SYMMETRY PACKED FOUR-INDEX LIST [I(PQ,RS) OR I(P<Q,RS) WITH 
C  A GAMMA(PQ)=GAMMA(RS) DPD] FROM A SOURCE VECTOR BELONGING TO 
C  A DIFFERENT DPD.  
C
C  SPECIFICALLY IT ALLOWS THE FOLLOWING:
C
C-------------------------------------------------------------------         
C         IN MEMORY                                    WRITTEN AS 
C          -------                                       -------
C I(PS,RQ), I(PR,SQ), I(PR,QS) OR I(PS,QR)               I(PQ,RS) 
C I(PS,RQ), I(PR,SQ), I(PR,QS) OR I(PS,QR)              I(P<Q,RS) 
C-------------------------------------------------------------------         
C
C THIS ROUTINE ACCESSES ALL IRREPS FOR A GIVEN LIST TO FILL ONE
C  IRREP OF THE ALTERNATIVE DPD LIST, AND THEREFORE SHOULD BE USED
C  ONLY WHEN THE REORDERING CANNOT BE DONE EASILY IN CORE.  THE
C  DISK ACCESS IS SEQUENTIAL WITH NON-SEQUENTIAL MEMORY ADDRESSING.
C
C INPUT:
C       Z - TARGET VECTOR FOR I/O OPERATIONS (AS IN GETLST).
C     BUF - SCRATCH ARRAY.  THE LENGTH OF THIS ARRAY MUST BE
C            AT LEAST THAT OF THE LARGEST LOGICAL RECORD FOR
C            THIS SPECIFIC LIST (MAX(IRPDPD(1..NIRREP,ISYTYP(1,LIST)).
C  PUTTYP - REFERS TO THE MODE OF STORAGE AND THE WAY THAT THE
C            ARRAY WILL BE READ. (CHARACTER*2) 
C
C            'FF' - ARRAY IS STORED ON DISK AS I(PQ,RS)
C            'PF' - ARRAY IS STORED ON DISK AS I(P<Q,RS)
C
C            ("P" REFERS TO PACKED, "F" TO FULL)
C
C    PERM - PERMUTATION TYPE (CHARACTER*4)
C            '1432' - WRITE I(PS,RQ) AS I(PQ,RS)
C            '1423' - WRITE I(PR,SQ) AS I(PQ,RS)
C            '1324' - WRITE I(PR,QS) AS I(PQ,RS)
C            '1342' - WRITE I(PS,QR) AS I(PQ,RS)
C
C    ANTI - ANTISYMMETRIZATION TYPE (CHARACTER*1) [IGNORED FOR PUTTYP 'FF']
C            'Y' - ANTISYMMETRIZE BRA INDICES AND WRITE I(P<Q,RS). 
C            'N' - DO NOT ANTISYMMETRIZE BRA INDICES AND WRITE I(P<Q,RS).
C
C  AUGTYP - HOW THE LIST WILL BE INCREMENTED (CHARACTER*1)
C            'S' - PERFORM SAXPY OPERATION INVOLVING INPUT ELEMENTS
C                   OF SOURCE VECTOR AND CORRESPONDING ELEMENTS ON LIST.
C            'O' - OVERWRITE EXISTING ELEMENTS.
C
C  FACTOR - COEFFICIENT TO USE IN SAXPY OPERATIONS SPECIFIED BY AUGTYP='S'
C           (DOUBLE PRECISION REAL). [IGNORED FOR AUGTYP='O'] 
C    NUMP - POPULATION BY IRREP CORRESPONDING TO FASTEST LEFT INDEX ON DISK.
C    NUMQ - POPULATION BY IRREP CORRESPONDING TO SLOWEST LEFT INDEX ON DISK.
C    NUMR - POPULATION BY IRREP CORRESPONDING TO FASTEST RIGHT INDEX ON DISK.
C    NUMS - POPULATION BY IRREP CORRESPONDING TO SLOWEST RIGHT INDEX ON DISK.
C
C  IFIRST - FIRST LOGICAL RECORD OF *ALTERNATIVE DPD* WHICH IS
C           REQUIRED (LIKE GETLST).
C  NUMREC - THE NUMBER OF LOGICAL RECORDS TO RETRIEVE (LIKE GETLST).
C  ICACHE - THE I/O CACHE TO BE USED (LIKE GETLST).
CIRREPGET - THE IRREP OF THE ALTERNATIVE DPD WHICH IS REQUESTED
C           (LIKE GETLST).
C    LIST - THE MOINTS LIST NUMBER (LIKE GETLST).
C
C ***IMPORTANT*** 
C   1. THIS ROUTINE HAS NOT BEEN EXTENSIVELY DEBUGGED.  DO NOT
C      ASSUME THAT IT IS CORRECT FOR ALL CASES.
C   2. PRESENTLY, THE ARGUMENTS IFIRST AND NUMREC ARE IGNORED, AND
C      THIS ROUTINE WILL RETRIEVE THE FULL LIST FOR A GIVEN IRREP.
C   3. CURRENTLY, THE ROUTINE CAN NOT HANDLE LISTS WITH PACKED KET
C      INDICES.  THIS MAY BE ADDED IF THE NEED ARISES.
C
C SPECIAL ROUTINE FOR SPIN-ADAPTED CODE IN W5ABIN
C
CEND
      IMPLICIT INTEGER (A-Z)
      LOGICAL RHF
      DOUBLE PRECISION Z(1),BUF(1),HALF,FACTOR,FACTOR2
      CHARACTER*2 PUTTYP
      CHARACTER*4 PERM
      CHARACTER*1 AUGTYP,ANTI
      EQUIVALENCE (ONDISK(1),ICOUNT1)
      EQUIVALENCE (ONDISK(2),ICOUNT2)
      EQUIVALENCE (ONDISK(3),ICOUNT3)
      EQUIVALENCE (ONDISK(4),ICOUNT4)
      DIMENSION DISDSK(8,8),DSZDSK(8,8),DISMEM(8,8),DSZMEM(8,8)
      DIMENSION NUMP(8),NUMQ(8),NUMR(8),NUMS(8),PQFULL(8),DSZPCK(8)
      DIMENSION NUMP2(8),NUMQ2(8),NUMR2(8),NUMS2(8),IPTR(4),IPTRI(4)
      DIMENSION IRPTAR(4),ONDISK(4)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPX(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      DATA HALF /0.5/
      INDX(I,J,N)=I+(J-1)*N
C
      FACTOR2=HALF*FACTOR
C
C COMPUTE OFFSETS INTO DISTRIBUTIONS AND OFFSETS INTO DISTRIBUTION
C  MEMBERS.
C
      IF(PERM.EQ.'1432')THEN
       CALL ICOPY(NIRREP,NUMP,1,NUMP2,1)
       CALL ICOPY(NIRREP,NUMQ,1,NUMS2,1)
       CALL ICOPY(NIRREP,NUMR,1,NUMR2,1)
       CALL ICOPY(NIRREP,NUMS,1,NUMQ2,1)
       IPTR(1)=1
       IPTR(2)=4
       IPTR(3)=3
       IPTR(4)=2
       IPTRI(1)=1
       IPTRI(2)=4
       IPTRI(3)=3
       IPTRI(4)=2
      ELSEIF(PERM.EQ.'1423')THEN
       CALL ICOPY(NIRREP,NUMP,1,NUMP2,1)
       CALL ICOPY(NIRREP,NUMQ,1,NUMS2,1)
       CALL ICOPY(NIRREP,NUMR,1,NUMQ2,1)
       CALL ICOPY(NIRREP,NUMS,1,NUMR2,1)
       IPTR(1)=1
       IPTR(2)=4
       IPTR(3)=2
       IPTR(4)=3
       IPTRI(1)=1
       IPTRI(2)=3
       IPTRI(3)=4
       IPTRI(4)=2
      ELSEIF(PERM.EQ.'1324')THEN
       CALL ICOPY(NIRREP,NUMP,1,NUMP2,1)
       CALL ICOPY(NIRREP,NUMQ,1,NUMR2,1)
       CALL ICOPY(NIRREP,NUMR,1,NUMQ2,1)
       CALL ICOPY(NIRREP,NUMS,1,NUMS2,1)
       IPTR(1)=1
       IPTR(2)=3
       IPTR(3)=2
       IPTR(4)=4
       IPTRI(1)=1
       IPTRI(2)=3
       IPTRI(3)=2
       IPTRI(4)=4
      ELSEIF(PERM.EQ.'1342')THEN 
       CALL ICOPY(NIRREP,NUMP,1,NUMP2,1)
       CALL ICOPY(NIRREP,NUMQ,1,NUMR2,1)
       CALL ICOPY(NIRREP,NUMR,1,NUMS2,1)
       CALL ICOPY(NIRREP,NUMS,1,NUMQ2,1)
       IPTR(1)=1
       IPTR(2)=3
       IPTR(3)=4
       IPTR(4)=2
       IPTRI(1)=1
       IPTRI(2)=4
       IPTRI(3)=2
       IPTRI(4)=3
      ENDIF
      DO 10 IRREPON=1,NIRREP
       IOFFDS1=1
       IOFFDZ1=1
       IOFFDS2=1
       IOFFDZ2=1
       DO 20 IRREPQS=1,NIRREP
        IRREPPR=DIRPRD(IRREPQS,IRREPON)
        DISDSK(IRREPQS,IRREPON)=IOFFDS1
        DSZDSK(IRREPQS,IRREPON)=IOFFDZ1
        DISMEM(IRREPQS,IRREPON)=IOFFDS2
        DSZMEM(IRREPQS,IRREPON)=IOFFDZ2
        IOFFDS1=IOFFDS1+NUMR(IRREPPR)*NUMS(IRREPQS)
        IOFFDZ1=IOFFDZ1+NUMP(IRREPPR)*NUMQ(IRREPQS)
        IOFFDS2=IOFFDS2+NUMR2(IRREPPR)*NUMS2(IRREPQS)
        IOFFDZ2=IOFFDZ2+NUMP2(IRREPPR)*NUMQ2(IRREPQS)
20     CONTINUE
       PQFULL(IRREPON)=IOFFDZ1-1
       DSZPCK(IRREPON)=IRPDPD(IRREPON,ISYTYP(1,LIST))
       IF(IRREPON.EQ.IRREPGET)THEN
        DSZTAR=IOFFDZ2-1
        DISTAR=IOFFDS2-1
       ENDIF
10    CONTINUE
C
C
C LOOP THROUGH ALL IRREPS FOR THIS LIST.
C
      DO 110 IRREPON=1,NIRREP
       LOGREC=0
       FULLSZ=PQFULL(IRREPON)
       PCKSIZ=DSZPCK(IRREPON)
C
C LOOP OVER ALL IRREP PAIRS SUCH THAT DIRPRD(IRREPR,IRREPS)=IRREPDO
C  AND DETERMINE IRREPS OF P2,Q2,R2,S2 (TARGET INDICES - IRPTAR) 
C  FROM THIS INFORMATION.
C
       DO 120 IRREPS=1,NIRREP
        CALL IZERO(IRPTAR,4)
        IRPTAR(IPTR(4))=IRREPS 
        IRREPR=DIRPRD(IRREPS,IRREPON)
        IRPTAR(IPTR(3))=IRREPR
        CALL FILIRP(IRPTAR,IRREPGET)
        IRREPQ=IRPTAR(IPTR(2))
        IRREPP=IRPTAR(IPTR(1))
        IOFFS2=DISMEM(IRPTAR(4),IRREPGET)-1
        IOFFQ2=DSZMEM(IRPTAR(2),IRREPGET)-1
C
C READ ALL LOGICAL RECORDS FOR THIS IRREPR,IRREPS SET
C
        DO 130 ICOUNT4=1,NUMS(IRREPS)
         DO 140 ICOUNT3=1,NUMR(IRREPR)
          LOGREC=LOGREC+1
          CALL GETLST(BUF,LOGREC,1,1,IRREPON,LIST)
          IF(PUTTYP(1:1).EQ.'P')THEN
           CALL SYMEXP2(IRREPON,NUMP,PQFULL(IRREPON),
     &                   DSZPCK(IRREPON),1,BUF,BUF)
           IF(ANTI.EQ.'Y')CALL SSCAL(FULLSZ,HALF,BUF,1)
          ENDIF
C
C NOW EACH RECORD HAS SOME STUFF THAT WE WANT.  LOCATE STARTING 
C  POINT FOR THE BLOCK OF VALUES WHICH WE WANT.
C
          ISTART=DSZDSK(IRREPQ,IRREPON)
          ISTART2=DSZDSK(IRREPP,IRREPON)
C
C NOW LOOP OVER THE Q VALUES, CALCULATE ADDRESSES IN THE
C  TARGET VECTOR AND PUT THE ELEMENTS IN POSITION.
C
          LENGTH2=NUMQ(IRREPQ)  
          DO 150 ICOUNT2=1,NUMQ(IRREPQ)
           LENGTH=NUMP(IRREPP)
C
C CALCULATE THE DISTRIBUTION NUMBER THAT THIS R2,S2 PAIR CORRESPONDS TO
C   WITHIN THE IRREPS2 BLOCK AND ABSOLUTE STARTING ADDRESS FOR THE
C   (XXR2S2) DISTRIBUTION.
C
           DISOFF1=IOFFS2+INDX(ONDISK(IPTRI(3)),ONDISK(IPTRI(4)),
     &                          NUMR2(IRPTAR(3)))
           IADDR1 =INDX(1,DISOFF1,DSZTAR)
C
C NOW CALCULATE RELATIVE OFFSET WITHIN THIS DISTRIBUTION FOR THIS
C  VALUE OF S.
C
           IADDR2 =IOFFQ2+INDX(1,ONDISK(IPTRI(2)),NUMP2(IRPTAR(1)))
C
C CALCULATE TOTAL ABSOLUTE ADDRESS AND PUT ELEMENTS IN POSITION 
C
           OFFSET=IADDR1+IADDR2-1
           IF(AUGTYP.EQ.'S')THEN
            CALL SAXPY(LENGTH,FACTOR,Z(OFFSET),1,BUF(ISTART),1)
            IF(RHF) THEN
             CALL SAXPY(LENGTH,FACTOR2,Z(OFFSET),1,BUF(ISTART2),LENGTH2)
            ENDIF 
           ELSE
            CALL SCOPY(LENGTH,Z(OFFSET),1,BUF(ISTART),1)
           ENDIF
           ISTART=ISTART+LENGTH
           ISTART2=ISTART2+1
150       CONTINUE
          IOFF=1
          IF(PUTTYP(1:1).EQ.'P')THEN
           IOFF=FULLSZ+1
           IF(ANTI.EQ.'N')THEN
            CALL SQSYM(IRREPON,NUMP,PCKSIZ,FULLSZ,1,BUF(IOFF),BUF)
           ELSE
            CALL ASSYM(IRREPON,NUMP,1,1,BUF(IOFF),BUF)
           ENDIF
          ENDIF
          CALL PUTLST(BUF(IOFF),LOGREC,1,ICACHE,IRREPON,LIST)
140      CONTINUE
130     CONTINUE
120    CONTINUE
110   CONTINUE
      RETURN
      END
