      SUBROUTINE TFIJKL(CMO,ICORE,MAXCOR,IUHF)
C
C DRIVER FOR THE TRANSFORMATION OF TWO-ELECTRON INTEGRALS
C OF TYPE IJKL (FULL TRANSFORMATION)
C
CEND
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL RHF,LAST,FIRST
      INTEGER DIRPRD,POP,VRT
      CHARACTER*80 FNAME
C
      DIMENSION CMO(1),ICORE(MAXCOR)
      DIMENSION ISIZE2(8,8),ISIZE3(8,8,8),ISIZE4(8,8,8)
      DIMENSION IOFFAO(8),IOFFMO(8,2),NSTART(8),NEND(8)
      DIMENSION IOFF2(8,8,8),ISIZE3T(8,8),ISIZE3TT(8)
      DIMENSION IOFF4(8,8)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/INFO2/NBASIS,NBAS(8),NMO(8),POP(8,2),VRT(8,2)
      COMMON/SYMINF/NDUMMY,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/AOOFST/INDOCC(8,2)
C
      NNP1O2(I)=(I*(I+1))/2
C
C CREATE LOOK UP VECTOR FOR IRREP OF EACH AO
C
      ISYMAO=1
      IREORD=ISYMAO+NBASIS
      ISTART=IREORD+(1+IUHF)*NBASIS
      IF(MOD(ISTART,2).NE.1) ISTART=ISTART+1
C 
      NCOMP=NOCCO(1)+NVRTO(1)
C
      IND=0
      DO 1 IRREP=1,NIRREP
       DO 1 I=1,NBAS(IRREP)
        IND=IND+1
        ICORE(IND)=IRREP
1     CONTINUE
C
C FILL REORDER VECTOR FOR WRITING INTEGRALS TO DISK
C
      IND=IREORD-1
      DO 5 ISPIN=1,MIN(IUHF+1,2)
       IOFFOCC=0
       IOFFVRT=NOCCO(ISPIN)
       DO 4 IRREP=1,NIRREP
        DO 2 IOCC=1,POP(IRREP,ISPIN)
         IND=IND+1
         ICORE(IND)=IOCC+IOFFOCC
2       CONTINUE
        IOFFOCC=IOFFOCC+POP(IRREP,ISPIN)
        DO 3 IVRT=1,VRT(IRREP,ISPIN)
         IND=IND+1
         ICORE(IND)=IOFFVRT+IVRT
3       CONTINUE
        IOFFVRT=IOFFVRT+VRT(IRREP,ISPIN)
4      CONTINUE
5     CONTINUE
C
C DETERMINE VALUE OF RHF FLAG
C
      RHF=IUHF.EQ.0
C
C LOOP OVER SPIN CASES
C
C      RHF : ONE LOOP
C      UHF : TWO LOOPS
C      ROHF : CURRENTLY TWO LOOPS, LATER THE HANDY-POPLE TRANSFORMATION
C             MIGHT BE IMPLEMENTED
C
      DO 1000 ISPIN=1,MIN(1+IUHF,2)
C
C OPEN INTEGRAL FILE IJIJ (UNIT 10)
C
       LUINT=10
       CALL GFNAME('IJKL     ',FNAME,ILENGTH)
       OPEN(LUINT,FILE=FNAME(1:ILENGTH),STATUS='OLD',
     &      FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       CALL LOCATE(LUINT,'TWOELSUP')
C
C FILL NSTART AND NEND FOR HANDLING OF OUT-OF-CORE
C
       DO 6 IRREP=1,NIRREP
        NSTART(IRREP)=1
        NEND(IRREP)=NMO(IRREP)
6      CONTINUE
C
C CALCULATE SIZE OF SYMMETRIC TRIANGULAR ARRAY FOR EACH IRREP
C
       CALL IZERO(IOFFAO,8)
       CALL IZERO(IOFFMO,16)
       CALL IZERO(ISIZE2,64)
       CALL IZERO(ISIZE3,512)
       CALL IZERO(ISIZE4,512)
       CALL IZERO(ISIZE3T,64)
       CALL IZERO(ISIZE3TT,8)
       CALL IZERO(IOFF2,512)
       CALL IZERO(IOFF4,64)
       NSIZE=0
       MSIZE=0
       IAO=0
       IMO=0
       DO 10 IRREP1=1,NIRREP
C
C OFFSET OF BASIS FUNCTIONS FOR EACH IRREP
C
        IOFFAO(IRREP1)=IAO
        IOFFMO(IRREP1,1)=IMO
        IOFFMO(IRREP1,2)=IMO+NCOMP
C
C SIZE OF TWO-DIMENSION INTEGRAL ARRAY OF IRREP IRREP
C
        IAO=IAO+NBAS(IRREP1)
        IMO=IMO+NMO(IRREP1)
C
C SIZE OF TWO-DIMENSION INTEGRAL ARRAY OF IRREP IRREP
C
        DO 11 IRREP2=1,NIRREP
         IRREP12=DIRPRD(IRREP1,IRREP2)
         IF(IRREP12.NE.1) THEN
          ISIZE2(IRREP1,IRREP2)=NBAS(IRREP1)*NBAS(IRREP2)
C
C SIZE OF TWO- AND THREE-DIMENSIONAL INTEGRAL ARRAY OF IRREP IRREP
C
          DO 12 IRREP3=1,NIRREP
           IF(IRREP1.NE.IRREP3.AND.IRREP2.NE.IRREP3) THEN
            IRREP4=DIRPRD(IRREP3,IRREP12)
            ISIZE3(IRREP3,IRREP1,IRREP2)=NBAS(IRREP3)
     &                 *ISIZE2(IRREP1,IRREP2)
            ISIZE4(IRREP3,IRREP1,IRREP2)=
     &                  ISIZE3(IRREP3,IRREP1,IRREP2)
     &                  *NMO(IRREP4)
           ENDIF
12        CONTINUE
         ENDIF
11      CONTINUE
C
10     CONTINUE
C
C DETERMINE REQUIRED OFFSET
C
       IOFF=0
C
C WE HAVE HERE IRREP4 > IRREP3, IRREP2 > IRREP1 AND
C              IRREP4 > IRREP2
C
C SLOWEST INDEX (OCCUPIED ORBITALS ONLY)
C
       DO 13 IRREP4=4,NIRREP
C
C SECOND SLOWEST INDEX
C
        DO 14 IRREP3=1,NIRREP
         IF(IRREP3.LT.IRREP4) THEN
          IRREP12=DIRPRD(IRREP3,IRREP4)
C
C THIRD SLOWEST INDEX (SHOULD ALWAYS BE GREATER THAN IRREP2)
C
          IOFF22=0
          IOFF4(IRREP4,IRREP3)=IOFF
          ITMP=0
          DO 15 IRREP2=2,NIRREP
           IF(IRREP2.NE.IRREP3.AND.IRREP2.LT.IRREP4) THEN
C
C FASTEST INDEX
C
            IRREP1=DIRPRD(IRREP12,IRREP2)
            IF(IRREP1.LT.IRREP2) THEN
             ITMP=ITMP+NBAS(IRREP1)*NBAS(IRREP2)
             IOFF2(IRREP3,IRREP2,IRREP1)=IOFF22
             IOFF2(IRREP3,IRREP1,IRREP2)=IOFF22
             IOFF22=IOFF22+NBAS(IRREP1)*NBAS(IRREP2)
             IOFF=IOFF+ISIZE4(IRREP3,IRREP1,IRREP2)
            ENDIF
           ENDIF
15        CONTINUE
          MSIZE=MAX(MSIZE,NBAS(IRREP3)*ITMP,
     &              2*NBAS(IRREP1)*NBAS(IRREP2),
     &              2*NBAS(IRREP3)*NBAS(IRREP4))
         ENDIF
14      CONTINUE
13     CONTINUE
C
       DO 16 IRREP1=4,NIRREP
        DO 17 IRREP2=1,NIRREP
         IF(IRREP2.GE.IRREP1) GO TO 17
         IRREP12=DIRPRD(IRREP1,IRREP2)
         DO 18 IRREP3=2,NIRREP
          IF(IRREP3.EQ.IRREP2.OR.IRREP3.GE.IRREP1) GO TO 18
           IRREP4=DIRPRD(IRREP3,IRREP12)
           IF(IRREP4.LT.IRREP3) THEN
            ISIZE3T(IRREP1,IRREP2)=ISIZE3T(IRREP1,IRREP2)+
     &                      ISIZE3(IRREP2,IRREP3,IRREP4)
            ISIZE3TT(IRREP1)=ISIZE3TT(IRREP1)+
     &                      ISIZE3(IRREP2,IRREP3,IRREP4)
           ENDIF
19        CONTINUE
18       CONTINUE
17      CONTINUE
16     CONTINUE
C
C DETERMINE TOTAL SIZE OF FOUR-DIMENSIONAL INTEGRAL ARRAY
C
       NSIZE=IOFF
C
C INTEGRAL BUFFER LENGTH
C
       ILNBUF=600
C
C WE IMPLEMENT AT THE MOMENT ONLY FULL IN CORE MBPT(2) TRANSFORNATIONS
C 
       I000=ISTART
C
C I000 .... FOUR DIMENSIONAL INTEGRAL ARRAY OF SIZE NSIZE
C
       I010=I000+NSIZE*IINTFP
       IEND2=I010+IINTFP*MSIZE
C
C I010 .... SCRATCH ARRAY OF MAXIMUM SIZE NBASIS
C
       I020=I010+IINTFP*NBASIS
C
C I020 .... INTEGRAL VALUES OF ONE RECORD READ IN FROM IIJJ
C
       I030=I020+ILNBUF*IINTFP
C 
C I030 .... INTEGRAL INDICES OF ONE RECORD READ IN FROM IIJJ
C
       IEND1=I030+ILNBUF
       IEND=MAX(IEND1,IEND2)
C
       IOFFSET=0 
C
C CHECK WHETHER FULL IN-CORE ALGORITHM IS POSSIBLE
C
       IF(IEND.GE.MAXCOR) THEN
C
        I010=ISTART
        I020=I010+IINTFP*NBASIS
        I030=I020+IINTFP*ILNBUF
        IEND1=I030+ILNBUF
        IEND2=I010+IINTFP*MSIZE
        IEND=MAX(IEND1,IEND2)
        IF(MOD(IEND,2).NE.1) IEND=IEND+1
        I000=IEND
C
C MXCOR IS MAXIMUM OF AVAILABLE CORE MEMORY
C
        MXCOR=MAXCOR-IEND
C 
C SET FIRST ORBITAL TO 1
C
        JSTART=NMO(1)+NMO(2)+NMO(3)+1
C
C SET JEND TO TOTAL NUMBER OF CORRELATED ORBITALS
C
        JEND=0
        DO 915 IRREP=1,NIRREP
         JEND=JEND+NMO(IRREP)
915     CONTINUE
C
       ELSE
C
        NSIZE1=NSIZE
        LAST=.TRUE.
        GO TO 950
C
       ENDIF
C
900    CONTINUE
C
C RESET MXCOR1
C
       MXCOR1=MXCOR
C
       NSIZE1=0
       CALL IZERO(NEND,8)
       DO 901 IRREP=1,NIRREP
        NSTART(IRREP)=1
901    CONTINUE
C
C DETERMINE NSTART
C
       JST=JSTART
C
       DO 910 IRREP=1,NIRREP
        JRREP=IRREP
        IF(JST.LE.NMO(IRREP)) GO TO 920
        JST=JST-NMO(IRREP)
910    CONTINUE
920    CONTINUE
C
       NSTART(JRREP)=JST
C
       FIRST=.TRUE.
C
930    CONTINUE
C
       NLEFT=NMO(JRREP)+1-NSTART(JRREP)   
C
       IF(ISIZE3TT(JRREP).NE.0) THEN
        NORB=MXCOR1/(ISIZE3TT(JRREP)*IINTFP)
        IF(NORB.EQ.0.AND.FIRST) THEN
         write(*,*) ' @TFIJKL-F, Sorry, there is not enough memory !'
         CALL ERREX
        ENDIF
       ELSE
        NORB=NMO(JRREP)+1
       ENDIF
       FIRST=.FALSE.
       NEND(JRREP)=NSTART(JRREP)+MIN(NLEFT,NORB)-1
       MXCOR1=MXCOR1-MIN(NLEFT,NORB)*ISIZE3TT(JRREP)*IINTFP
       JSTART=JSTART+MIN(NLEFT,NORB)
       NSIZE1=NSIZE1+MIN(NLEFT,NORB)*ISIZE3TT(JRREP)
       LAST=.FALSE.
       IF(JSTART.GT.JEND) LAST=.TRUE.
C
       IF((.NOT.LAST).AND.(NORB.GT.NLEFT)) THEN
        JRREP=JRREP+1
        NSTART(JRREP)=1
        GO TO 930
       ENDIF
C
C DETERMINE REQUIRED OFFSET
C
       IOFF=0
C
C WE HAVE HERE IRREP4 > IRREP3, IRREP2 > IRREP1 AND
C              IRREP4 > IRREP2
C
C SLOWEST INDEX (OCCUPIED ORBITALS ONLY)
C
       DO 26 IRREP4=4,NIRREP
C
C SECOND SLOWEST INDEX
C
        DO 27 IRREP3=1,NIRREP
         IF(IRREP3.LT.IRREP4) THEN
          IRREP12=DIRPRD(IRREP3,IRREP4)
C
C THIRD SLOWEST INDEX (SHOULD ALWAYS BE GREATER THAN IRREP2)
C
          IOFF22=0
          IOFF4(IRREP4,IRREP3)=IOFF
          DO 28 IRREP2=2,NIRREP
           IF(IRREP2.NE.IRREP3.AND.IRREP2.LT.IRREP4) THEN
C
C FASTEST INDEX
C
            IRREP1=DIRPRD(IRREP12,IRREP2)
            IF(IRREP1.LT.IRREP2) THEN
             IOFF=IOFF+ISIZE3(IRREP3,IRREP1,IRREP2)*
     &            (NEND(IRREP4)-NSTART(IRREP4)+1)
            ENDIF
           ENDIF
28        CONTINUE
         ENDIF
27      CONTINUE
26     CONTINUE
C
       NSIZE1=IOFF
C
950    CONTINUE
C
       CALL LOAD4F(CMO,ICORE(I000),ICORE(I010),ICORE(I020),
     &             ICORE(I030),ICORE(ISYMAO),NBAS,
     &             NSTART,NEND,IOFFAO,IOFF4,ISIZE3T,
     &             IOFF2,NSIZE1,ILNBUF,ISPIN,LUINT,
     &             LAST,IOFFSET)
C
       CALL INTRN4F(CMO,ICORE(I000),ICORE(I010),RHF,
     &              NBASIS,NBAS,NMO,NSTART,NEND,ISIZE3,
     &              IOFFMO,NSIZE1,ISPIN,ICORE(IREORD))
       IOFFSET=IOFFSET+NSIZE1
       IF(.NOT.LAST) GO TO 900
1000   CONTINUE        
       RETURN 
       END
