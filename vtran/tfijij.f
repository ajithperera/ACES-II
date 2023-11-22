      SUBROUTINE TFIJIJ(CMO,ICORE,MAXCOR,IUHF)
C
C DRIVES THE TRANSFORMATIONOF IJIJ INTEGRALS IN CASE
C OF A FULL AO->MO INTEGRAL TRANSFORMATION.
C
CEND
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL RHF,LAST,FIRST
      INTEGER DIRPRD,POP,VRT
      CHARACTER*80 FNAME
C
      DIMENSION CMO(1),ICORE(MAXCOR)
      DIMENSION ISIZE2(8,8),ISIZE3(8,8),ISIZE4(8,8)
      DIMENSION IOFFAO(8),IOFFMO(8,2),ISIZE3T(8)
      DIMENSION IOFF2(8,8),IOFF3(8,8)
      DIMENSION IOFF4(8,8)
      DIMENSION NSTART(8),NEND(8)
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
       CALL GFNAME('IJIJ    ',FNAME,ILENGTH)
       OPEN(LUINT,FILE=FNAME(1:ILENGTH),STATUS='OLD',
     &      FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       CALL LOCATE(LUINT,'TWOELSUP')
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
       CALL IZERO(ISIZE3,64)
       CALL IZERO(ISIZE3T,8)
       CALL IZERO(ISIZE4,64)
       CALL IZERO(IOFF2,64)
       CALL IZERO(IOFF3,64)
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
        IAO=IAO+NBAS(IRREP1)
        IMO=IMO+NMO(IRREP1)
C
C SIZE OF TWO-DIMENSION INTEGRAL ARRAY OF IRREP IRREP
C
        DO 11 IRREP2=1,IRREP1-1
         ISIZE2(IRREP2,IRREP1)=NBAS(IRREP1)*NBAS(IRREP2)
         ISIZE3(IRREP2,IRREP1)=NBAS(IRREP1)*ISIZE2(IRREP2,IRREP1)
         ISIZE4(IRREP2,IRREP1)=ISIZE3(IRREP2,IRREP1)
     &                         *NMO(IRREP2)
11      CONTINUE
C
10     CONTINUE
C
       DO 13 IRREP1=1,NIRREP
        DO 13 IRREP2=1,NIRREP
         ISIZE3T(IRREP1)=ISIZE3T(IRREP1)+ISIZE3(IRREP1,IRREP2)
13     CONTINUE
C
C DETERMINE REQUIRED OFFSET
C
       IOFF=0
       DO 15 IRREP1=1,NIRREP
        DO 16 IRREP2=IRREP1+1,NIRREP
C
         IOFF4(IRREP1,IRREP2)=IOFF
         IOFF=IOFF+ISIZE4(IRREP1,IRREP2)
         MSIZE=MAX(MSIZE,ISIZE3(IRREP1,IRREP2),
     &             2*NBAS(IRREP1)*NBAS(IRREP2))
C
16      CONTINUE
15     CONTINUE
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
C MULTIPLE PASSES ARE REQUIRED (OUT-OF-CORE ALGORITHM)
C
C ALLOCATE CORE MEMORY FOR OUT-OF-CORE
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
C SET FIRST ORBITAL (JSTART) TO 1
C
       JSTART=1
C
C SET JEND TO TOTAL NUMBER OF ALL ORBITALS
C
       JEND=0
       DO 915 IRREP=1,NIRREP-1
        JEND=JEND+NMO(IRREP) 
915    CONTINUE
C
       ELSE
C
        NSIZE1=NSIZE
        LAST=.TRUE.
        GO TO 950
C
       ENDIF
C
900   CONTINUE
C
C RESET MXCOR1
C
      MXCOR1=MXCOR
C
      NSIZE1=0
      CALL IZERO(NEND,8)
      DO 901 IRREP=1,NIRREP
       NSTART(IRREP)=1
901   CONTINUE
C
      JST=JSTART 
C
      DO 910 IRREP=1,NIRREP
       JRREP=IRREP
       IF(JST.LE.NMO(IRREP)) GO TO 920 
       JST=JST-NMO(IRREP)
910   CONTINUE
920   CONTINUE
C
      NSTART(JRREP)=JST
C
      FIRST=.TRUE.
C
930   CONTINUE
C
      NLEFT=NMO(JRREP)+1-NSTART(JRREP)
C
      IF(ISIZE3T(JRREP).NE.0) THEN
       NORB=MXCOR1/(ISIZE3T(JRREP)*IINTFP)
       IF(NORB.EQ.0.AND.FIRST) THEN
        write(*,*) '@TFIJIJ-F, Sorry, there is not enough memory !'
        CALL ERREX
       ENDIF
      ELSE
       NORB=NMO(JRREP)+1
      ENDIF
      FIRST=.FALSE.
      NEND(JRREP)=NSTART(JRREP)+MIN(NLEFT,NORB)-1
      MXCOR1=MXCOR1-MIN(NLEFT,NORB)*ISIZE3T(JRREP)*IINTFP
      JSTART=JSTART+MIN(NLEFT,NORB)
      NSIZE1=NSIZE1+MIN(NLEFT,NORB)*ISIZE3T(JRREP)
      LAST=.FALSE.
      IF(JSTART.GT.JEND) LAST=.TRUE.
C
C CHECK IF MORE IRREPS HAVE TO BE CONSIDERED IN THIS PASS
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
       DO 17 IRREP1=1,NIRREP
        DO 18 IRREP2=IRREP1+1,NIRREP
C
         IOFF4(IRREP1,IRREP2)=IOFF
         IOFF=IOFF+ISIZE3(IRREP1,IRREP2)*
     &            (NEND(IRREP1)-NSTART(IRREP1)+1)
C
18      CONTINUE
17     CONTINUE
C
C DETERMINE TOTAL SIZE OF FOUR-DIMENSIONAL INTEGRAL ARRAY
C
       NSIZE1=IOFF
C
950   CONTINUE
C
       CALL LOAD3F(CMO,ICORE(I000),ICORE(I010),ICORE(I020),
     &             ICORE(I030),ICORE(ISYMAO),NBAS,
     &             NMO,NSTART,NEND,IOFFAO,IOFF4,ISIZE3,
     &             NSIZE1,ILNBUF,ISPIN,LUINT,LAST,
     &             IOFFSET)
C
       CALL INTRN3F(CMO,ICORE(I000),ICORE(I010),RHF,
     &              NBASIS,NBAS,NMO,NSTART,NEND,
     &              ISIZE3,IOFFMO,NSIZE1,ISPIN,
     &              ICORE(IREORD))
C
       IOFFSET=IOFFSET+NSIZE1
C
       IF(.NOT.LAST) GO TO 900
C
1000   CONTINUE        
       RETURN 
       END
