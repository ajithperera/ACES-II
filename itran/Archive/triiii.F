      SUBROUTINE TRIIII(CMO,ICORE,MAXCOR,IUHF)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL RHF,DOALL,NOABCD,ABIJ,LAST,FIRST,SINGLE_PASS
      INTEGER DIRPRD,POP,VRT,HalfGB
      CHARACTER*80 FNAME
C
#include <maxbasfn.par>
C  
C
      DIMENSION CMO(maxbasfn*maxbasfn),ICORE(MAXCOR)
      DIMENSION ISIZE2(8),ISIZE3(8),ISIZE4(8),IOFFI(8),
     &          IOFFAO(8),IOFFMO(8,2)
      DIMENSION NFIRST(8),NSTART(8),NEND(8),
     &          CMO_SORTED(maxbasfn*maxbasfn)
C   
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/INFO2/NBASIS,NBAS(8),NMO(8),POP(8,2),VRT(8,2)
      COMMON/SYMINF/NDUMMY,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/AOOFST/INDOCC(8,2)
      COMMON/DOINTS/DOALL,ABIJ,NOABCD
C
      NNP1O2(I)=(I*(I+1))/2
C
C CREATE LOOK UP VECTOR FOR IRREP OF EACH AO
C
      MAXCOR = 4400
      ISYMAO=1
      IREORD=ISYMAO+NBASIS
      ISTART=IREORD+(1+IUHF)*NBASIS
      IF(MOD(ISTART,2).NE.1) ISTART=ISTART+1
      SINGLE_PASS = .TRUE. 
            IPASS = 0
C
      NCOMP=NOCCO(1)+NVRTO(1)
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
         ICORE(IND)=IVRT+IOFFVRT
3       CONTINUE
        IOFFVRT=IOFFVRT+VRT(IRREP,ISPIN)
4      CONTINUE
5     CONTINUE
C
#ifdef _NO_NOSKIP
      Write(6,*) "Unsorted MO vector"
      call output(CMO, 1, Nbasis, 1, Nbasis, Nbasis, 1)
      IMO_COUNT = 1
      Do I = 1, Nbasis
         IMO_INDEX = I
         Do J = 1, Nbasis
            CMO_SORTED(IMO_COUNT) = CMO(IMO_INDEX)
                        IMO_INDEX = IMO_INDEX + Nbasis
                        IMO_COUNT = IMO_COUNT + 1
         End Do
C         IMO_COUNT = IMO_COUNT + NBasis
      End Do
      Write(6,*) "Sorted MO vector"
      call output(CMO_SORTED, 1, Nbasis, 1, Nbasis, Nbasis, 1)
#endif
C
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
C OPEN INTEGRAL FILE IIII (UNIT 10)
C
       LUINT=10
       CALL GFNAME('IIII    ',FNAME,ILENGTH)
       OPEN(LUINT,FILE=FNAME(1:ILENGTH),STATUS='OLD',
     &      FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       CALL LOCATE(LUINT,'TWOELSUP')
C
C FILL NFIRST
C
       IF(NOABCD) THEN
        CALL ICOPY(NIRREP,POP(1,ISPIN),1,NFIRST,1)
       ELSE
        CALL ICOPY(NIRREP,NMO,         1,NFIRST,1)
       ENDIF
C
C   NSTART AND NEND PARAMETER FOR FULL IN-CORE TRANSFORMATION
C
       DO 6 IRREP=1,NIRREP
        NSTART(IRREP)=1
        NEND(IRREP)=NFIRST(IRREP)
6      CONTINUE
C
C CALCULATE SIZE OF SYMMETRIC TRIANGULAR ARRAY FOR EACH IRREP
C
       CALL IZERO(IOFFI,8)
       CALL IZERO(ISIZE2,8)
       CALL IZERO(ISIZE3,8)
       CALL IZERO(ISIZE4,8)
       CALL IZERO(IOFFAO,8)
       CALL IZERO(IOFFMO,16)
            NSIZE = 0
            MSIZE = 0
       MSIZE_AOMO = 0
              IAO = 0
              IMO = 0
       DO 10 IRREP=1,NIRREP
C
C OFFSET FOR BASIS FUNCTIONS
C
        IOFFAO(IRREP)=IAO
        IOFFMO(IRREP,1)=IMO
        IOFFMO(IRREP,2)=IMO+NCOMP
C
C OFFSET FOR INTEGRAL ARRAY
C
        IOFFI(IRREP)=NSIZE
C
C SIZE OF A TRINAGULAR INTEGRAL ARRAY WITHIN IRREP IRREP
C
        ISIZE2(IRREP)=NNP1O2(NBAS(IRREP))
C
C SIZE OF A THREE_DIMENSIONAL INTEGRAL ARRAY WITHIN IRREP IRREP
C
        ISIZE3(IRREP)=ISIZE2(IRREP)*NBAS(IRREP)
C
C SIZE OF FOUR-DIMENSIONAL INTEGRAL ARRAY WITHIN IRREP IRREP
C
        ISIZE4(IRREP)=ISIZE3(IRREP)*NFIRST(IRREP)
C
#ifdef _READ_ONLY
C If Load1 only do reading the ints then ISIZE4 must only be
C ISIZE2(IRREP)**2.
C
        ISIZE4(IRREP)=ISIZE2(IRREP)*ISIZE2(IRREP)
C
#endif         
C
C TOTAL SIZE OF INTEGRAL ARRAY
C
        NSIZE=NSIZE+ISIZE4(IRREP)
C
C MAXIMUM SIZE OF THREE DIMSNIONAL ARRAY (CHECK ALL IRREPS)

        MSIZE=MAX(MSIZE,ISIZE3(IRREP),2*NBAS(IRREP)*NBAS(IRREP))
C
#ifdef _READ_ONLY
C Since this location is used to keep a NBAS(IRREP)*NBAS(IRREP) for
C a given NNP1O2(NBAS(IRREP)) and then one index transformed integral
C which are of ISIZE3(IRREP)*NFIRST(IRREP) for all the irreps, the 
C MSIZE need to be increased. Also, this location is also used to store 
C a copy of the AO integrals while building  the full integral file
C (note that what is read in load1.F is mu>nu  for a given lam>sig
C and this need to be expanded and that requires a copy of the 
C AO integralr read in load1.F (see trnsfrm_1sindex.F). So, the 
CSSS
CSSS        MAX_BFNS_IRREP = MAX(MAX_BFNS_IRREP, NBAS(IRREP))
        MSIZE = MSIZE + ISIZE4(IRREP) + NBAS(IRREP)*NBAS(IRREP) 
        Print*, "The MSIZE:", MSIZE
#endif
C
C COUNTER FOR BASIS FUNCTIONS
C
        IAO=IAO+NBAS(IRREP)
        IMO=IMO+NMO(IRREP)
C
10     CONTINUE

#ifdef _READ_ONLY
C Since this location is used to keep a NBAS(IRREP)*NBAS(IRREP) for
C a given NNP1O2(NBAS(IRREP)) and then one index transformed integral
C which are of ISIZE3(IRREP)*NFIRST(IRREP) for all the irreps, the
C MSIZE need to be increased. Also, this location is also used to store
C a copy of the AO integrals while building  the full integral file
C (note that what is read in load1.F is mu>nu  for a given lam>sig
C and this need to be expanded and that requires a copy of the
C AO integralr read in load1.F (see trnsfrm_1sindex.F). So, the
     
        MSIZE = NSIZE + NBASIS*NBASIS 
        Print*, "The MSIZE:", MSIZE
#endif
C     
C
C INTEGRAL BUFFER LENGTH
C
       ILNBUF=600
C
C TRY FIRST FULL IN CORE (PARTIAL OR FULL) TRANSFORNATIONS
C
       I000=ISTART
C
C I000 .... INTEGRAL ARRAY OF SIZE NSIZE
C
       I010=I000+NSIZE*IINTFP
C
C I010 .... SCRATCH ARRAY OF MAXIMUM SIZE NBASIS
C
       I020=I010+IINTFP*NBASIS
C
C I020 .... INTEGRAL VALUES FOR EACH RECORD READ FROM IIII
C
       I030=I020+ILNBUF*IINTFP
C
C I030 .... INTEGRAL INDICES FOR EACH RECORD READ FROM IIII
C
C
       I040=I030+ILNBUF
C         
C I040 ....Temporarily keep the one index transformed MOs.
C
       IEND1=I030+ILNBUF
       IEND2=I010+IINTFP*MSIZE
       IEND=MAX(IEND1,IEND2)
C
#ifdef _READ_ONLY
       I050=I040+NBASIS*NBASIS*IINTFP
       IEND=I050+MSIZE*IINTFP
#endif
C
C CHECK IF FULL IN-CORE TRANSFORMATION IS POSSIBLE
C
#ifdef _DEBUG_LVL0
       Write(6,*) "Choose Betweem in core, out of core, 
     & IEND, MAXCOR", IEND, MAXCOR
       Write(6,*)
#endif
C 
       IF (IEND.GE.MAXCOR) THEN
          SINGLE_PASS = .FALSE. 
#ifdef _DEBUG_LVL0
       Write(6,*), "Multipass Algorithm is being used"
#endif
C
C MULTIPLE PASSES ARE REQUIRED (OUT-OF-CORE ALGORITHM)
C
C ALLOCATE CORE MEMORY FOR OUT-OF CORE
C
C  WE NEED DEFINITELY MSIZE*IINTFP FOR THE SECOND STEP
C  AND NBASIS*IINTFP, ILNBUF*(1+IINTFP) FOR THE
C  FIRST STEP, SUBTRACT FROM TOTAL MEMORY AND THEN
C  DETERMINE MAXIMUM NUMBER OF ORBITALS POSSIBLE
C  TO USE IN FIRST STEP
C
        I010=ISTART
        I020=I010+IINTFP*NBASIS
        I030=I020+IINTFP*ILNBUF
        IEND1=I030+ILNBUF
        IEND2=I010+IINTFP*MSIZE
        IEND=MAX(IEND1,IEND2)
C
#ifdef _READ_ONLY
       I040=I030
       I050=I040+NBASIS*NBASIS*IINTFP
       IEND=I050+MSIZE*IINTFP
#endif 
C
        IF(MOD(IEND,2).NE.1) IEND=IEND+1
        I000=IEND
C
C MXCOR IS MAXIMUM OF AVAILABLE CORE MEMORY
C
        MXCOR=MAXCOR-IEND
        If (MXCOR .LE. 0) CALL INSMEM("@-TRIIII", NSIZE, MAXCOR)
#ifdef _DEBUG_LVLM1
        Print*, "This is what is left for Bufs", MXCOR
#endif
C
C SET FIRST ORBITAL TO 1
C
        JSTART=1
C
C SET JEND TO TOTAL NUMBER OF ALL (OCCUPIED) ORBITALS
C
        JEND=0
        DO 915 IRREP=1,NIRREP
         JEND=JEND+NFIRST(IRREP)
915     CONTINUE
C
       ELSE
C
C FULL INCORE ALGORITHM IS POSSIBLE, SET LAST TO TRUE
C
        Write(6,*), "Single pass Algorithm is being used"
        NSIZE1=NSIZE
        LAST=.TRUE.
C
C In this case all the I<J; K<L handled in a one passs. Set the
C number of KL pairs handlled to NAO(NAO+1)/2
C
        KLPAIRS = NBASIS*(NBASIS+1)/2 
        GO TO 950

       ENDIF
900    CONTINUE
C
C RESET MXCOR1
C
       MXCOR   = 1100
       MXCOR1  = MXCOR 
C
C DETERMINE NSTART AND NEND FOR OUT-OF-CORE TRANSFORMATION. NOTE
C THE OUT-OF-CORE TRANSFORMATION IS DONE AS A MULTIPLE PASS
C ALGORITHM WITH SEVERAL READS OF THE IIII FILE
C
C SET NSIZE1 TO ZERO
C
       NSIZE1=0
C
C  ZERO NEND, NSTART SHOULD STILL BE 1
C
       CALL IZERO(NEND,8)
       CALL IZERO(IOFFI,8)
       DO 901 IRREP=1,NIRREP
        NSTART(IRREP)=1
901    CONTINUE
C
C  DETERMINE STARTING ADDRESS USING JSTART
C
       JST=JSTART
       DO 910 IRREP=1,NIRREP
        JRREP=IRREP
        IF(JST.LE.NFIRST(IRREP)) GO TO 920
        JST=JST-NFIRST(IRREP)
910    CONTINUE
C
920    CONTINUE
C
       NSTART(JRREP)=JST
C
       FIRST=.TRUE.
C
C WITHIN THIS IRREP THERE ARE NEND+1-NSTART ORBITALS LEFT
C DETERMINE HOW MANY FIT INTO CORE
C
930    CONTINUE
       
       IPASS = IPASS + 1
       NLEFT=NFIRST(JRREP)+1-NSTART(JRREP)
C
#ifdef _READ_ONLY
C Instead of the occupied (or totals) we need to go over all.
C
       NLEFT=NBAS(JRREP)+1-NSTART(JRREP)
#endif 
C       
C
C DETERMINE HOW MUCH CORE MEMORY IS REQUIRED
C
C TAKE CARE ABOUT CASES WHERE ISIZE3 IS ZERO
C
#ifdef _DEBUG_LVL0
       Print*, "The MXCOR1 and ISIZE3(JRREP)",MXCOR1,ISIZE3(JRREP)
#endif
       IF(ISIZE3(JRREP).NE.0) THEN
        NORB=MXCOR1/(ISIZE3(JRREP)*IINTFP)
        IF(NORB.EQ.0.AND.FIRST) THEN
         write(*,*) ' @TRIIII-F, Sorry, there is not enough memory !'
         CALL ERREX
        ENDIF
       ELSE
        NORB=NMO(JRREP)+1
       ENDIF
       FIRST=.FALSE.
       NEND(JRREP)=NSTART(JRREP)+MIN(NLEFT,NORB)-1
       MXCOR1=MXCOR1-MIN(NLEFT,NORB)*IINTFP*ISIZE3(JRREP)
       JSTART=JSTART+MIN(NLEFT,NORB)
       IOFFI(JRREP)=NSIZE1
       NSIZE1=NSIZE1+MIN(NLEFT,NORB)*ISIZE3(JRREP)
       LAST=.FALSE.
       IF(JSTART.GT.JEND) LAST=.TRUE.
C
C
C In this case all the I<J; K<L handled in a one passs. Set the
C number of KL pairs handlled to NAO(NAO+1)/2
C
        KLPAIRS = NBAS(JRREP)*(NORB+1)/2 

C CHECK IF MORE IRREPS HAVE TO BE CONSIDERED IN THIS PATH
C
       IF((.NOT.LAST).AND.(NORB.GT.NLEFT)) THEN
        JRREP=JRREP+1
        NSTART(JRREP)=1
        GO TO 930
       ENDIF
C
950    CONTINUE
#ifdef _DEBUG_LVL0
       Write(6,*)
       Write(6,*) "NSIZE1, NORB, NLEFT", NSIZE1, NORB, NLEFT
       Write(6,*) "The KLPAIRS:", KLPAIRS
       Write(6,*)
#endif
#ifdef _DEBUG_LVL0
       Print*, "Reading AO ints"
       Call Start_time
#endif 
       CALL LOAD1(CMO,ICORE(I000),ICORE(I010),ICORE(I020),
     &            ICORE(I030),ICORE(ISYMAO),NBAS,
     &            NFIRST,NSTART,NEND,ISIZE3,IOFFAO,
     &            IOFFI,NSIZE1,ILNBUF,ISPIN,LUINT,LAST)
#ifdef _DEBUG_LVL0
       Print*, "Transforming the first index"
C       Print*, "The icore addresses, MSIZE1, I000, I010, 1040"
C       Print*, MSIZE1, I000, I010, I040 
       Call End_time
       Call Start_time
#endif
C
C This directive is only for debugging and it need to be out
C after all is tested. At the moment we can do both ways. 
C
#ifdef _READ_ONLY
       IF (SINGLE_PASS) THEN
          CALL TRNSFORM_1INDX_1PASS(CMO,ICORE(I000),ICORE(I050),
     &                              ICORE(I040),KLPAIRS,NFIRST,
     &                              NSTART,NEND,ISPIN,NBAS,NMO,
     &                              NBASIS)
       ELSE
          IF (IPASS .EQ. 1) 
     &        CALL LOAD1(CMO,ICORE(I000),ICORE(I010),ICORE(I020),
     &                   ICORE(I030),ICORE(ISYMAO),NBAS,
     &                   NFIRST,NSTART,NEND,ISIZE3,IOFFAO,
     &                   IOFFI,NSIZE1,ILNBUF,ISPIN,LUINT,LAST)

          CALL TRNSFORM_1INDX_MPASS(CMO,ICORE(I000),ICORE(I050),
     &                              ICORE(I040),KLPAIRS,NFIRST,                      &                              NSTART,NEND,ISPIN,NBAS,NMO,
     &                              NORB,NBASIS)
       END IF

#endif
C
#ifdef _DEBUG_LVL0
       Call End_time
       Print*, "Other 3 indices transformation"
       Call Start_time
#endif
C
C The starting address of I010 must be the largest of 2*NBAS(IRREP)*NBAS(IRREP)
C NBAS(IRREP)*(NBAS(IRREP)+1)/2*NBAS(IRREP) (which is ISIZE3(IRREP).
C
       CALL INTRN1(CMO,ICORE(I000),ICORE(I010),
     &             RHF,NBASIS,NBAS,NMO,POP(1,ISPIN),
     &             POP(1,3-ISPIN),NSTART,NEND,ISIZE3,
     &             IOFFMO,NSIZE,ISPIN,LAST,ICORE(IREORD))
#ifdef _DEBUG_LVL0
       Call End_time
#endif
C
       IF(.NOT.LAST) GO TO 900
1000   CONTINUE
       RETURN
       END
