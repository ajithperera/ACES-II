      SUBROUTINE TRIIII(CMO,ICORE,MAXCOR,IUHF)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL RHF,DOALL,NOABCD,ABIJ,LAST,FIRST,ONE_PASS_FONCE,
     &        TWO_PASS_FONCE, MUL_PASS_PONCE
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
C
            NSIZE     = 0
            MSIZE     = 0
              IAO     = 0
              IMO     = 0
#ifdef _2HOT
       MAX_AOS_IRREP = 0
       MAX_MOS_IRREP  = 0
       ISIZE_AUX_AOMO = 0
       MAX_AOFMO_IRREP= 0
       MAX_AOPMO_IRREP= 0
#endif
C
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
        MSIZE = NSIZE + NBAS(IRREP)*NBAS(IRREP)
        Print*, "The NSIZE and MSIZE:", NSIZE, MSIZE
        Write(6,*)
#endif
C
#ifdef _2HOT
C a copy of the AO integrals while building  the full integral file
C (note that what is read in load1.F is mu>nu  for a given lam>sig
C and this need to be expanded and that requires a copy of the 
C AO integralr read in load1.F (see trnsfrm_1sindex.F). So, the 

        MAX_AOS_IRREP = MAX(MAX_BFNS_IRREP, NBAS(IRREP))
        MAX_MOS_IRREP = MAX(MAX_MOS_IRREP,  NFIRST(IRREP))
        ISIZE_AOFMO   = ISIZE_AOFMO + ISIZE3(IRREP)*NFIRST(IRREP)
        ISIZE_AOPMO   = ISIZE_AOPMO + ISIZE3(IRREP) 
C
        ITMP_MAX1       = MAX(ISIZE3(IRREP)*NFIRST(IRREP), MSIZE)
        ITMP_MAX2       = MAX(ISIZE3(IRREP), MSIZE)
        MAX_AOFMO_IRREP = MAX(MAX_AOFMO_IRREP, ITMP_MAX1)
        MAX_AOPMO_IRREP = MAX(MAX_AOPMO_IRREP, ITMP_MAX2)
        ISIZE_AUX_AOMO  = ISIZE_AUX_AOMO + MAX_AOFMO_IRREP 
C
        Print*, "Maximums for memory allocations"
        write(*,"(5(1x,i5))"), MAX_AOS_IRREP,MAX_MOS_IRREP,
     &        MAX_AOFMO_IRREP, MAX_AOPMO_IRREP, ISIZE_AUX_AOMO
#endif  
C
C COUNTER FOR BASIS FUNCTIONS
C
        IAO=IAO+NBAS(IRREP)
        IMO=IMO+NMO(IRREP)
C
10     CONTINUE

#ifdef _2HOT
C Since this location is used to keep a NBAS(IRREP)*NBAS(IRREP) for
C a given NNP1O2(NBAS(IRREP)) and then one index transformed integral
C which are of ISIZE3(IRREP)*NFIRST(IRREP) for all the irreps, the
C MSIZE need to be increased. Also, this location is also used to store
C a copy of the AO integrals while building  the full integral file
C (note that what is read in load1.F is mu>nu  for a given lam>sig
C and this need to be expanded and that requires a copy of the
C AO integralr read in load1.F (see trnsfrm_1sindex.F). So, the
C
C INTEGRAL BUFFER LENGTH
C
      ILNBUF=600
       
      ISIZE_OF_AOINTS = NSIZE
      ISIZE_4_1FINDEX = MAX_BFNS_IRREP*MAX_BFNS_IRREP + 
     &                  ISIZE_AOFMO
      ISIZE_4_1PINDEX = MAX_BFNS_IRREP*MAX_BFNS_IRREP + 
     &                  ISIZE_AOPMO
C
      IMEM_OF_AUX       = ILNBUF*IINTFP + ILNBUF + MAX_BFNS_IRREP*
     &                    MAX_BFNS_IRREP*IINTFP  + MAX_BFNS_IRREP*
     &                    MAX_MOS_IRREP*IINTFP
C
      IMEM_OF_AOS_1PASS = (NSIZE + MAX(NSIZE,ISIZE_AOFMO))*IINTFP +
     &                     MAX_AOS_IRREP*MAX_AOS_IRREP*IINTFP     +
     &                     IMEM_OF_AUX

      IMEM_OF_AOS_2PASS = (NSIZE + ISIZE_AOFMO)*IINTFP + IMEM_OF_AUX

      IMEM_OF_AOS_MPASS = (NSIZE + ISIZE_AOPMO)*IINTFP + IMEM_OF_AUX

      IF ((IMEM_OF_AOS_1PASS + ISTART) .LT. MAXCOR) THEN
          ONE_PASS_FONCE = .TRUE.
      ELSE IF ((IMEM_OF_AOS_2PASS + ISTART) .LT. MAXCOR) THEN
          TWO_PASS_FONCE = .TRUE.
      ELSE IF ((IMEM_OF_AOS_MPASS + ISTART) .LT. MAXCOR) THEN
          MUL_PASS_PONCE = .TRUE.
          MEM_LEFT_4TRNF = MAXCOR - IMEM_OF_AOS_MPASS
      ELSE
          CALL INSMEM("@-TRIIII", (IMEM_OF_AOS_MPASS + ISTART), 
     &                 MAXCOR)
      ENDIF
#ifdef _DEBUG_LVL0
      Write(6,*)
      Print*, "Memory requirements for one, two or mul pass"
      write(*,"(3(1x,i5))"),IMEM_OF_AOS_1PASS, IMEM_OF_AOS_2PASS,
     &                      IMEM_OF_AOS_MPASS 
      Print*, "The choosen path, one, two or mul", ONE_PASS_FONCE,
     &         TWO_PASS_FONCE, MUL_PAS_PONCE
      Write(6,*)
#endif
#endif
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
       I050=I040+MAX_BFNS_IRREP*MAX_BFNS_IRREP*IINTFP
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
       I050=I040+MAX_BFNS_IRREP*MAX_BFNS_IRREP*IINTFP
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
#ifndef _2HOT
       CALL LOAD1(CMO,ICORE(I000),ICORE(I010),ICORE(I020),
     &            ICORE(I030),ICORE(ISYMAO),NBAS,
     &            NFIRST,NSTART,NEND,ISIZE3,IOFFAO,
     &            IOFFI,NSIZE1,ILNBUF,ISPIN,LUINT,LAST)
#endif
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
#ifndef _2HOT
#ifdef _READ_ONLY
          CALL TRNSFORM_1INDX_1PASS(CMO,ICORE(I000),ICORE(I050),
     &                              ICORE(I040),KLPAIRS,NFIRST,
     &                              NSTART,NEND,ISPIN,NBAS,NMO,
     &                              NBASIS)
#endif
#endif
#ifdef _2HOT
          ONE_PASS_FONCE = .FALSE.
          TWO_PASS_FONCE = .TRUE.
          MUL_PASS_PONCE = .FLASE.
          IF (ONE_PASS_FONCE) THEN
             CALL READ_AO_INTS(ICORE(I000),ICORE(I020),ICORE(I030),
     &                         ICORE(ISYMAO),IOFFAO,NBAS,NSIZE1,
     &                         ILNBUF,LUINT,LAST,.TRUE.,.FALSE.)
C
             CALL COMPLETE_AO_INTS(ICORE(I000),ICORE(I050),KLPAIRS,
     &                             NBAS,NBASIS)
C
             CALL TRNSFORM_1INDX_FONCE(CMO,ICORE(I000),ICORE(I050),
     &                                 ICORE(I040),KLPAIRS,NFIRST,
     &                                 NSTART,NEND,ISPIN,NBAS,NMO,
     &                                 NBASIS)

             CALL INTRN1(CMO,ICORE(I000),ICORE(I010),
     &                   RHF,NBASIS,NBAS,NMO,POP(1,ISPIN),
     &                   POP(1,3-ISPIN),NSTART,NEND,ISIZE3,
     &                   IOFFMO,NSIZE,ISPIN,.TRUE.,ICORE(IREORD))
C
          ELSE IF (TWO_PASS_FONCE) THEN

             CALL READ_AO_INTS(ICORE(I000),ICORE(I020),ICORE(I030),
     &                         ICORE(ISYMAO),IOFFAO,NBAS,NSIZE1,
     &                         ILNBUF,LUINT,.FALSE.,.TRUE.,.FALSE.)

             CALL READ_AO_INTS(ICORE(I000),ICORE(I020),ICORE(I030),
     &                         ICORE(ISYMAO),IOFFAO,NBAS,NSIZE1,
     &                         ILNBUF,LUINT,.TRUE.,.FALSE.,.TRUE.)

             CALL TRNSFORM_1INDX_FONCE(CMO,ICORE(I000),ICORE(I050),
     &                                 ICORE(I040),KLPAIRS,NFIRST,
     &                                 NSTART,NEND,ISPIN,NBAS,NMO,
     &                                 NBASIS)

             CALL INTRN1(CMO,ICORE(I000),ICORE(I010),
     &                   RHF,NBASIS,NBAS,NMO,POP(1,ISPIN),
     &                   POP(1,3-ISPIN),NSTART,NEND,ISIZE3,
     &                   IOFFMO,NSIZE,ISPIN,.TRUE.,ICORE(IREORD))
C
          ELSE IF (MUL_PASS_PONCE) THEN


             CALL READ_AO_INTS(ICORE(I000),ICORE(I020),ICORE(I030),
     &                         ICORE(ISYMAO),IOFFAO,NBAS,NSIZE1,
     &                         ILNBUF,LUINT,.FALSE.,.TRUE.,.FALSE.)

             CALL READ_AO_INTS(ICORE(I000),ICORE(I020),ICORE(I030),
     &                         ICORE(ISYMAO),IOFFAO,NBAS,NSIZE1,
     &                         ILNBUF,LUINT,.TRUE.,.FALSE.,.TRUE.)
C 
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
             

             CALL TRNSFORM_1INDX_PONCE(CMO,ICORE(I000),ICORE(I050),
     &                                 ICORE(I040),KLPAIRS,NFIRST,
     &                                 NSTART,NEND,ISPIN,NBAS,NMO,
     &                                 IMEM_OF_AOS_MPASS)
             
          ENDIF 
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
