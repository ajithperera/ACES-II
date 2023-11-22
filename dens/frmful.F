      SUBROUTINE FRMFUL(DOO,DVV,DVO,DFULL,
     &                  ISCR,MAXCOR,NBAS,IUHF,KEY)
C
C THIS ROUTINE ACCEPTS OCCUPIED-OCCUPIED, VIRTUAL-VIRTUAL AND
C  OCCUPIED-VIRTUAL SYMMETRY-PACKED SUBBLOCKS OF A TWO-INDEX
C  QUANTITY AND CONSTRUCTS THE FULL UNPACKED MATRIX.
C
C INPUT:
C       DOO   - THE SYMMETRY-PACKED OCCUPIED-OCCUPIED BLOCK 
C                OF THE DENSITY MATRIX.
C       DVV   - THE SYMMETRY-PACKED OCCUPIED-OCCUPIED BLOCK 
C                OF THE DENSITY MATRIX.
C       DVO   - THE SYMMETRY-PACKED OCCUPIED-OCCUPIED BLOCK 
C                OF THE DENSITY MATRIX.
C       ISPIN - THE SPIN CASE.
C       NBAS  - NUMBER OF BASIS FUNCTIONS.
C       MAXCOR- AVAILABLE LENGTH FOR ISCR.
C       KEY   - OPTION TO DECIDE WHETHER DENSITY (DEN), INTERMEDIATE
C               MATRIX (INT) OR SOME ADDITIONAL CONTRIBUTION TO THE
C               DENSITY (REL) IS FORMED
C
C SCRATCH:
C       ISCR  - USED FOR SYMMETRY VECTORS AND OTHER THINGS.
C                MUST BE OF LENGTH (1+IINTFP)*MAX(NVRT,NOCC)
C      
C OUTPUT:
C       DFULL - THE FULL UNPACKED M.O. DENSITY MATRIX.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER POP,VRT
      CHARACTER*1 SPTYPE(2)
      CHARACTER*3 KEY
CJDW 10/4/95. Added for dropped core triple logic.
      LOGICAL DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,TRIP2,GABCD,
     &        RELAXED,TRULY_NONHF
CJDW END
      DIMENSION DOO(1),DVV(1),DVO(1),DFULL(NBAS,NBAS),
     &          ISCR(MAXCOR)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
CJDW 10/4/95. Added for dropped core triple logic.
      COMMON /DERIV/ DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,
     &               TRIP2,GABCD,RELAXED,TRULY_NONHF
CJDW END
CJDW KKB stuff
C     COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
C     COMMON/SYM/POP(8,2),VRT(8,2),NDOVA,NDOVB,NDOOA,NDOOB,NDVVA,NDVVB
      COMMON /SYMPOP2/ IRPDPD(8,22)
      COMMON /SYM2/    POP(8,2),VRT(8,2),NDOVA,NDOVB,NDOOA,NDOOB,NDVVA,
     &                                                           NDVVB
      COMMON /DROPGEO/ NDRGEO
CJDW END
      DATA SPTYPE /'A','B'/
C
C A bug fix; In prior versions IONE does not get initilized for NONHF
C options. 06/06, Ajith Perera.

      DATA IONE /1/
C
      SCFOCC=2.0
      IF(IUHF.EQ.1)SCFOCC=1.0
      IOFFO=1
      IOFFV=1
      IOFFOV=1
C
C IF SECOND DERIVATIVES ARE CALCULATED, SAVE DENSITY MATRICES
C AND INTERMEDIATES ON GAMLAM
C
CJDW 10/4/95. Problems here for dropped core. In the KKB code (and
C             presumably earlier vintages of the standard) we have
C             "      IF(IFLAGS(3).GE.2.AND.KEY.EQ.'DEN') THEN"
C
C             But in the latest standard code we have
C             "      IF(KEY.EQ.'DEN') THEN"
C
C             The latter causes a crash in dropped core triples, so we
C             not call it if we have dropped core triples. But why does
C             it not
C             have a problem with CCSD ?? Because the lists already exist.
C             Moreover, when they were created, they had dropped core
C             dimensions, which makes sense, but now they need full
C             dimensions. So it makes sense.
C
      IF(NDRGEO.EQ.0.OR.(.NOT.TRIP1.AND..NOT.TRIP2))THEN
        IF(KEY.EQ.'DEN') THEN
          CALL SAVD(DOO,DOO(1+NDOOA),DVV,DVV(1+NDVVA),DVO,
     &              DVO(1+NDOVA),IUHF)
        ENDIF
      ENDIF
C
      IF(IFLAGS(3).GE.2.AND.KEY.EQ.'INT') THEN
       CALL SAVI(DOO,DOO(1+NDOOA),DVV,DVV(1+NDVVA),DVO,
     &           DVO(1+NDOVA),ISCR,IUHF)
      ENDIF

      DO 100 ISPIN=1,1+IUHF
C
       CALL ZERO(DFULL,NBAS*NBAS)
C
C FIRST, PUT THE OO, VV AND VO BLOCKS INTO PLACE.
C
       LENOO=NOCCO(ISPIN)*NOCCO(ISPIN)
       LENVV=NVRTO(ISPIN)*NVRTO(ISPIN)
       LENVO=NVRTO(ISPIN)*NOCCO(ISPIN)
#ifdef _DEBUG_LVLM
       Write(6,*) "Occ-Occ and Vrt-Vrt blocks of the D matrix"
       loo = NOCCO(ISPIN)
       lvv = NVRTO(ISPIN)
       lvo = LENVO
CSS       Write(6,*) "The DOO at entry"
cSS       Write(6,"(5(1x,F10.7))")(Doo(i),i=1,(Nocco(1)*(Nocco(1)+1))/2)
CSS       call output(doo(IOFFO), 1, loo, 1, loo, loo, loo, 1)
C       call output(dvv, 1, lvv, 1, lvv, lvv, lvv, 1)
C       call output(dvo, 1, lvo, 1, lv0, lvo, lvo, 1)
#endif
       I000=1
       I010=I000+LENOO*IINTFP
       I020=I010+LENOO
       IF(I020.GT.MAXCOR)CALL INSMEM('FILFUL',I020,MAXCOR)
       CALL PUTBLK(DFULL,DOO(IOFFO),'OO',ISCR(I000),ISCR(I010),NBAS,
     &            IRPDPD(1,20+ISPIN),LENOO,ISPIN)
       I000=1
       I010=I000+LENVV*IINTFP
       I020=I010+LENVV
       IF(I020.GT.MAXCOR)CALL INSMEM('FILFUL',I020,MAXCOR)
       CALL PUTBLK(DFULL,DVV(IOFFV),'VV',ISCR(I000),ISCR(I010),NBAS,
     &             IRPDPD(1,18+ISPIN),LENVV,ISPIN)
       I000=1
       I010=I000+LENVO*IINTFP
       I020=I010+LENVO
       IF(I020.GT.MAXCOR)CALL INSMEM('FILFUL',I020,MAXCOR)
       CALL PUTBLK(DFULL,DVO(IOFFOV),'VO',ISCR(I000),ISCR(I010),NBAS,
     &             IRPDPD(1,8+ISPIN),LENVO,ISPIN)
#ifdef _DEBUG_LVL0
       If (KEY .EQ. "DEN") THEN
          Write(6,"(a)") "Relaxed or response density matrix (no ref.)"
          Call output(Dfull,1,Nbas,1,Nbas,Nbas,Nbas,1)
       Elseif (KEY .EQ. "INT") THEN
          Write(6,"(a)") "Rotation gradients"
          Call output(Dfull,1,Nbas,1,Nbas,Nbas,Nbas,1)
       Endif
#endif 
C
C NOW PUT THE OV BLOCK INTO PLACE BY SYMMETRIZATION.
C
       DO 10 IROW=1,NOCCO(ISPIN)
*VOCL LOOP,NOVREC
CDIR$ IVDEP
        DO 20 ICOL=NOCCO(ISPIN)+1,NBAS
         DFULL(IROW,ICOL)=DFULL(ICOL,IROW)
20      CONTINUE
10     CONTINUE
       IF(IUHF.EQ.0) THEN
        CALL SSCAL(NBAS*NBAS,SCFOCC,DFULL,1)
       ENDIF
       IF(KEY.EQ.'DEN') THEN
        DO 30 IORB=1,NOCCO(ISPIN)
         DFULL(IORB,IORB)=DFULL(IORB,IORB)+SCFOCC
30      CONTINUE
CJDW KKB stuff
CJDW END
        IF(KEY.EQ.'DEN') THEN
         CALL PUTREC(20,'JOBARC','RELDENS'//SPTYPE(ISPIN),
     &               NBAS*NBAS*IINTFP,DFULL)
CJDW KKB stuff. ???????
         CALL PUTREC(20,'JOBARC','RELDSIZ'//SPTYPE(ISPIN),
     &               IONE,NBAS*NBAS*IINTFP)
CJDW END
        ENDIF
       ELSE IF(KEY.EQ.'INT') THEN
        CALL PUTREC(20,'JOBARC','IINTERM'//SPTYPE(ISPIN),
     &              NBAS*NBAS*IINTFP,DFULL)
CJDW KKB stuff ??
        CALL PUTREC(20,'JOBARC','IINTSIZ'//SPTYPE(ISPIN),
     &              IONE,NBAS*NBAS*IINTFP)
CJDW END
       ELSE IF(KEY.EQ.'REL') THEN
        CALL PUTREC(20,'JOBARC','NHFDENS'//SPTYPE(ISPIN),
     &              NBAS*NBAS*IINTFP,DFULL)
CJDW KKB stuff
        CALL PUTREC(20,'JOBARC','NHFDSIZ'//SPTYPE(ISPIN),
     &              IONE,NBAS*NBAS*IINTFP)
CJDW END
       ENDIF
       IOFFO=IOFFO+NDOOA
       IOFFV=IOFFV+NDVVA
       IOFFOV=IOFFOV+NDOVA
100   CONTINUE
      RETURN
      END 
