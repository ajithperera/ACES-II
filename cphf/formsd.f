
      SUBROUTINE FORMSD(ICORE,MAXCOR)
C
C  THIS ROUTINE TRANSFORMS THE OVERLAP MATRIX DERIVATIVES
C  FROM THE AO TO THE MO BASIS
C
CEND
C
C CODED JAN/91 JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,POP,VRT
      LOGICAL SCF,NONHF,FIELD,GEOM,THIRD,MAGN,SPIN
C
      DIMENSION ICORE(MAXCOR)
      DIMENSION ICMO(2)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NJUNK(18)
      COMMON/BASSYM/NBAS(8),NBASIS,NBASSQ,NBASTT
      COMMON/METHOD/IUHF,SCF,NONHF
      COMMON/OFFSET/IOFFC(8)
      COMMON/PERT/NTPERT,NPERT(8),IPERT(8),IXPERT,IYPERT,IZPERT,
     &            IYZPERT,IXZPERT,IXYPERT,ITRANSX,ITRANSY,ITRANSZ,
     &            NUCIND
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/LSYM/NLENQ(8),NLENT(8)
      COMMON/PERT2/FIELD,GEOM,THIRD,MAGN,SPIN
C
C DETERMINE FACTOR WHICH CORRESPONDS TO SYMMETRIC (ANTISYMMETRIC) MATRICES
C
      IF(GEOM) THEN
       IANTI=0
      ELSE IF(MAGN) THEN
       IANTI=1
      ENDIF
C
C  ALLOCATE MEMORY FOR MO COEFFICIENTS
C
      ICMO(1)=1
      ICMO(2)=ICMO(1)+NBASSQ*IUHF*IINTFP
      ISCR=ICMO(2)+NBASIS*NBASIS*IINTFP
      ISD=ISCR+NBASIS*NBASIS*IINTFP
      LENSD=0
      DO 10 IRREP=1,NIRREP
        LENSD=MAX(LENSD,NLENT(IRREP))
10    CONTINUE
      IEND=ISD+IINTFP*NBASIS*NBASIS
      MXSCR=NBASIS*NBASIS
C
C GET EIGEN VECTORS AND SYMMETRY PACK THEM
C
      CALL GETREC(20,'JOBARC','SCFEVCA0',NBASIS*NBASIS*IINTFP,
     &            ICORE(ISCR))
      CALL SYMC(ICORE(ISCR),ICORE(ICMO(1)),NBASIS,NBAS,.FALSE.,1)
      IF(IUHF.EQ.1) THEN
       CALL GETREC(20,'JOBARC','SCFEVCB0',
     &             NBASIS*NBASIS*IINTFP,ICORE(ISCR))
       CALL SYMC(ICORE(ISCR),ICORE(ICMO(2)),NBASIS,NBAS,.FALSE.,2)
      ENDIF
C
C FILL IOFFC
C
      IOFFC(1)=1
      DO 1 IRREP=1,NIRREP-1
       IOFFC(IRREP+1)=IOFFC(IRREP)+NBAS(IRREP)*NBAS(IRREP)
1     CONTINUE 
C
C ALLOCATE MEMORY FOR THE OVERLAP DERIVATIVES IN THE AO BASIS
C
C
C ALLOCATE MEMORY FOR THE TRANSFORMED OVERLAP MATRIX DERIVATIVES
C
      ISDEXP=IEND
      ISDMO=ISDEXP+NBASIS*NBASIS*IINTFP
C
C  LOOP OVER ALL IRREPS
C
      DO 100 IRREPP=1,NIRREP
C
C  CHECK IF THERE ANY PERTURBATIONS BELONGING TO THIS IRREP
C 
       IF(GEOM) THEN
C
C NUMBER OF GEOMETRICAL PERTURBATIONS
C 
        NP=NPERT(IRREPP)
       ELSE IF(MAGN) THEN
C
C NUMBER OF MAGNETIC FIELD COMPONENTS
C
        NP=0
        IF(IYZPERT.EQ.IRREPP) NP=NP+1 
        IF(IXZPERT.EQ.IRREPP) NP=NP+1 
        IF(IXYPERT.EQ.IRREPP) NP=NP+1
       ENDIF
C 
       IF(NP.EQ.0) GO TO 100
C
C  LOOP NOW OVER ALL PERTURBATIONS IN THIS IRREP
C
       DO 150 IP=1,NP
C
C GET TRIANGULAR OVERLAP DERIVATIVES IN AO BASIS
C
       CALL GETLST(ICORE(ISD),IP,1,1,IRREPP,101)
C
C EXPAND TO A FULL SYMMETRY PACKED SQUARE MATRIX
C
       CALL MATEXP(IRREPP,NBAS,ICORE(ISD),ICORE(ISDEXP),IANTI)
       NOCCSQ=IRPDPD(IRREPP,21)
       NVRTSQ=IRPDPD(IRREPP,19)
       NOVSQ=IRPDPD(IRREPP,9)
C
C  COMPUTE THE OCCUPIED-OCCUPIED BLOCK
C
       IOFF=ISDMO
       CALL FORMIJ(IRREPP,ICORE(ISDEXP),ICORE(ICMO(1)),NBASSQ,
     &             ICORE(ISCR),MXSCR,IUHF,ICORE(IOFF),NOCCSQ,
     &             .FALSE.)
C
       CALL PUTLST(ICORE(IOFF),IP,1,1,IRREPP,170)
       IF(IUHF.EQ.1) THEN
        IOFF2=IOFF+IINTFP*NOCCSQ
        CALL PUTLST(ICORE(IOFF2),IP,1,1,IRREPP,171)
       ENDIF
C
C  COMPUTE THE VIRTUAL-VIRTUAL BLOCK
C
       IOFF=IOFF+IINTFP*(IUHF*NOCCSQ+IRPDPD(IRREPP,22))
       CALL FORMAB(IRREPP,ICORE(ISDEXP),ICORE(ICMO(1)),NBASSQ,
     &             ICORE(ISCR),MXSCR,IUHF,ICORE(IOFF),NVRTSQ,
     &             .FALSE.)
C
       CALL PUTLST(ICORE(IOFF),IP,1,1,IRREPP,172)
       IF(IUHF.EQ.1) THEN
        IOFF2=IOFF+IINTFP*NVRTSQ
        CALL PUTLST(ICORE(IOFF2),IP,1,1,IRREPP,173)
       ENDIF
C
C  COMPUTE THE VIRTUAL-OCCUPIED BLOCK
C
       IOFF=IOFF+IINTFP*(IUHF*NVRTSQ+IRPDPD(IRREPP,20))
       CALL FORMAI(IRREPP,ICORE(ISDEXP),ICORE(ICMO(1)),NBASSQ,
     &             ICORE(ISCR),MXSCR,IUHF,ICORE(IOFF),NOVSQ,
     &             .FALSE.)
C
       CALL PUTLST(ICORE(IOFF),IP,1,1,IRREPP,174)
       IF(IUHF.EQ.1) THEN
        IOFF2=IOFF+IINTFP*NOVSQ
        CALL PUTLST(ICORE(IOFF2),IP,1,1,IRREPP,175)
       ENDIF
C
150    CONTINUE
C
100   CONTINUE
C
C ALL DONE, RETURN 
C
      RETURN
      END
