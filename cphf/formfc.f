      SUBROUTINE FORMFC(ICORE,MAXCOR)
C
C  THIS ROUTINE TRANSFORMS THE FERMI-CONTACT INTEGRALS
C  FROM AO TO MO BASIS
C
CEND
C
C CODED JAN/91 JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,POP,VRT
      LOGICAL SCF,NONHF,ROHF,SEMI,GEOM,FIELD,THIRD,MAGN,SPIN
C
      DIMENSION ICORE(MAXCOR)
      DIMENSION ICMO(2)
C
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NJUNK(18)
      COMMON/BASSYM/NBAS(8),NBASIS,NBASSQ,NBASTT
      COMMON/METHOD/IUHF,SCF,NONHF
      COMMON/OPENSH/ROHF,SEMI
      COMMON/OFFSET/IOFFC(8)
      COMMON/PERT/NTPERT,NPERT(8),IPERT(8),IXPERT,IYPERT,IZPERT,
     &            IXYPERT,IXZPERT,IYZPERT,ITRANSX,ITRANSY,ITRANSZ,
     &            NUCIND
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/LSYM/NLENQ(8),NLENT(8)
      COMMON/PERT2/FIELD,GEOM,THIRD,MAGN,SPIN
C
      DATA ONEM/-1.D0/
C
      IANTI=0
C
C  ALLOCATE MEMORY FOR MO COEFFICIENTS
C
      ICMO(1)=1
      ICMO(2)=ICMO(1)+NBASSQ*IUHF*IINTFP
      ISCR=ICMO(2)+NBASIS*NBASIS*IINTFP
      IFD=ISCR+IINTFP*NBASIS*NBASIS
      IFDEXP=ISCR+IINTFP*NBASSQ
      IFDMO=IFDEXP+IINTFP*NBASSQ
      IEND=IFDMO+IINTFP*NBASSQ*(IUHF+1)
      IF(IEND.GE.MAXCOR) CALL ERREX
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
C  LOOP OVER ALL IRREPS
C
      DO 100 IRREPP=1,NIRREP
C
C  CHECK IF THERE ANY PERTURBATIONS BELONGING TO THIS IRREP
C 
       NP=NPERT(IRREPP)
C
       IF(NP.EQ.0) GO TO 100
C
C  LOOP OVER ALL PERTURBATIONS IN THIS IRREP
C
       DO 150 IP=1,NP
C
C READ IN THE FOCK MATRIX DERIVATIVE INTEGRALS
C
       CALL GETLST(ICORE(ISCR),IP,1,1,IRREPP,102)
C
       IFDEXP2=IFDEXP
       CALL MATEXP(IRREPP,NBAS,ICORE(ISCR),ICORE(IFDEXP),IANTI)
C
       NOCCSQ=IRPDPD(IRREPP,21)
       NVRTSQ=IRPDPD(IRREPP,19)
       NOVSQ=IRPDPD(IRREPP,9)
C
C  COMPUTE THE OCCUPIED-OCCUPIED BLOCK
C
       IOFF=IFDMO
       CALL FORMIJ(IRREPP,ICORE(IFDEXP),ICORE(ICMO(1)),NBASSQ,
     &             ICORE(ISCR),MXSCR,IUHF,ICORE(IOFF),NOCCSQ,
     &             .FALSE.)
C
       CALL PUTLST(ICORE(IOFF),IP,1,1,IRREPP,184)

       IF(IUHF.EQ.1) THEN
        CALL SSCAL(IRPDPD(IRREPP,22),ONEM,ICORE(IOFF+IINTFP*NOCCSQ),1)
        CALL PUTLST(ICORE(IOFF+IINTFP*NOCCSQ),IP,1,1,IRREPP,185)
       ENDIF
C
C  COMPUTE THE VIRTUAL-VIRTUAL BLOCK
C
       IOFF=IOFF+IINTFP*(IUHF*NOCCSQ+IRPDPD(IRREPP,22))
       CALL FORMAB(IRREPP,ICORE(IFDEXP),ICORE(ICMO(1)),NBASSQ,
     &             ICORE(ISCR),MXSCR,IUHF,ICORE(IOFF),NVRTSQ,
     &             .FALSE.)
C
       CALL PUTLST(ICORE(IOFF),IP,1,1,IRREPP,178)
C
       IF(IUHF.EQ.1) THEN
        CALL SSCAL(IRPDPD(IRREPP,20),ONEM,ICORE(IOFF+IINTFP*NVRTSQ),1)
        CALL PUTLST(ICORE(IOFF+IINTFP*NVRTSQ),IP,1,1,IRREPP,179)
       ENDIF
C
C  COMPUTE THE VIRTUAL-OCCUPIED BLOCK
C
       IOFF=IOFF+IINTFP*(IUHF*NVRTSQ+IRPDPD(IRREPP,20))
       CALL FORMAI(IRREPP,ICORE(IFDEXP),ICORE(ICMO(1)),NBASSQ,
     &             ICORE(ISCR),MXSCR,IUHF,ICORE(IOFF),NOVSQ,
     &             .FALSE.)
C
       CALL PUTLST(ICORE(IOFF),IP,1,1,IRREPP,180)
C
       IF(IUHF.EQ.1) THEN
        CALL SSCAL(IRPDPD(IRREPP,10),ONEM,ICORE(IOFF+IINTFP*NOVSQ),1)
        CALL PUTLST(ICORE(IOFF+IINTFP*NOVSQ),IP,1,1,IRREPP,181)
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
