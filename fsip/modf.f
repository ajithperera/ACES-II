      SUBROUTINE MODF(ICORE,MAXCOR,IUHF,ITYPE)
C
C THIS ROUTINE ADDS OR SUBTRACTS THE DIAGONAL ELEMENTS FROM
C THE ONE-PARTICLE PART OF THE EFFECTIVE HAMILTONIAN
C
C ITYPE > 0 - f(pp) added
C ITYPE < 0 - f(pp) subtracted
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,FACT
      CHARACTER*8 LABEL(2)
      DIMENSION ICORE(MAXCOR)
      COMMON/SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      DATA ONE /1.0/
      DATA ZILCH /0.0/
      DATA ONEM/-1.0/
      DATA LABEL /'SCFEVALA','SCFEVALB'/
C
      NBAs=NOCCO(1)+NVRTO(1)
C
      FACT=ZILCH
      IF(ITYPE.GT.0)FACT=ONE
      IF(ITYPE.LT.0)FACT=ONEM
C
      DO 10 ISPIN=1,1+IUHF
       I000=1
       I010=I000+IINTFP*NBAS
       I020=I010+IINTFP*NFMI(ISPIN)
       I030=I020+IINTFP*NFEA(ISPIN)
       CALL GETREC(20,'JOBARC',LABEL(ISPIN),IINTFP*NBAS,ICORE(I000))
       CALL GETLST(ICORE(I010),1,1,1,ISPIN,91)
       CALL GETLST(ICORE(I020),1,1,1,ISPIN,92)
       IOFFO=I000
       IOFFV=NOCCO(ISPIN)*IINTFP+I000
       IOFFMI=I010
       IOFFEA=I020 
       DO 20 IRREP=1,NIRREP
        NOCC=POP(IRREP,ISPIN)
        NVRT=VRT(IRREP,ISPIN)
        NUMMI=NOCC*NOCC
        NUMEA=NVRT*NVRT
        CALL SAXPY(NOCC,FACT,ICORE(IOFFO),1,ICORE(IOFFMI),NOCC+1)
        CALL SAXPY(NVRT,FACT,ICORE(IOFFV),1,ICORE(IOFFEA),NVRT+1)
        IOFFO=IOFFO+NOCC*IINTFP
        IOFFV=IOFFV+NVRT*IINTFP
        IOFFMI=IOFFMI+NUMMI*IINTFP
        IOFFEA=IOFFEA+NUMEA*IINTFP
20     CONTINUE
       CALL PUTLST(ICORE(I010),1,1,1,ISPIN,91)
       CALL PUTLST(ICORE(I020),1,1,1,ISPIN,92)
10    CONTINUE
C
      RETURN
      END
