      SUBROUTINE L2INL1_LEOM(ICORE,MAXCOR,IUHF)
C
C  THIS ROUTINE IS THE DRIVER FOR THE CONTRIBUTION OF THE DOUBLES
C  TO THE SINGLES. THE CONTRIBUTION IS ACTUALLY VERY SIMILAR TO
C  THE FOURTH-ORDER SINGLES CONTRIBUTION
C
CEND
C
C CODED AUGUST/90  JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LAMBDA
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      INTEGER DIRPRD,POP1,POP2,VRT1,VRT2
C
      DIMENSION ICORE(MAXCOR)
C
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      COMMON/SYM/POP1(8),POP2(8),VRT1(8),VRT2(8),NT(2),NF1(2),NF2(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/FLAGS/IFLAGS(100)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
C
      DATA ONE /1.0D0/
C
      LAMBDA=.TRUE.
C
      NOCCA=NOCCO(1)
      NOCCB=NOCCO(2)
      NVRTA=NVRTO(1)
      NVRTB=NVRTO(2)
      NWAA=NT(1)
      NWBB=NT(2)
C
      ITAA=MAXCOR+1-NWAA*IINTFP
      MXCOR=MAXCOR+1-NWAA*IINTFP
      IF(IUHF.EQ.0) THEN
       ITBB=ITAA
      ELSE
       ITBB=ITAA-NWBB*IINTFP
       MXCOR=MXCOR-NWBB*IINTFP
      ENDIF
C
      LSTOFF=0
      IF(CCSD) THEN
       LSTOFF=100
       IF(IFLAGS(3).EQ.2) LSTOFF=200
      ENDIF

C
C  GET LAMBDA 1 INCREMENTS
C
      CALL GETLST(ICORE(ITAA),1,1,1,3,90)

      IF(IUHF.EQ.1) CALL GETLST(ICORE(ITBB),1,1,1,4,90)
C
C  WITHIN THE SPIN ADAPTED RHF CODE NO CALL TO T2T1AA2 IS NECESSARY
C
      IF(IUHF.EQ.1) THEN
       CALL T2T1AA2(ICORE(ITAA),ICORE,MXCOR,POP1,VRT1,1,LAMBDA,LSTOFF)
C
       CALL T2T1AA2(ICORE(ITBB),ICORE,MXCOR,POP2,VRT2,2,LAMBDA,LSTOFF)
      ENDIF
C
      CALL T2T1AB2(ICORE(ITAA),ICORE,MXCOR,POP1,POP2,VRT1,VRT2,1,
     &             IUHF,LAMBDA,LSTOFF)
C
      IF(IUHF.EQ.1) THEN
      CALL T2T1AB2(ICORE(ITBB),ICORE,MXCOR,POP2,POP1,VRT2,VRT1,2,
     &             IUHF,LAMBDA,LSTOFF)
      ENDIF
C
C  WITHIN THE SPIN ADAPTED RHF CODE NO CALL TO T2T1AA1 IS NECESSARY
C
      IF(IUHF.EQ.1) THEN
       CALL T2T1AA1(ICORE(ITAA),ICORE,MXCOR,POP1,VRT1,1,LAMBDA,LSTOFF)
C
       CALL T2T1AA1(ICORE(ITBB),ICORE,MXCOR,POP2,VRT2,2,LAMBDA,LSTOFF)
      ENDIF
C
      CALL T2T1AB1(ICORE(ITAA),ICORE,MXCOR,POP1,POP2,VRT1,VRT2,1,
     &             IUHF,LAMBDA,LSTOFF)
C
      IF(IUHF.EQ.1) THEN
      CALL T2T1AB1(ICORE(ITBB),ICORE,MXCOR,POP2,POP1,VRT2,VRT1,2,
     &             IUHF,LAMBDA,LSTOFF)
      ENDIF
C
C  FOR CCSD ADD HERE F(I,A) TO THE T1 INCREMENTS
C
      IF(CCSD) THEN
C
C  ADD F(A,I) TO Z(A,I)
C
       CALL GETLST(ICORE(1),1,1,1,1,93)
       CALL SAXPY(NWAA,ONE,ICORE,1,ICORE(ITAA),1)
C
C  FOR UHF ONLY, ADD F(a,i) TO Z(a,i) 
C
       IF(IUHF.EQ.1) THEN
        CALL GETLST(ICORE,1,1,1,2,93)
        CALL SAXPY(NWBB,ONE,ICORE,1,ICORE(ITBB),1)
       ENDIF
      ENDIF
C
C     CALCULATE NEW AMPLITUDES
C
      CALL NEWL1_leom(ICORE(ITAA),ICORE,MXCOR,1)
      CALL PUTLST(ICORE(ITAA),1,1,1,3,90)
      IF(IUHF.NE.0) THEN
       CALL NEWL1_leom(ICORE(ITBB),ICORE,MXCOR,2)
       CALL PUTLST(ICORE(ITBB),1,1,1,4,90)
      ENDIF
      RETURN
      END
