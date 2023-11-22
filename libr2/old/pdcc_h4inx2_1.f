











      SUBROUTINE PDCC_H4INX2_1(ICORE,MAXCOR,IUHF)
C
C     SUBROUTINE CALLS PDCC_H4X2ALL FOR THE THREE SPIN CASES
C     
CEND
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      DIMENSION ICORE(MAXCOR)
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
C
      IBOT=1
      IF(IUHF.EQ.0) IBOT=3
      Write(6,*) "Entering PDCC_H4INX2_1"
C
      DO 100 ISPIN=1,3
       CALL PDCC_H4X2ALL_1(ICORE,MAXCOR,ISPIN,IUHF,.TRUE.)
100   CONTINUE
C
      RETURN
      END
