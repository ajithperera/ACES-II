       SUBROUTINE INTO5(AIOO,ICORE,MAXCOR,IUHF)     
C
C
C THIS ROUTINE CALCULATES THE FIFTH TERM OF THE OCCUPIED-OCCUPIED BLOCK
C OF THE INTERMEDIATE I
C
C THE GENERAL FORMULA IS
C
C 1. TERM :  - SUM M,N SUM E <IM//NE> G(JM,NE)
C
C          SPIN TYPES : AA   AAAA AAAA    (UHF)
C                            ABAB ABAB    (UHF+RHF)
C                            ABBA ABBA    (UHF+RHF)
C                       BB   BBBB BBBB    (UHF)
C                            BABA BABA    (UHF)
C                            BAAB BAAB    (UHF)
C
C  FOR RHF A SPIN ADAPTED CODE IS USED SO ONLY ONE TERM
C  HAS TO BE CALCULATED.
C
C  THIS TERM IS ONLY REQUIRED FOR METHODS WHICH INCLUDE
C  SINGLE EXCITATION, E.G. MBPT(4), QCISD, CCSD
C
C THIS SUBROUTINE USES EXPLICITELY SYMMETRY
C
CEND
C
C CODED  AUGUST/90 JG
C 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYT,DISSYW,POP,VRT
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      DIMENSION ICORE(MAXCOR),AIOO(1)
      COMMON /SYM2/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /METH/ MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC 
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &DIRPRD(8,8)
      COMMON/SYMPOP2/IRPDPD(8,22)
      COMMON/SYMPOP/IRP_DM(8,22),ISYTYP(2,500),NTOT(18) 
      COMMON /SHIFT/ ISHIFT
C
      DATA ONE,ONEM /1.0D+0,-1.0D+0/
C
      MXCOR=MAXCOR
C
      DO 1000 ISPIN=1,IUHF+1 
C
       IF(ISPIN.EQ.1) THEN
        IOFF=1
       ELSE
        IOFF=NF1AA+1
       ENDIF
C
       IF(IUHF.EQ.1) THEN
C
C      AA AND BB SPIN CASES
C
C  LISTG :    GAMMA AMPLITUDES
C  LISTW :    INTEGRALS 
C  FACT :     PREFACTOR
C
C
       LISTW=6+ISPIN + ISHIFT 
       LISTG=106+ISPIN
C
       FACT=ONEM
C
C LOOP OVER IRREPS OF MN BLOCK (THE SAME IRREPS AS THE OI AND OA BLOCKS 
C HAVE
C
       DO 100 IRREP=1,NIRREP
C
C
        NOCCSQ=0
        DO 110 IRREPJ=1,NIRREP
         NOCCSQ=NOCCSQ+POP(IRREPJ,ISPIN)
     &            *POP(DIRPRD(IRREP,IRREPJ),ISPIN)
110     CONTINUE
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTG))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTG))
        DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
        NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))  
        I001=1
        I002=I001+IINTFP*NOCCSQ*NUMSYT
        I003=I002+IINTFP*NOCCSQ*MAX(NUMSYW,NUMSYT) 
        I004=I003+IINTFP*MAX(NUMSYW*DISSYW,NUMSYT*DISSYT) 
        IF(MIN(NUMSYT,DISSYT,NUMSYW,DISSYW).NE.0) THEN
         I005=I004+3*IINTFP*MAX(DISSYT,NUMSYT,NUMSYW,DISSYW)
         IF(I005.LT.MXCOR) THEN
C  
C         IN CORE VERSION
C
       CALL IO5AA(ICORE(I001),ICORE(I002),ICORE(I003),AIOO(IOFF),
     &             FACT,ISPIN,POP(1,ISPIN),VRT(1,ISPIN),DISSYT,
     &             NUMSYT,DISSYW,NUMSYW,LISTG,LISTW,IRREP,ICORE(I004))
         ELSE
          CALL INSMEM('I05AA',I005,MXCOR)
         ENDIF
        ELSE
        ENDIF 
100    CONTINUE
       ENDIF
C
C       AB SPIN CASE
C
       LISTW=11-ISPIN + ISHIFT 
       LISTG=111-ISPIN
       FACT=ONEM
C
C      LOOP OVER IRREPS.
C
       DO 200 IRREP=1,NIRREP
C
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTG))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTG))
        DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
        NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))
        I001=1
        I002=I001+IINTFP*MAX(NUMSYT*DISSYT,NUMSYW*DISSYW*(1-IUHF))
        I003=I002+IINTFP*NUMSYW*MAX(DISSYW,NUMSYT) 
        I004=I003+IINTFP*MAX(NUMSYT*DISSYT,NUMSYW*DISSYW)
        IF(MIN(NUMSYT,DISSYT,NUMSYW,DISSYW).NE.0) THEN
         I005=I004+3*IINTFP*MAX(NUMSYT,DISSYT,NUMSYW,DISSYW)
         IF(I005.LT.MXCOR) THEN
C
C         IN CORE VERSION
C
          CALL IO5AB(ICORE(I001),ICORE(I002),ICORE(I003),AIOO(1),
     &               AIOO(1+NF1AA*IUHF),FACT,ISPIN,POP(1,1),POP(1,2),
     &               VRT(1,1),VRT(1,2),DISSYT,NUMSYT,DISSYW,NUMSYW,
     &               LISTG,LISTW,IRREP,ICORE(I004),IUHF)
         ELSE
          CALL INSMEM('I05AB',I005,MXCOR)
         ENDIF
        ELSE
C
C
        ENDIF
200   CONTINUE
1000  CONTINUE
C
      RETURN

      END