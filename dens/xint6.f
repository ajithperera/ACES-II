       SUBROUTINE XINT6(XIA,ICORE,MAXCOR,IUHF)     
C
C
C THIS ROUTINE CALCULATES THE SIXTH TERM OF THE OCCUPIED-VIRTUAL BLOCK
C OF THE INTERMEDIATE XIA
C
C THE GENERAL FORMULA IS
C
C 1. TERM : + 1/2 SUM E,F,G G(EF,GI) <EF//GA>
C
C          SPIN TYPES : AA   AAAA AAAA    (UHF)
C                            ABAB ABAB    (UHF + RHF)
C                       BB   BBBB BBBB    (UHF)
C                            BABA BABA    (UHF)
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
      DIMENSION ICORE(MAXCOR),XIA(1)
      COMMON /SYM2/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /METH/ MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC 
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &DIRPRD(8,8)
      COMMON /SYMPOP2/ IRPDPD(8,22)
      COMMON /SYMPOP/ IRP_DM(8,22),ISYTYP(2,500),NTOT(18)
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
        IOFF=NTAA+1
       ENDIF
       IF(IUHF.EQ.1) THEN
C
C      AA AND BB SPIN CASES
C
C  LISTG :    GAMMA AMPLITUDES
C  LISTW :    INTEGRALS 
C  FACT :     PREFACTOR
C
C
       LISTW=230+ISPIN
       LISTG=126+ISPIN
C
       FACT=ONE
C
C LOOP OVER IRREPS OF EF BLOCK (THE SAME IRREPS AS THE AM AND IM BLOCKS 
C HAVE
C
       DO 100 IRREP=1,NIRREP
C
C DETERMINE LENGTH OF EXPANDED OCCUPIED-OCCUPIED BLOCK
C
        NVRTSQ=0 
        DO 110 IRREPJ=1,NIRREP
         NVRTSQ=NVRTSQ+VRT(IRREPJ,ISPIN)
     &                    *VRT(DIRPRD(IRREPJ,IRREP),ISPIN)
110     CONTINUE 
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTG))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTG))
        DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
        NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))  
        I001=1
        I002=I001+IINTFP*DISSYT*NUMSYT
        IF(MIN(NUMSYT,DISSYT,NUMSYW,DISSYW).NE.0) THEN
         I003=I002+3*IINTFP*MAX(DISSYT,NUMSYT,NUMSYW,DISSYW)
         MAXSIZE=(MXCOR-I003)/IINTFP
         IF(MAXSIZE.GT.NVRTSQ) THEN
C  
C         IN CORE VERSION
C
       CALL XIA6AA(ICORE(I001),ICORE(I003),ICORE(I003),MAXSIZE,
     &             XIA(IOFF),FACT,ISPIN,POP(1,ISPIN),VRT(1,ISPIN),
     &             DISSYT,NUMSYT,DISSYW,NUMSYW,NVRTSQ,LISTG,LISTW,
     &             IRREP,ICORE(I002))
         ELSE
          CALL INSMEM('XIA6AA',MAXSIZE,NVRTSQ)
         ENDIF
        ELSE
        ENDIF 
100    CONTINUE
       ENDIF
C
C       AB SPIN CASE
C
       LISTW=233
       LISTG=128+ISPIN+1-IUHF
       FACT=ONE
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
        I002=I001+IINTFP*NUMSYT*DISSYT
        IF(MIN(NUMSYT,DISSYT,NUMSYW,DISSYW).NE.0) THEN
         I003=I002+3*IINTFP*MAX(NUMSYT,DISSYT,NUMSYW,DISSYW)
         MAXSIZE=(MXCOR-I003)/IINTFP
         IF(MAXSIZE.GT.DISSYW) THEN
C
C         IN CORE VERSION
C
          CALL XIA6AB(ICORE(I001),ICORE(I003),MAXSIZE,XIA(IOFF),
     &                FACT,ISPIN,POP(1,ISPIN),POP(1,3-ISPIN),
     &                VRT(1,ISPIN),VRT(1,3-ISPIN),DISSYT,NUMSYT,
     &                DISSYW,NUMSYW,LISTG,LISTW,IRREP,ICORE(I002),IUHF)
         ELSE
          CALL INSMEM('XIA6AB',MAXSIZE,NVRTSQ)
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