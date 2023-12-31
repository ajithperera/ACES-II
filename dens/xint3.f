       SUBROUTINE XINT3(XIA,ICORE,MAXCOR,IUHF)     
C
C
C THIS ROUTINE CALCULATES THE FIRST TERM OF THE OCCUPIED-OCCUPIED BLOCK
C OF THE INTERMEDIATE I
C
C THE GENERAL FORMULA IS
C
C 1. TERM :  SUM M,E,F   G(IM,EF) <AM//EF>
C
C          SPIN TYPES : AA   AAAA AAAA
C                            ABAB ABAB
C                       BB   BBBB BBBB
C                            BABA BABA
C
C  NOTE THAT FOR MBPT(2) GAMMA IS EXPLIITELY GIVEN BY
C
C  G(JM,EF) = 1/2 T1(JM,EF)
C
C  FOR MBPT(3) GAMMA IS EXPLICITELY
C
C  G(JM,EF) = 1/2 T2(JM,EF)
C
C  FOR MBPT(4) GAMMA IS GIVEN AS 
C
C  G(JM,EF) = 1/2 T3(JM,EF) + 1/2 X(JM,EF)
C
C  FOR CCD GAMMA IS GIVEN BY
C
C  G(JM,EF) = 1/4 T(JM,EF) + 1/4 T(JM,EF) + 1/4 X(JM,EF)
C
C
C THIS SUBROUTINE USES EXPLICITELY SYMMETRY
C
CEND
C
C CODED  JULY/90  JG
C 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYT,DISSYW,POP,VRT
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      LOGICAL DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,TRIP2,GABCD,
     &        RELAXED,TRULY_NONHF
      DIMENSION ICORE(MAXCOR),XIA(1)
      COMMON/SYM2/POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC 
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYMPOP2/IRPDPD(8,22)
      COMMON/SYMPOP/IRP_DM(8,22),ISYTYP(2,500),NTOT(18)
      COMMON/DERIV/DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,TRIP2,
     &             GABCD,RELAXED,TRULY_NONHF
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
        NT=NTAA
       ELSE
        IOFF=NTAA+1
        NT=NTBB
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
       LISTW=26+ISPIN + ISHIFT 
       IF(MBPT2) THEN 
        LISTG=43+ISPIN
        FACT=ONEM
       ELSE IF(MBPT3.AND.(.NOT.ROHF)) THEN
        LISTG=60+ISPIN
        FACT=ONEM
       ELSE
        LISTG=113+ISPIN
        FACT=ONEM
       ENDIF
C
C LOOP OVER IRREPS OF EF BLOCK (THE SAME IRREPS AS THE IN AND MN BLOCKS 
C HAVE
C
       DO 100 IRREP=1,NIRREP
C
C DETERMINE LENGTH OF EXPANDED OCCUPIED-OCCUPIED BLOCK
C
        NOCC2SQ=0 
        DO 110 IRREPJ=1,NIRREP
         NOCC2SQ=NOCC2SQ+POP(IRREPJ,ISPIN)
     &                    *POP(DIRPRD(IRREPJ,IRREP),ISPIN)
110     CONTINUE 
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTG))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTG))
        DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
        NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))  
        I001=1
        I002=I001+IINTFP*NOCC2SQ*DISSYT
        I003=I002+3*IINTFP*MAX(NUMSYT,NUMSYW,DISSYT,DISSYW)
        IF(MIN(NUMSYT,DISSYT,NUMSYW,DISSYW).NE.0) THEN
         MAXSIZE=(MXCOR-I003)/IINTFP
         IF(MAXSIZE.GT.DISSYW) THEN
C  
C         IN CORE VERSION
C
          CALL XIA3AA(ICORE(I001),ICORE(I003),MAXSIZE,XIA(IOFF),FACT,
     &                ISPIN,POP(1,ISPIN),VRT(1,ISPIN),DISSYT,NUMSYT,
     &                DISSYW,NUMSYW,LISTG,LISTW,IRREP,ICORE(I002))
         ELSE
          CALL INSMEM('XIA3AA',MAXSIZE,DISSYW)
         ENDIF
        ELSE
        ENDIF 
100    CONTINUE
       ENDIF
C
C       AB SPIN CASE
C
       LISTW=31-ISPIN + ISHIFT 
       IF(MBPT2) THEN
        LISTG=46
        FACT=ONE
       ELSE IF(MBPT3.AND.(.NOT.ROHF)) THEN
        LISTG=63
        FACT=ONE
       ELSE
        LISTG=116
        FACT=ONE
       ENDIF
C
C      LOOP OVER IRREPS.
C
       DO 200 IRREP=1,NIRREP
C
        DISSYT=IRPDPD(IRREP,ISYTYP(1,46))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,46))
        DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
        NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))
        I001=1
        I002=I001+IINTFP*NUMSYT*DISSYT
        IF(MIN(NUMSYT,DISSYT).NE.0) THEN
         I003=I002+3*IINTFP*MAX(NUMSYT,DISSYT,NUMSYW,DISSYW)
         MAXSIZE=(MXCOR-I003)/IINTFP
         IF(MAXSIZE.GT.DISSYW) THEN
C
C         IN CORE VERSION
C
          CALL XIA3AB(ICORE(I001),ICORE(I003),MAXSIZE,
     &                XIA(IOFF),FACT,ISPIN,POP(1,ISPIN),
     &                POP(1,3-ISPIN),VRT(1,ISPIN),VRT(1,3-ISPIN),
     &                DISSYT,NUMSYT,DISSYW,NUMSYW,LISTG,LISTW,
     &                IRREP,ICORE(I002),IUHF)
         ELSE
          CALL INSMEM('XIA3AB',MAXSIZE,DISSYW)
         ENDIF
        ELSE
C
C
        ENDIF
200   CONTINUE
C
1000  CONTINUE
C
      RETURN

      END
