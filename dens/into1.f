       SUBROUTINE INTO1(AIOO,ICORE,MAXCOR,IUHF)     
C
C
C THIS ROUTINE CALCULATES THE FIRST TERM OF THE OCCUPIED-OCCUPIED BLOCK
C OF THE INTERMEDIATE I
C
C THE GENERAL FORMULA IS
C
C 1. TERM :  - SUM M,E,F   G(JM,EF) <IM//EF>
C
C          SPIN TYPES : AA   AAAA AAAA
C                            ABAB ABAB
C                       BB   BBBB BBBB
C                            BABA BABA
C
C  NOTE THAT FOR MBPT(2) GAMMA IS EXPLICITELY GIVEN BY
C
C   2 G(JM,EF) =  T1(JM,EF)
C
C  FOR MBPT(3) GAMMA IS EXPLICITELY
C
C   2G(JM,EF) =  T2(JM,EF)
C
C  FOR MBPT(4) GAMMA IS GIVEN AS
C
C   2G(JM,EF) = T3(JM,EF) + X(JM,EF)
C
C  FOR CCD GAMMA IS GIVEN AS
C
C   4G(JM,EF) = T(JM,EF)+L(JM,EF)+X(JM,EF)
C
C NOTE THAT IN RHF, A SPIN ADAPTED CODE IS USED
C
C THIS SUBROUTINE USES EXPLICITELY SYMMETRY
C
CEND
C
C CODED  JULY/90  JG
C 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYT,POP,VRT
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      LOGICAL DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,TRIP2,GABCD,
     &        RELAXED,TRULY_NONHF
      DIMENSION ICORE(MAXCOR),AIOO(1)
      COMMON/SYM2/POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC 
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/DERIV/DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,TRIP2,
     &             GABCD,RELAXED,TRULY_NONHF
      COMMON/SYMPOP2/IRPDPD(8,22)
      COMMON/SYMPOP/IRP_DM(8,22),ISYTYP(2,500),NTOT(18) 
      COMMON /SHIFT/ ISHIFT 
C
      DATA ONE,ONEM /1.0D+0,-1.0D+0/
C
      MXCOR=MAXCOR
C
       DO 1000 ISPIN=1,IUHF+1

       IF(ISPIN.EQ.1) THEN
        IOFF=1
       ELSE
        IOFF=NF1AA+1
       ENDIF
C
C  SKIP FIRST PART IN RHF
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
       LISTW=13+ISPIN + ISHIFT 
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
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTW))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTW))  
        I001=1
        I002=I001+IINTFP*NOCC2SQ*DISSYT
        I003=I002+IINTFP*NOCC2SQ*DISSYT
        IF(MIN(NUMSYT,DISSYT).NE.0) THEN
         IF(I003.LT.MXCOR) THEN
C  
C         IN CORE VERSION
C
          CALL IO1AA(ICORE(I001),ICORE(I002),AIOO(IOFF),FACT,ISPIN,
     &               POP(1,ISPIN),VRT(1,ISPIN),DISSYT,NUMSYT,LISTG,
     &               LISTW,IRREP,1)
         ELSE
          CALL INSMEM('I01AA',I003,MXCOR)
         ENDIF
        ELSE
        ENDIF 
100    CONTINUE
       ENDIF
C
C       AB SPIN CASE
C
       LISTW=16 + ISHIFT 
       IF(MBPT2) THEN
        LISTG=46
        FACT=ONEM
       ELSE IF(MBPT3.AND.(.NOT.ROHF)) THEN
        LISTG=63
        FACT=ONEM
       ELSE
        LISTG=116
        FACT=ONEM
       ENDIF
C
C      LOOP OVER IRREPS.
C
       DO 200 IRREP=1,NIRREP
C
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTW))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTW))
        I001=1
        I002=I001+IINTFP*NUMSYT*DISSYT
        I003=I002+IINTFP*NUMSYT*DISSYT
        IF(MIN(NUMSYT,DISSYT).NE.0) THEN
         I004=I003+3*IINTFP*MAX(NUMSYT,DISSYT)
         IF(I004.LT.MXCOR) THEN
C
C         IN CORE VERSION
C
          CALL IO1AB(ICORE(I001),ICORE(I002),AIOO(IOFF),
     &               FACT,ISPIN,POP(1,ISPIN),POP(1,3-ISPIN),
     &               VRT(1,ISPIN),VRT(1,3-ISPIN),DISSYT,NUMSYT,
     &               LISTG,LISTW,IRREP,ICORE(I003),IUHF,1)
         ELSE
          CALL INSMEM('IO1AB',I004,MXCOR)
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
