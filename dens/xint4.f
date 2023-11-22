C
C THIS ROUTINE CALCULATES THE SECOND TERM OF THE OCCUPIED-VIRTUAL BLOCK
C OF THE INTERMEDIATE I
C
C THE GENERAL FORMULA IS
C
C + 2  SUM M,N,O   G(IM,NO) <AM//NO>
C
C          SPIN TYPES : AA   AAAA AAAA
C                            ABAB ABAB
C                       BB   BBBB BBBB
C                            BABA BABA
C
C ACTUALLY 4*G(IJ,KL) IS STORED, SO THAT ALL FACTORS REDUCE TO ONE
C
C THIS SUBROUTINE USES EXPLICITELY SYMMETRY
C
C CODED  AUGUST/90  JG
C 
      SUBROUTINE XINT4(XIA,ICORE,MAXCOR,IUHF)     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYT,DISSYW,POP,VRT
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      DIMENSION ICORE(MAXCOR),XIA(1)
      COMMON /SYM2/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /METH/ MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC 
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON /SYMPOP2/ IRPDPD(8,22)
      COMMON /SYMPOP/ IRP_DM(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /CONTROL/ IPRNT,IXXX,IXXX2
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
       LISTW=6+ISPIN + ISHIFT 
       LISTG=110+ISPIN
       FACT=ONE
C
C LOOP OVER IRREPS OF NO BLOCK (THE SAME IRREPS AS THE IN AND MN BLOCKS 
C HAVE
C
       DO 100 IRREP=1,NIRREP
C
C DETERMINE LENGTH OF EXPANDED OCCUPIED-OCCUPIED BLOCK
C
        NOCCSQ=0
        DO 101 IRREPJ=1,NIRREP
         IRREPI=DIRPRD(IRREP,IRREPJ)
         NOCCSQ=NOCCSQ+POP(IRREPJ,ISPIN)*POP(IRREPI,ISPIN)
101     CONTINUE
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTG))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTG))
        DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
        NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))  
        I001=1
        I002=I001+IINTFP*NOCCSQ*DISSYT
        I003=I002+IINTFP*NUMSYW*DISSYW
        IF(MIN(NUMSYT,DISSYT,NUMSYW,DISSYW).NE.0) THEN
         I004=I003+3*IINTFP*MAX(NUMSYW,NUMSYT,DISSYW,DISSYT)
         IF(MXCOR.GT.I004) THEN
C  
C         IN CORE VERSION
C
          CALL XIA4AA(ICORE(I001),ICORE(I002),XIA(IOFF),FACT,
     &                ISPIN,POP(1,ISPIN),VRT(1,ISPIN),DISSYT,NUMSYT,
     &                DISSYW,NUMSYW,LISTG,LISTW,IRREP,ICORE(I003))
         ELSE
          STOP 'XIA4AA'
         ENDIF
        ELSE
        ENDIF 
100    CONTINUE
       ENDIF
C
C       AB SPIN CASE
C
       LISTW=8+ISPIN+1-IUHF + ISHIFT 
       LISTG=113
       FACT=ONE
C
C      LOOP OVER IRREPS.
C
       DO 200 IRREP=1,NIRREP
C
        DISSYT=IRPDPD(IRREP,ISYTYP(1,13))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,13))
        DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
        NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))
        I001=1
        I002=I001+IINTFP*NUMSYT*DISSYT
        I003=I002+IINTFP*MAX(NUMSYW*DISSYW,NUMSYT*DISSYT)
        IF(MIN(NUMSYT,DISSYT,NUMSYW,DISSYW).NE.0) THEN
         I004=I003+3*IINTFP*MAX(NUMSYT,DISSYT,NUMSYW,DISSYW)
         IF(MXCOR.GT.I004) THEN
C
C         IN CORE VERSION
C
          CALL XIA4AB(ICORE(I001),ICORE(I002),
     &                XIA(IOFF),FACT,ISPIN,POP(1,ISPIN),
     &                POP(1,3-ISPIN),VRT(1,ISPIN),VRT(1,3-ISPIN),
     &                DISSYT,NUMSYT,DISSYW,NUMSYW,LISTG,LISTW,
     &                IRREP,ICORE(I003),IUHF)
         ELSE
          STOP 'XIA4AB'
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
