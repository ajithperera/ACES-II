C
C THIS ROUTINE CALCULATES THE FIFTH TERM OF THE OCCUPIED-VIRTUAL BLOCK
C OF THE INTERMEDIATE X
C
C THE GENERAL FORMULA IS
C
C   SUM N,E,F <AE//NF> G(IE,NF)
C
C          SPIN TYPES : AA   AAAA AAAA
C                            ABAB ABAB
C                            ABBA ABBA
C                       BB   BBBB BBBB
C                            BABA BABA
C                            BAAB BAAB
C
C THIS SUBROUTINE USES EXPLICITELY SYMMETRY
C
C NOTE THAT THE GAMMA AMPLITUDES HAVE TO BE RESORTED BEFORE
C CALLING XINT5. THIS IS DONE IN SORTGAM.
C
C CODED  AUGUST/90  JG
C 
      SUBROUTINE XINT5(XIA,ICORE,MAXCOR,IUHF)     
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
C
CSSS       call dzero(xia,ntaa)
CSSS       call checksum("xia0:",xia,ntaa)
       IF(IUHF.NE.0) THEN
C
C      AA AND BB SPIN CASES
C
C  LISTG :    GAMMA AMPLITUDES
C  LISTW :    INTEGRALS 
C  FACT :     PREFACTOR
C
C
       LISTW=26+ISPIN + ISHIFT 
       LISTG=122+ISPIN
C
C NOTE SINCE WE ARE USING <EA//NF> INSTEAD OF <AE//NF> SWITCH THE SIGN
C
       FACT=ONE
C
C LOOP OVER IRREPS OF EF BLOCK (THE SAME IRREPS AS THE IN AND MN BLOCKS 
C HAVE
C
       DO 100 IRREP=1,NIRREP
C
C DETERMINE LENGTH OF EXPANDED VIRTUAL-VIRTUAL BLOCK
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
        I002=I001+IINTFP*NUMSYT*DISSYT
        I003=I002+IINTFP*NUMSYW*NVRTSQ
        I004=I003+3*IINTFP*MAX(NUMSYT,NUMSYW,DISSYT,DISSYW)
        IF(MIN(NUMSYT,DISSYT,NUMSYW,DISSYW).NE.0) THEN
         IF(MXCOR.GE.I004) THEN
C  
C         IN CORE VERSION
C
          CALL XIA5AA(ICORE(I001),ICORE(I002),XIA(IOFF),FACT,
     &                ISPIN,POP(1,ISPIN),VRT(1,ISPIN),DISSYT,NUMSYT,
     &                DISSYW,NUMSYW,LISTG,LISTW,IRREP,ICORE(I003))
         ELSE
          CALL INSMEM('XIA5AA',I004,MXCOR)
         ENDIF
        ELSE
        ENDIF 
100    CONTINUE
C
       ENDIF
CSSS       call checksum("xia1:",xia,ntaa)
CSSS       call dzero(xia,ntaa)
C
C       AB SPIN CASE
C
       LISTW=28+ISPIN+1-IUHF + ISHIFT 
       LISTG=127-ISPIN-1+IUHF
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
        I003=I002+IINTFP*MAX(NUMSYW*DISSYW,NUMSYT*DISSYT)
        IF(MIN(NUMSYT,DISSYT,NUMSYW,DISSYW).NE.0) THEN
         I004=I003+3*IINTFP*MAX(NUMSYT,DISSYT,NUMSYW,DISSYW)
         IF(MXCOR.GE.I004) THEN 
C
C         IN CORE VERSION
C
          CALL XIA5AB(ICORE(I001),ICORE(I002),
     &                XIA(IOFF),FACT,ISPIN,POP(1,ISPIN),
     &                POP(1,3-ISPIN),VRT(1,ISPIN),VRT(1,3-ISPIN),
     &                DISSYT,NUMSYT,DISSYW,NUMSYW,LISTG,LISTW,
     &                IRREP,ICORE(I003),IUHF)
         ELSE
          CALL INSMEM('XIA5AB',I004,MXCOR)
         ENDIF
        ELSE
C
C
        ENDIF
200   CONTINUE
CSSS       call checksum("xia2:",xia,ntaa)
CSSS       call dzero(xia,ntaa)
C
C       BA SPIN CASE
C
       LISTW=31-ISPIN + ISHIFT 
       LISTG=119-ISPIN
C
       FACT=ONE
C
C      LOOP OVER IRREPS.
C
       DO 300 IRREP=1,NIRREP
C
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTG))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTG))
        DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
        NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))
        I001=1
        I002=I001+IINTFP*NUMSYT*DISSYT
        I003=I002+IINTFP*MAX(NUMSYW*DISSYW,NUMSYT*DISSYT)
        IF(MIN(NUMSYT,DISSYT,NUMSYW,DISSYW).NE.0) THEN
         I004=I003+3*IINTFP*MAX(NUMSYT,DISSYT,NUMSYW,DISSYW)
         IF(MXCOR.GE.I004) THEN
C
C         IN CORE VERSION
C
          CALL XIA5BA(ICORE(I001),ICORE(I002),
     &                XIA(IOFF),FACT,ISPIN,POP(1,ISPIN),
     &                POP(1,3-ISPIN),VRT(1,ISPIN),VRT(1,3-ISPIN),
     &                DISSYT,NUMSYT,DISSYW,NUMSYW,LISTG,LISTW,
     &                IRREP,ICORE(I003),IUHF)
         ELSE
          CALL INSMEM('XIA5BA',I004,MXCOR) 
         ENDIF
        ELSE
C
C
        ENDIF
300   CONTINUE
CSSS       call checksum("xia3:",xia,ntaa)
1000  CONTINUE
C
      RETURN

      END