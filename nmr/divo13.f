       SUBROUTINE DIVO13(DIVO,ICORE,MAXCOR,IUHF,IP,ANTI)
C
C
C THIS ROUTINE CALCULATES THE FIRST TERM OF THE OCCUPIED-OCCUPIED BLOCK
C OF THE INTERMEDIATE DI
C
C THE GENERAL FORMULA IS
C
C 1. TERM : - SUM M,E,F   G(IM,EF) d<EF||AM>/dchi
C
C          SPIN TYPES : AA   AAAA AAAA
C                            ABAB ABAB
C                       BB   BBBB BBBB
C                            BABA BABA
C
C  NOTE THAT FOR MBPT(2) GAMMA IS EXPLICITELY GIVEN BY
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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYG,DISSYW,POP,VRT
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,ANTI
      DIMENSION ICORE(MAXCOR),DIVO(1)
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
      COMMON/DSYM/IRREPX,IPERT,NDT(2),NDF1(2),NDF2(2),
     &            IOFFIJ(8,2),IOFFAB(8,2),IOFFAI(8,2)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /METH/ MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD 
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
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
        IOFF=NDT(1)+1
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
       LISTW=326+ISPIN
       IF(MBPT2) THEN 
        LISTG=43+ISPIN
        FACT=ONE
       ELSE IF(MBPT3) THEN
        LISTG=60+ISPIN
        FACT=ONE
       ELSE
        LISTG=113+ISPIN
        FACT=ONE
       ENDIF
C
C LOOP OVER IRREPS OF EF BLOCK (THE SAME IRREPS AS THE IN AND MN BLOCKS 
C HAVE
C
       DO 100 IRREPR=1,NIRREP
C
        IRREPL=DIRPRD(IRREPX,IRREPR)
C
C DETERMINE LENGTH OF EXPANDED OCCUPIED-OCCUPIED BLOCK
C
        NO2SQ=0 
        DO 110 IRREPJ=1,NIRREP
         NO2SQ=NO2SQ+POP(IRREPJ,ISPIN)
     &                    *POP(DIRPRD(IRREPJ,IRREPR),ISPIN)
110     CONTINUE 
        DISSYG=IRPDPD(IRREPR,ISYTYP(1,43+ISPIN))
        NUMSYG=IRPDPD(IRREPR,ISYTYP(2,43+ISPIN))
        DISSYW=DISSYG
        NUMSYW=IRPDPD(IRREPL,ISYTYP(2,LISTW))  
        I001=1
        I002=I001+IINTFP*NO2SQ*DISSYG
        I003=I002+3*IINTFP*MAX(NUMSYG,NUMSYW,DISSYG,DISSYW)
        IF(MIN(NUMSYG,DISSYG,NUMSYW,DISSYW).NE.0) THEN
         MAXSIZE=(MXCOR-I003)/IINTFP
         IF(MAXSIZE.GT.DISSYW) THEN
C  
C         IN CORE VERSION
C
          CALL DIVO3AA(ICORE(I001),ICORE(I003),MAXSIZE,DIVO(IOFF),
     &                 FACT,ISPIN,POP(1,ISPIN),VRT(1,ISPIN),DISSYG,
     &                 NUMSYG,DISSYW,NUMSYW,LISTG,LISTW,IRREPL,
     &                 IRREPR,ICORE(I002))
         ELSE
          CALL INSMEM('DIVO3AA',MAXSIZE,DISSYW)
         ENDIF
        ELSE
        ENDIF 
100    CONTINUE
       ENDIF
C
C       AB SPIN CASE
C
       LISTW=331-ISPIN
       IF(MBPT2) THEN
        LISTG=46
        FACT=ONEM
       ELSE IF(MBPT3) THEN
        LISTG=63
        FACT=ONEM
       ELSE
        LISTG=116
        FACT=ONEM
       ENDIF
C
C      LOOP OVER IRREPS.
C
       DO 200 IRREPR=1,NIRREP
C
        IRREPL=DIRPRD(IRREPX,IRREPR)
C
        DISSYG=IRPDPD(IRREPR,ISYTYP(1,46))
        NUMSYG=IRPDPD(IRREPR,ISYTYP(2,46))
        DISSYW=DISSYG
        NUMSYW=IRPDPD(IRREPL,ISYTYP(2,LISTW))
        I001=1
        I002=I001+IINTFP*NUMSYG*DISSYG
        IF(MIN(NUMSYG,DISSYG).NE.0) THEN
         I003=I002+3*IINTFP*MAX(NUMSYG,DISSYG,NUMSYW,DISSYW)
         MAXSIZE=(MXCOR-I003)/IINTFP
         IF(MAXSIZE.GT.DISSYW) THEN
C
C         IN CORE VERSION
C
          CALL DIVO3AB(ICORE(I001),ICORE(I003),MAXSIZE,
     &                 DIVO(IOFF),FACT,ISPIN,POP(1,ISPIN),
     &                 POP(1,3-ISPIN),VRT(1,ISPIN),VRT(1,3-ISPIN),
     &                 DISSYG,NUMSYG,DISSYW,NUMSYW,LISTG,LISTW,
     &                 IRREPL,IRREPR,ICORE(I002),IUHF)
         ELSE
          CALL INSMEM('DIVO3AB',MAXSIZE,DISSYW)
         ENDIF
        ELSE
C
C
        ENDIF
200   CONTINUE
C
      call checksum('divo13  ',divo(ioff),ndt(ispin))
1000  CONTINUE
C
      RETURN

      END
