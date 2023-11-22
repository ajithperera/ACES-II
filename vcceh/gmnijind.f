C
      SUBROUTINE GMNIJIND(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF)
C
C This routine calculates the hole-hole ladder contribution
C to the doubles equation. Spin orbital formula for this 
C contribution is
C
C  + 1/2 SUM T(BET) * G(ALP)    ; formulas for G(ALP) are already given
C        M,N (AB,MN)  (MN,IJ)                  (MN,IJ)
C
C and calculated. Spin integrated formulas can be written as (here
C alpha and beta is implicit)
C
C UHF
C
C   Q(MN,IJ) = SUM M<N T(AB,MN)*G(MN,IJ)  [AAAA]
C    
C   Q(mn,ij) = SUM m<n T(ab,mn)*G(mn,ij)  [BBBB]
C           
C   Q(Mn,Ij) = SUM M,n T(Ab,Mn)*G(Mn,Ij)  [ABAB]
C
C RHF
C
C   Q(Mn,Ij) = SUM M,n T(Ab,Mn)*G(Mn,Ij)  [ABAB]
C               
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER DIRPRD,DISSYG,DISSYT,POP,VRT
      CHARACTER*2 SPCASE(3)
      DIMENSION ICORE(MAXCOR)
C
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYM/POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NJUNK(18)
      COMMON /FLAGS/ IFLAGS(100)
C   
C Common blocks used in the quadratic term
C 
      COMMON /QADINTI/ INTIMI, INTIAE, INTIME
      COMMON /QADINTG/ INGMNAA, INGMNBB, INGMNAB, INGABAA,
     &                 INGABBB, INGABAB, INGMCAA, INGMCBB,
     &                 INGMCABAB, INGMCBABA, INGMCBAAB,
     &                 INGMCABBA
      COMMON /APERTT1/ IAPRT1AA             
      COMMON /BPERTT1/ IBPRT1AA
      COMMON /APERTT2/ IAPRT2AA1, IAPRT2BB1, IAPRT2AB1, IAPRT2AB2, 
     &                 IAPRT2AB3, IAPRT2AB4, IAPRT2AA2, IAPRT2BB2,
     &                 IAPRT2AB5 
C
      DATA SPCASE /'AA', 'BB', 'AB'/
C     
      MXCOR = MAXCOR
C
      IF (IUHF .EQ. 0) THEN
         IBOT = 3
      ELSE
         IBOT = 1
      ENDIF
C     
      DO 1000 ISPIN = IBOT, 3
C     
         IF (ISPIN .LT. 3) THEN
C     
C AAAA or BBBB cases.
C     
            LISTT = ISPIN + (IAPRT2AA2 - 1)
            LISTG = ISPIN + (INGMNAA - 1)
            LISTZ = ISPIN + 60
C
            DO 100 IRREPMN = 1, NIRREP
               IRREPTMN = IRREPMN
               IRREPTAB = DIRPRD(IRREPTMN, IRREPX)
               IRREPGMN = IRREPTMN
               IRREPGIJ = DIRPRD(IRREPGMN, IRREPX)
               IRREPQIJ = IRREPGIJ
               IRREPQAB = IRREPQIJ
C
               DISSYG = IRPDPD(IRREPGMN, ISYTYP(1, LISTG))
               NUMSYG = IRPDPD(IRREPGIJ, ISYTYP(2, LISTG))
               DISSYT = IRPDPD(IRREPTAB, ISYTYP(1, LISTT))
               NUMSYT = IRPDPD(IRREPTMN, ISYTYP(2, LISTT))
C
               I000 = 1
               I010 = I000 + MAX(DISSYT*NUMSYT, DISSYT*NUMSYG)*IINTFP
               I020 = I010 + MAX(DISSYT*NUMSYT, DISSYT*NUMSYG)*IINTFP
C
               MAXSIZE = (MXCOR - I020)/IINTFP
C     
C Can we do it with the current algorithm available.
C     
               IF (MIN (NUMSYT, NUMSYG, DISSYT, DISSYG) .NE. 0) THEN
                  IF (MAXSIZE .GT. DISSYG) THEN
C
                     CALL MKDBLAA2(ICORE(I020), ICORE(I000), 
     &                             ICORE(I010), MAXSIZE, DISSYG, 
     &                             NUMSYG, DISSYT, NUMSYT, LISTG, 
     &                             LISTT, IRREPTAB, IRREPTMN,
     &                             IRREPGMN, IRREPGIJ, IRREPQAB, 
     &                             IRREPQIJ, IRREPX, POP(1,ISPIN), 
     &                             VRT(1,ISPIN), LISTZ, ISPIN)
                  ELSE
                     CALL INSMEM('GMNIJIND', I020, MXCOR)
                  ENDIF
C
               ENDIF
 100        CONTINUE
         ELSE
C     
C AB spin case. Again first try for in-core algorithm.
C     
            LISTT = IAPRT2AB5
            LISTG = INGMNAB
            LISTZ = 63
C
            DO 110 IRREPMN = 1, NIRREP
               IRREPTMN = IRREPMN
               IRREPTAB = DIRPRD(IRREPTMN, IRREPX)
               IRREPGMN = IRREPTMN
               IRREPGIJ = DIRPRD(IRREPGMN, IRREPX) 
               IRREPQIJ = IRREPGIJ
               IRREPQAB = IRREPQIJ
C     
C Get information about the size of the various arrays and allocate core memory.
C     
               DISSYG = IRPDPD(IRREPGMN, ISYTYP(1, LISTG))
               NUMSYG = IRPDPD(IRREPGIJ, ISYTYP(2, LISTG))
               DISSYT = IRPDPD(IRREPTAB, ISYTYP(1, LISTT))
               NUMSYT = IRPDPD(IRREPTMN, ISYTYP(2, LISTT))
C
               I000 = 1
               I010 = I000 + MAX(NUMSYT*DISSYT, DISSYT*NUMSYG)*IINTFP
               I020 = I010 + MAX(NUMSYT*DISSYT, DISSYT*NUMSYG)*IINTFP
               MAXSIZE = (MXCOR - I020)/IINTFP
C     
C Can we do it with the memory available
C     
               IF (MIN(NUMSYT, NUMSYG, DISSYT, DISSYG) .NE. 0) THEN
C
                  IF(MAXSIZE .GT. DISSYG) THEN
C     
                     CALL MKDBLAB2(ICORE(I020), ICORE(I000), 
     &                             ICORE(I010), MAXSIZE, DISSYG,
     &                             NUMSYG, DISSYT, NUMSYT, LISTG, 
     &                             LISTT, IRREPTAB, IRREPTMN,
     &                             IRREPGMN, IRREPGIJ, IRREPQAB, 
     &                             IRREPQIJ, IRREPX, POP(1,1),
     &                             POP(1,2), VRT(1,1), VRT(1,2),
     &                             LISTZ, IUHF)
C     
C No algorithm available
C
                  ELSE
                     CALL INSMEM('GMNIJIND', I020, MXCOR)
                  ENDIF
               ENDIF
C     
 110        CONTINUE
         ENDIF
 1000 CONTINUE
C
      RETURN
      END
