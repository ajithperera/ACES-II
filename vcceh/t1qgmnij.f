C
      SUBROUTINE T1QGMNIJ(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF)
C
C This routine calculates the following contribution to the W(MN,IJ)
C intermediates
C     Z(MN,IJ) = P_(IJ) SUM E T(E,J)*Hbar(MN,IE)
C
C Spin integrated expressions for RHF and UHF are given below.
C
C  RHF:
C      SUM e Hbar(Mn,Ie)*T(e,j) + TRANSPOSITION OF I AND J
C
C  UHF:
C
C AAAA:  SUM E Hbar(MN,IE)*T(E,J) + ANTISYMMETRIZATION
C
C BBBB:  SUM e Hbar(mn,ie)*T(e,j) + ANTISYMMETRIZATION
C
C BABA:  SUM e Hbar(Mn,Ie)*T(e,j) + SUM E Hbar(Mn,EJ)*T(E,I)
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER DIRPRD,DISSYW,DISSYZ,DISSYWA,DISSYWB,POP,VRT
      DIMENSION ICORE(MAXCOR)
C     
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYM/POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,
     &            NF2AA,NF2BB
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTIN,IINTFP,IALONE,IBITWO
      COMMON /FILES/ LUOUT,MOINTS
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
C Allocate core memory for T1 list and I(MI).
C     
       MXCOR = MAXCOR
C
       LIST1  = IAPRT1AA
C     
       NFVOAA = IRPDPD(IRREPX, 9)
       NFVOBB = IRPDPD(IRREPX, 10)
C     
       I0TA  = MXCOR + 1 - NFVOAA*IINTFP
       MXCOR = MXCOR - NFVOAA*IINTFP
C       
       CALL GETLST(ICORE(I0TA), 1, 1, 1, 1, LIST1)
C
       IF (IUHF .EQ. 0) THEN
          I0TB = I0TA
       ELSE
          I0TB = I0TA - NFVOBB*IINTFP
          MXCOR = MXCOR - NFVOBB*IINTFP
          CALL GETLST(ICORE(I0TB), 1, 1, 1, 2, LIST1)
       ENDIF 
C     
C AAAA and BBBB spin cases (UHF only).
C     
       IF (IUHF .EQ. 1) THEN
C     
          DO 1000 ISPIN = 1, 2
C     
             IF (ISPIN .EQ. 1) THEN
                IT  = I0TA
                NFT = NFVOAA
             ELSE
                IT  = I0TB
                NFT = NFVOBB
             ENDIF
C
             LISTW  = 6 + ISPIN
             LSTGMN = (INGMNAA - 1) + ISPIN
C     
C Loop over all irreps.
C     
             DO 100 IRREPMN = 1, NIRREP
C
                IRREPWMN = IRREPMN
                IRREPWEI = IRREPWMN
                IRREPGMN = IRREPWMN
                IRREPGIJ = DIRPRD(IRREPGMN, IRREPX)
C     
                NOCCSQ = 0
                DO 90 IRREPJ = 1, NIRREP
                   IRREPI = DIRPRD(IRREPJ, IRREPGIJ)
                   NOCCSQ = NOCCSQ + POP(IRREPI, ISPIN)*
     &                      POP(IRREPJ, ISPIN)
 90             CONTINUE
C
                DISSYW = IRPDPD(IRREPWMN, ISYTYP(1, LISTW))
                NUMSYW = IRPDPD(IRREPWEI, ISYTYP(2, LISTW))
                DISSYZ = IRPDPD(IRREPGMN, ISYTYP(1, LSTGMN))
                NUMSYZ = IRPDPD(IRREPGIJ, ISYTYP(2, LSTGMN))
C
                I001 = 1
                I002 = I001 + IINTFP*MAX(DISSYW*NUMSYW, NUMSYZ*DISSYZ)
                I003 = I002 + IINTFP*DISSYW*NOCCSQ
                I004 = I003 + MAX(NUMSYW, DISSYW, NUMSYZ, DISSYZ)
C
                IF (MIN (NUMSYW, DISSYW) .NE. 0) THEN
C
                   IF(I004 .LT. MXCOR) THEN
C     
C In core version
C     
                      CALL MKT1GMNIJAA(ICORE(I001), ICORE(I002), 
     &                                 ICORE(IT), POP(1,ISPIN),
     &                                 VRT(1,ISPIN), DISSYZ, DISSYW,
     &                                 NUMSYZ, NUMSYW, NOCCSQ, NFT,
     &                                 LISTW, LSTGMN, IRREPWEI, 
     &                                 IRREPGIJ, IRREPX, ISPIN,
     &                                 ICORE(I003))
                   ELSE
C     
C Out of core algorithm
C     
                      CALL INSMEM('T1GMNIJAB', I004, MXCOR)
                   ENDIF
                ENDIF
C     
 100         CONTINUE
 1000     CONTINUE
C
       ENDIF
C     
C AB spin case
C     
       LISTWA = 10
       LISTWB = 9
       LSTGMN = INGMNAB
C     
C Loop over irreps
C     
       DO 200 IRREPMN = 1, NIRREP
C
          IRREPWAMN = IRREPMN
          IRREPWAEI = IRREPWAMN
          IRREPWBMN = IRREPWAMN
          IRREPWBEJ = IRREPWBMN
          IRREPGMN  = IRREPWAMN
          IRREPGIJ  = DIRPRD(IRREPGMN, IRREPX)
C     
          DISSYZ  = IRPDPD(IRREPGMN, ISYTYP(1, LSTGMN))
          NUMSYZ  = IRPDPD(IRREPGIJ, ISYTYP(2, LSTGMN))
          DISSYWA = IRPDPD(IRREPWAMN, ISYTYP(1, LISTWA))
          NUMSYWA = IRPDPD(IRREPWAEI, ISYTYP(2, LISTWA))
C
          IF (IUHF .EQ. 0) THEN
             DISSYWB = 0
             NUMSYWB = 0
          ELSE
             DISSYWB = IRPDPD(IRREPWBMN, ISYTYP(1, LISTWB))
             NUMSYWB = IRPDPD(IRREPWBEJ, ISYTYP(2, LISTWB))
          ENDIF
C
          I001 = 1
          I002 = I001 + IINTFP*NUMSYZ*DISSYZ
          I003 = I002 + IINTFP*MAX(NUMSYWA*DISSYWA, NUMSYWB*DISSYWB,
     &           NUMSYZ*DISSYZ)
          I004 = I003 + 3*IINTFP*MAX(NUMSYZ, NUMSYWA, NUMSYWB, 
     &           DISSYZ, DISSYWA, DISSYWB)
C
          IF (I004 .LT. MXCOR) THEN
C     
             CALL MKT1GMNIJAB(ICORE(I001), ICORE(I002), 
     &                        ICORE(I0TA), ICORE(I0TB),
     &                        POP(1,1), POP(1,2), VRT(1,1),
     &                        VRT(1,2), DISSYZ, DISSYWA, 
     &                        DISSYWB, NUMSYZ, NUMSYWA,
     &                        NUMSYWB, NFVOAA, NFVOBB, LISTWA, 
     &                        LISTWB, LSTGMN, IRREPWAEI,
     &                        IRREPWBEJ, IRREPGMN, IRREPGIJ,
     &                        IRREPX, IUHF, ICORE(I003))
          ELSE
C
             CALL INSMEM('T1GMNIJAB', I004, MXCOR)
C     
C Out of core algorithm
C
          ENDIF
C
 200   CONTINUE
C     
       RETURN
       END
