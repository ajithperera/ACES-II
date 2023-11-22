C
      SUBROUTINE GABEFIND(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF)
C
C This routine calculates the term 
C
C                M               IJ
C  - P_(AB) SUM T(ALP) [1/2(SUM T(BET)*Hbar(AM,EF)]
C                B          E,F  EF
C
C to the doubles amplitude using symmetry adapted integrals and
C amplitudes. By doing this we avoid the storage G(AB,EF) intermediate.
C
C IN RHF
C
C  - SUM m T(b,m) [ SUM E,f T(Ef,Ij)*Hbar(Am,Ef)] + Transpose of A and B
C
C In UHF
C
C  - SUM M T(B,M) [SUM E<F T(EF,IJ)*Hbar(AM,EF)] + Antisymetrization in A,B
C
C  - SUM m T(b,m) [SUM e<f T(ef,ij)*Hbar(am,ef)] + Antisymetrization in a,b
C
C  - SUM m T(b,m) [SUM E,f T(Ef,Ij)*Hbar(Am,Ef)] 
C
C  - SUM M T(A,M) [SUM E,f T(Ef,Ij)*Hbar(Mb,Ef)]
C
C There is in core and out core algorithm to deal with Hbar(AB,CI) 
C integrals. Here T1(alpha) and T2(beta) is implicit.
C
C Originally coded by JG JUNE/90 and modified for the quadratic
C term by ajith 06/94
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER DIRPRD, DISSYWA, DISSYWB, DISSYZ, POP, VRT, DISSYW,
     &        DISSYT
      DIMENSION ICORE(MAXCOR)
C
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWO
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON/SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,
     &            NF1AA,NF1BB,NF2AA,NF2BB
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
      MXCOR = MAXCOR
      LIST1 = IAPRT1AA
C
C Determine the size of the T1(ALP) vector, now T1(ALP) vector is not 
C totally symmetric.
C  
      NFVOAA = IRPDPD(IRREPX, 9)
      NFVOBB = IRPDPD(IRREPX, 10)
C
      I0TA  = MXCOR + 1 - NFVOAA*IINTFP
      MXCOR = MXCOR - NFVOAA*IINTFP
C
      IF (IUHF .EQ. 0) THEN
         I0TB = I0TA
      ELSE
         I0TB  = I0TA - NFVOBB*IINTFP
         MXCOR = MXCOR - NFVOBB*IINTFP
      ENDIF
C
C Get the T1(ALP) list to the bottom of the core
C
      CALL GETLST(ICORE(I0TA), 1, 1, 1, 1, LIST1)
C
      IF (IUHF .EQ. 1) THEN
         CALL GETLST(ICORE(I0TB), 1, 1, 1, 2, LIST1)
      ENDIF
C
C AAAA AND BBBB spin cases
C     
      IF (IUHF .EQ. 1) THEN
C     
         DO 100 ISPIN = 1, 2
C     
            IF(ISPIN .EQ. 1) THEN
               I000 = I0TA
                 NT = NFVOAA
            ELSE
               I000 = I0TB
                 NT = NFVOBB
            ENDIF
C
            NOCC = NOCCO(ISPIN)
            NVRT = NVRTO(ISPIN)
C
            LISTZ = ISPIN + 60
            LISTW = ISPIN + 26
            LISTT = ISPIN + (IAPRT2AA2 - 1)
C
C Note that now T2(alpha) amplitudes are overwritten by T2(beta)
C amplitudes. But this is not the case for T1 amplitudes, since
C both need to be kept in the disk. Notice that Q(IJ,AB) contribution
C overwrites the CC D(IJ, AB) contribution which is of course now
C possible, since the double amplitudes are now totally symmetric.
C
            DO 50 IRREPEF = 1, NIRREP 
C
               IRREPWEF = IRREPEF
               IRREPWAM = IRREPWEF
               IRREPTEF = IRREPWEF
               IRREPTIJ = DIRPRD(IRREPTEF, IRREPX)
               IRREPQIJ = IRREPTIJ
               IRREPQAB = IRREPQIJ
C
               NVRTSQ = 0
               DO 45 IRREPBO = 1, NIRREP
                  IRREPAO = DIRPRD(IRREPBO, IRREPQAB)
                  NVRTSQ = NVRTSQ + VRT(IRREPBO, ISPIN)*
     &                     VRT(IRREPAO, ISPIN)
 45            CONTINUE
C     
C Retrive T2 amplitudes and calculate new ones
C     
               DISSYZ = IRPDPD(IRREPQAB, ISYTYP(1, LISTZ))
               NUMSYZ = IRPDPD(IRREPQIJ, ISYTYP(2, LISTZ))
               DISSYT = IRPDPD(IRREPTEF, ISYTYP(1, LISTT))
               NUMSYT = IRPDPD(IRREPTIJ, ISYTYP(2, LISTT))
               DISSYW = IRPDPD(IRREPWEF, ISYTYP(1, LISTW))
               NUMSYW = IRPDPD(IRREPWAM, ISYTYP(2, LISTW))
C
               I001 = 1
               I002 = I001 + IINTFP*MAX(NUMSYZ*DISSYZ, NUMSYW*NUMSYT,
     &                NUMSYZ*NVRTSQ, NUMSYT*DISSYT)
               I003 = I002 + IINTFP*MAX(NUMSYT*DISSYT, NUMSYW*NUMSYT,
     &                NUMSYZ*NVRTSQ, NUMSYZ*DISSYZ)
               I004 = I003 + 3*IINTFP*MAX(NUMSYW,NUMSYZ,DISSYZ,DISSYW,
     &                DISSYT,NUMSYT)
C     
C For the integrals an in-core version is usually not possible
C try to allocate as much as possible
C     
               IF (MIN(NUMSYZ, DISSYZ) .NE. 0) THEN
                  MAXSIZE = (MXCOR - I004)/IINTFP
                  IF(MAXSIZE .GT. DISSYW) THEN
C
C In core version 
C     
                     CALL MKDBLAA1(ICORE(I001), ICORE(I001),
     &                             ICORE(I002), ICORE(I002),
     &                             ICORE(I004), MAXSIZE,
     &                             ICORE(I000), POP(1,ISPIN),
     &                             VRT(1,ISPIN), DISSYZ,
     &                             DISSYW, DISSYT, NUMSYZ, NUMSYW,
     &                             NUMSYT, NVRTSQ, NT, LISTT, LISTZ, 
     &                             LISTW, IRREPTEF, IRREPTIJ,
     &                             IRREPWEF, IRREPWAM, IRREPQAB,
     &                             IRREPQIJ, IRREPX, IUHF, ISPIN,
     &                             ICORE(I003))
                  ELSE
                     CALL INSMEM('GABEFIND', DISSYW, MAXSIZE)
                  ENDIF
               ENDIF
 50         CONTINUE
 100     CONTINUE     
      ENDIF
C     
C AB spin case
C     
      LISTZ  = 63
      LISTWA = 30
      LISTWB = 29
      LISTT  = IAPRT2AB5
C     
C Loop over irreps
C     
      DO 200 IRREPEF = 1, NIRREP
         IRREPWAEF = IRREPEF
         IRREPWAAM = IRREPWAEF
         IRREPWBEF = IRREPWAEF
         IRREPWBBM = IRREPWBEF
         IRREPTEF  = IRREPWAEF
         IRREPTIJ  = DIRPRD(IRREPTEF, IRREPX)
         IRREPQIJ  = IRREPTIJ
         IRREPQAB  = IRREPQIJ
C     
C Retrieve T2 amplitudes and calculate Z amplitudes
C     
         DISSYZ  = IRPDPD(IRREPQAB, ISYTYP(1, LISTZ))
         NUMSYZ  = IRPDPD(IRREPQIJ, ISYTYP(2, LISTZ))
         DISSYT  = IRPDPD(IRREPTEF, ISYTYP(1, LISTT))
         NUMSYT  = IRPDPD(IRREPTIJ, ISYTYP(2, LISTT))
         DISSYWA = IRPDPD(IRREPWAEF, ISYTYP(1, LISTWA))
         NUMSYWA = IRPDPD(IRREPWAAM, ISYTYP(2, LISTWA))
C
         IF (IUHF .EQ. 0) THEN
            DISSYWB = 0
            NUMSYWB = 0
         ELSE 
            DISSYWB = IRPDPD(IRREPWBEF, ISYTYP(1, LISTWB))
            NUMSYWB = IRPDPD(IRREPWBBM, ISYTYP(2, LISTWB))
         ENDIF
C     
C Holds Z-vector and T2-amplitudes in-core for results in-core
C     
         I001 = 1
         I002 = I001 + IINTFP*MAX(NUMSYZ*DISSYZ, NUMSYWA*NUMSYT,
     &          NUMSYWB*NUMSYT, NUMSYT*DISSYT)
         I003 = I002 + IINTFP*MAX(NUMSYT*DISSYT, NUMSYWA*NUMSYT,
     &          NUMSYWB*NUMSYT, NUMSYZ*DISSYZ)
         I004 = I003 + 3*IINTFP*MAX(NUMSYWA, NUMSYWB, NUMSYZ, 
     &          DISSYZ, DISSYWA, DISSYWB, DISSYT, NUMSYT)
C     
C For the integrals an in-core version is usually not possible.
C Try to allocate as much as possible
C     
         IF (MIN(NUMSYZ, DISSYZ) .NE. 0) THEN
            MAXSIZE = (MXCOR - I004)/IINTFP
            IF(MAXSIZE .GT. MAX(DISSYWA, DISSYWB)) THEN
C     
C In core algorithm and partial out of core algorithm
C     
               CALL MKDBLAB1(ICORE(I001), ICORE(I001), ICORE(I002),
     &                       ICORE(I002), ICORE(I004), MAXSIZE,
     &                       ICORE(I0TA), ICORE(I0TB), POP(1,1),
     &                       POP(1,2), VRT(1,1), VRT(1,2), DISSYZ,
     &                       DISSYWA, DISSYWB, DISSYT, NUMSYZ, NUMSYWA,
     &                       NUMSYWB, NUMSYT, NFVOAA, NFVOBB, LISTT,
     &                       LISTZ, LISTWA, LISTWB, IRREPWAEF, 
     &                       IRREPWAAM, IRREPWBEF, IRREPWBBM, IRREPTEF,
     &                       IRREPTIJ, IRREPQAB, IRREPQIJ, IRREPX,
     &                       IUHF, ICORE(I003))
            ELSE
C     
C Out of core algorithm
C     
               CALL INSMEM('GABEFIND', MAX(DISSYWA, DISSYWB), MAXSIZE)
            ENDIF
C
         ENDIF
 200  CONTINUE
C
      RETURN
      END
