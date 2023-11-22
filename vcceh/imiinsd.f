C
      SUBROUTINE IMIINSD(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF)
C
C This routine calculates the term 
C
C  Z(IJ,AB) = - P_(IJ) SUM M T(IM,AB) I(MJ)
C
C to the doubles equation using symmetry packed arrays
C Spin-integrated formulas for this contribution
C
C In RHF 
C
C  - SUM m T(Im,Ab) I(mj) + SUM M T(Mj,Ab) I(MJ)
C
C IN UHF
C
C  - SUM M T(IM,AB) I(MJ) + SUM M T(MJ,AB) I(MI) [AAAA]
C
C  - SUM m T(Im,Ab) I(mj) + SUM M T(Mj,Ab) I(MI) [ABAB]
C
C  - SUM m T(im,ab) I(mj) + SUM m T(mj,ab) I(mi) [BBBB]
C
C In addition to those this subroutine also calculate
C
C  - SUM M T(M,A) I(MI) [AA]
C
C  - SUM m T(m,a) I(mi) [BB]
C
C contribution to the singles equation.
C
C Initially the term, 1/2 SUM E Hbar(E,M) T(E,I) is added to the 
C to the I(MI) intermediate before contraction with T2 
C amplitudes is carried out (here alpha and beta is implicit).
C
C Coded bu JG June/90 and July/90 and modified for quadratic term
C by Ajith 07/94
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER DIRPRD,DISSYT,DISSYZ,POP,VRT
      DIMENSION ICORE(MAXCOR), ITOFF1(8), ITOFF2(8)
      CHARACTER*6 SPCASE(2)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWO  
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,
     &             NF1AA,NF1BB,NF2AA,NF2BB
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /FLAGS/IFLAGS(100)
      COMMON /FILES/ LUOUT, MOINTS
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

      DATA ONE,ONEM,HALF/1.0D0,-1.0D0,0.5D0/
      DATA SPCASE /'AA =  ', 'BB =  '/
C
C Calculate size of I(M,I) array and load these arrays in to the 
C memory.
C
      NFOOAA = IRPDPD(IRREPX, 21)
      NFOOBB = IRPDPD(IRREPX, 22)
C
      I0AA  = MAXCOR + 1 - NFOOAA*IINTFP 
      MXCOR = MAXCOR - NFOOAA*IINTFP
C
      IF( IUHF .EQ. 0) THEN
         I0BB = I0AA
      ELSE
         I0BB  = I0AA -  NFOOBB*IINTFP
         MXCOR = MXCOR - NFOOBB*IINTFP
      ENDIF
C
      CALL GETLST(ICORE(I0AA), 1, 1, 1, 1, INTIMI)
C
         IF (IFLAGS(1) .GE. 50) THEN
            CALL HEADER('Checksum of I(MI) @-IMINSD', 0, LUOUT)
            
            WRITE(LUOUT, *) SPCASE(1), SDOT(NFOOBB, ICORE(I0AA),
     &                      1, ICORE(I0AA), 1)
         ENDIF
C
      IF (IUHF .EQ. 1) THEN
         CALL GETLST(ICORE(I0BB), 1, 1, 1, 2, INTIMI)
C
         IF (IFLAGS(1) .GE. 50) THEN
            CALL HEADER('Checksum of I(MI) @-IMINSD', 0, LUOUT)
            
            WRITE(LUOUT, *) SPCASE(2), SDOT(NFOOBB, ICORE(I0BB),
     &                      1, ICORE(I0BB), 1)
         ENDIF
C
      ENDIF
C
      NFVOAAX = IRPDPD(IRREPX, 9)
      NFVOBBX = IRPDPD(IRREPX, 10)
      NFVOAA1 = NTAA
      NFVOBB1 = NTBB
C
      I0TA  = I0BB - NFVOAAX*IINTFP
      I0ZA  = I0TA - NFVOAA1*IINTFP
C
      MXCOR = MXCOR - (NFVOAA1 + NFVOAAX)*IINTFP
C
C Get the T1(BET) amplitudes and Q(AI) singles already in the disk.
C
      CALL GETLST(ICORE(I0TA), 1, 1, 1, 1, IBPRT1AA)
      CALL GETLST(ICORE(I0ZA), 1, 1, 1, 3, 90)
C
      IF(IUHF.EQ.0) THEN
         I0TB = I0TA
         I0ZB = I0ZA
      ELSE
         I0TB = I0ZA - NFVOBBX*IINTFP
         I0ZB = I0TB - NFVOBB1*IINTFP
C
         MXCOR = MXCOR - (NFVOBBX + NFVOBB1)*IINTFP
C
         CALL GETLST(ICORE(I0TB), 1, 1, 1, 2, IBPRT1AA)
         CALL GETLST(ICORE(I0ZB), 1, 1, 1, 4, 90)
      ENDIF
C
C Now we are in a position to evaluate I(M,I) contribution to the
C singles and doubles equation. First do the singles. 
C
      DO 300 ISPIN = 1, IUHF + 1
C     
         IF (ISPIN .EQ. 1) THEN
            ITOFF1(1) = I0TA
            DO 5000 IRREPM = 2, NIRREP
               IRREP = IRREPM - 1
               IRREPAO = DIRPRD(IRREP, IRREPX)
               ITOFF1(IRREPM) = ITOFF1(IRREPM - 1) + VRT(IRREPAO, 1)*
     &                           POP(IRREP, 1)*IINTFP
 5000       CONTINUE
            IOFFF = I0AA
            IOFFZ = I0ZA
            I0Z   = I0ZA
         ELSE  
            ITOFF1(1) = I0TB
            DO 5100 IRREPM = 2, NIRREP
               IRREP = IRREPM - 1
               IRREPAO = DIRPRD(IRREP, IRREPX)
               ITOFF1(IRREPM) = ITOFF1(IRREPM - 1) + VRT(IRREPAO, 2)*
     &                           POP(IRREP, 2)*IINTFP
 5100       CONTINUE
            IOFFF = I0BB
            IOFFZ = I0ZB
            I0Z   = I0ZB
         ENDIF
C
         DO 250 IRREPI = 1, NIRREP
C            
            IRREPII = IRREPI
            IRREPIM = DIRPRD(IRREPII, IRREPX)
            IRREPTM = IRREPIM
            IRREPTA = DIRPRD(IRREPTM, IRREPX)
            IRREPQA = IRREPTA
            IRREPQI = IRREPQA
C            
            NOCCTM = POP(IRREPTM, ISPIN)
            NVRTTA = VRT(IRREPTA, ISPIN)
            NOCCIM = POP(IRREPIM, ISPIN)
            NOCCII = POP(IRREPII, ISPIN)
C
            IOFFT = ITOFF1(IRREPTM)
C
            IF (MIN(NVRTTA, NOCCTM, NOCCIM, NOCCII) .GT. 0) THEN 
C
               CALL XGEMM('N', 'N', NVRTTA, NOCCII, NOCCIM, ONEM, 
     &                     ICORE(IOFFT), NVRTTA, ICORE(IOFFF), NOCCIM, 
     &                     ONE, ICORE(IOFFZ), NVRTTA)
            ENDIF
C
            IOFFF = IOFFF + NOCCIM*NOCCII*IINTFP
            IOFFZ = IOFFZ + NOCCII*NVRTTA*IINTFP
C
 250     CONTINUE
C
C Update the  Q(AI) and Q(ai) singles lists
C
         CALL PUTLST(ICORE(I0Z), 1, 1, 1, ISPIN + 2, 90)
C
         IF (IFLAGS(1) .GE. 20) THEN
C
            IF (ISPIN .EQ. 1) THEN
               NSIZE = NFVOAA1
            ELSE
               NSIZE = NFVOBB1
            ENDIF
C
            CALL HEADER('Checksum @-IMIINSD', 0, LUOUT)
            
            WRITE(LUOUT, *) SPCASE(ISPIN), SDOT(NSIZE, ICORE(I0Z),
     &                      1, ICORE(I0Z), 1)
         ENDIF
C
 300  CONTINUE
C
C Before we calculate the I(M,J) contributions to the doubles equation we 
C need to add  1/2 SUM M Hbar(E,M)*T(E,J) to the I(M,J)) intermediate.
C Note the difference from the CC equations. In CC equations F(E,M)
C intermediate is used instead of Hbar(E,M) matrix elements.
C
      IOFFMJ = I0AA
      IOFT   = I0TA
C     
C Overwrite the T1(beta) by T1(alpha) and get Hbar(E,M)) and
C overwrite the allocaton for the Q(A,I).
C
      CALL GETLST(ICORE(I0TA), 1, 1, 1, 1, IAPRT1AA)
      CALL GETLST(ICORE(I0ZA), 1, 1, 1, 1, 93)
C
      ITOFF2(1) = I0ZA
C
      DO 5200 IRREPM = 2, NIRREP
         IRREP = IRREPM - 1 
         IRREPE = DIRPRD(IRREP, 1)
         ITOFF2(IRREPM) = ITOFF2(IRREPM - 1) + VRT(IRREPE, 1)*
     &                    POP(IRREP, 1)*IINTFP
 5200 CONTINUE
C      
      DO 20 IRREPJ = 1, NIRREP
C
         IRREPTJ = IRREPJ
         IRREPTE = DIRPRD(IRREPTJ, IRREPX)
         IRREPWE = IRREPTE
         IRREPWM = IRREPWE
C     
         NOCCM = POP(IRREPWM, 1)
         NVRTE = VRT(IRREPWE, 1)
         NOCCJ = POP(IRREPTJ, 1)
C
         IOFFEM = ITOFF2(IRREPWM)
C
         CALL XGEMM('T', 'N', NOCCM, NOCCJ, NVRTE, HALF,
     &               ICORE(IOFFEM), NVRTE, ICORE(IOFT), NVRTE, ONE,
     &               ICORE(IOFFMJ), NOCCM)
C
         IOFFMJ = IOFFMJ + NOCCM*NOCCJ*IINTFP
         IOFT   = IOFT   + NOCCJ*NVRTE*IINTFP
C
 20   CONTINUE
C     
      IF (IUHF .EQ. 1) THEN
C     
         IOFFMJ = I0BB
         IOFT   = I0TB
C
         CALL GETLST(ICORE(I0TB), 1, 1, 1, 2, IAPRT1AA)
         CALL GETLST(ICORE(I0ZB), 1, 1, 1, 2, 93)
C
         ITOFF2(1) = I0ZB
C     
         DO 5300 IRREPM = 2, NIRREP
            IRREP = IRREPM - 1 
            IRREPE = DIRPRD(IRREP, 1)
            ITOFF2(IRREPM) = ITOFF2(IRREPM - 1) + VRT(IRREPE, 2)*
     &                       POP(IRREP, 2)*IINTFP
 5300    CONTINUE
C
         DO 21 IRREPJ = 1, NIRREP
C
            IRREPTJ = IRREPJ
            IRREPTE = DIRPRD(IRREPTJ, IRREPX)
            IRREPWE = IRREPTE
            IRREPWM = IRREPWE
C     
            NOCCM = POP(IRREPWM, 2)
            NVRTE = VRT(IRREPWE, 2)
            NOCCJ = POP(IRREPTJ, 2)
C
            IOFFEM = ITOFF2(IRREPWM)
C
            CALL XGEMM('T', 'N', NOCCM, NOCCJ, NVRTE, HALF, 
     &                  ICORE(IOFFEM), NVRTE, ICORE(IOFT), NVRTE, ONE,
     &                  ICORE(IOFFMJ), NOCCM)
C
            IOFFMJ = IOFFMJ + NOCCM*NOCCJ*IINTFP
            IOFT   = IOFT   + NOCCJ*NVRTE*IINTFP
C
 21      CONTINUE
C
      ENDIF
C     
C Lets do the contribution to the doubles equation 
C First AAAA and BBBB spin cases.
C
      IF (IUHF .EQ. 1) THEN
C     
C These cases are only necessary in the UHF case and in the RHF
C AAAA amplitudes are calculated from ABAB.
C     
         DO 100 ISPIN = 1, 2
C     
            IF (ISPIN .EQ. 1) THEN
               I000  = I0AA
               NFSIZ = NFOOAA
            ELSE
               I000  = I0BB
               NFSIZ = NFOOBB
            ENDIF
C
            LISTT = ISPIN + (IAPRT2AA2 - 1)
            LISTZ = ISPIN + 60
C     
            DO 50 IRREPIM = 1, NIRREP 
C     
C Retrieve T2 amplitudes and calculate contribution to doubles
C     
               IRREPTIM = IRREPIM
               IRREPTAB = DIRPRD(IRREPTIM, IRREPX)
               IRREPQAB = IRREPTAB
               IRREPQIJ = IRREPQAB
C     
               NOCCSQ1 = 0
               NOCCSQ2 = 0
C
               DO 45 IRREPI = 1, NIRREP
                  IRREPJ = DIRPRD(IRREPI, IRREPQIJ)
                  IRREPM = DIRPRD(IRREPI, IRREPTIM)
                  NOCCSQ1 = NOCCSQ1 + POP(IRREPI, ISPIN)*
     &                      POP(IRREPJ, ISPIN)
                  NOCCSQ2 = NOCCSQ2 + POP(IRREPI, ISPIN)*
     &                      POP(IRREPM, ISPIN)
 45            CONTINUE
C
               DISSYT = IRPDPD(IRREPTAB, ISYTYP(1, LISTT)) 
               NUMSYT = IRPDPD(IRREPTIM, ISYTYP(2, LISTT))
               DISSYZ = IRPDPD(IRREPQAB, ISYTYP(1, LISTZ))
               NUMSYZ = IRPDPD(IRREPQIJ, ISYTYP(2, LISTZ))
C
               I001 = 1
               I002 = I001 + IINTFP*MAX(NOCCSQ1*DISSYT, DISSYT*NUMSYT,
     &                DISSYZ*NUMSYZ, DISSYT*NOCCSQ2, NOCCSQ1*DISSYZ)
               I003 = I002 + IINTFP*MAX(NOCCSQ1*DISSYZ, DISSYT*NUMSYT,
     &                DISSYZ*NUMSYZ, DISSYT*NOCCSQ2, NOCCSQ1*DISSYT)
C
               IF (MIN(NUMSYT, NUMSYZ, DISSYT, DISSYZ) .NE. 0) THEN
C
                  I004 = I003 + IINTFP*MAX(DISSYT, DISSYZ)
C
                  IF (I004 .LT. MXCOR) THEN
C     
C In core version
C     
                     CALL MKDBLAA4(ICORE(I001), ICORE(I002),
     &                             ICORE(I000), POP(1,ISPIN), 
     &                             NOCCSQ1, NOCCSQ2, DISSYT, DISSYZ,
     &                             NUMSYT, NUMSYZ, NFSIZ,  LISTT,
     &                             LISTZ, IRREPTIM, IRREPQIJ, IRREPX,
     &                             ISPIN, ICORE(I003))
                  ELSE
                     CALL INSMEM('IMIINSD', I004, MXCOR)
                  ENDIF
C     
               ENDIF
C
 50         CONTINUE
 100     CONTINUE
      ENDIF
C     
C AB spin case
C     
      LISTT = IAPRT2AB5
      LISTZ = 63
C     
C Loop over irreps
C     
      DO 200 IRREPIM = 1, NIRREP
C
         IRREPTIM = IRREPIM
         IRREPTAB = DIRPRD(IRREPTIM, IRREPX)
         IRREPQAB = IRREPTAB
         IRREPQIJ = IRREPQAB
C     
C Retrive T2(beta) amplitudes and calculate doubles contribution
C     
         DISSYT = IRPDPD(IRREPTAB, ISYTYP(1, LISTT))
         NUMSYT = IRPDPD(IRREPTIM, ISYTYP(2, LISTT))
         DISSYZ = IRPDPD(IRREPQAB, ISYTYP(1, LISTZ))
         NUMSYZ = IRPDPD(IRREPQIJ, ISYTYP(2, LISTZ))
C
         I001 = 1
         I002 = I001 + IINTFP*MAX(NUMSYT*DISSYT, NUMSYZ*DISSYZ)
         I003 = I002 + IINTFP*MAX(NUMSYZ*DISSYZ, NUMSYT*DISSYT)
C
         IF (MIN(NUMSYT, NUMSYZ, DISSYT, DISSYZ) .NE. 0) THEN
C
            I004 = I003 + IINTFP*MAX(DISSYT, DISSYZ, NUMSYT,
     &              NUMSYZ)*3
C
            IF (I004 .LT. MXCOR) THEN
C     
C In core algorithm
C     
               CALL MKDBLAB4(ICORE(I001), ICORE(I002), 
     &                       ICORE(I0AA),ICORE(I0BB),
     &                       POP(1,1), POP(1,2), VRT(1,1), VRT(1,2),
     &                       DISSYT, DISSYZ, NUMSYT, NUMSYZ, NFOOAA, 
     &                       NFOOBB, LISTT, LISTZ, IRREPTAB, IRREPTIM,
     &                       IRREPQAB, IRREPQIJ, IRREPX, IUHF,
     &                       ICORE(I003))
            ELSE
C     
C Out of core algorithm
C     
               CALL INSMEM('IMIINSD', I004, MXCOR)
            ENDIF
C
         ENDIF
 200  CONTINUE
C
      RETURN
      END
