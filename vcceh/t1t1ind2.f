C
      SUBROUTINE T1T1IND2(ICORE, MAXCOR, IRREPX, IOFFSET, ITYPE, IUHF)
C
C This subroutine calculates 
C               A        B
C + P_(AB) SUM T(ALP) * T(BET) * Hbar(MN, IJ) and 
C          M<N  M        N
C
C               E        F
C + P_(AB) SUM T(ALP) * T(BET) * Hbar(AB, EF) to the
C          E<F  I        J
C doubles equation. In the CC code these two terms were included
C in Tau(MN,AB) in the doubles equation. For the  quadratic term
C it is not possible to include this term in the Tau
C intemediates. Spin integrated formulas are as follows. Here
C alpha and beta is implicit.
C
C RHF
C      
C  Z(Ij,Ab) = T(M,A)*T(n,b)*Hbar(Mn,Ij) [ABAB]
C
C UHF
C
C  Z(IJ,AB) = [T(M,A)*T(N,B) - T(M,B)*T(N,A)]*Hbar(MN,IJ) [AAAA]
C
C  Z(IJ,AB) = [T(m,a)*T(n,b) - T(m,b)*T(n,a)]*Hbar(mn,ij) [BBBB]
C
C  Z(Ij,Ab) = T(M,A)*T(n,b)*Hbar(Mn,Ij) [ABAB]
C 
C Similar expressions can be written for the second term.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER DIRPRD, POP, VRT, DISSYH, DISSYQ, DISTMP
C
      DIMENSION ICORE(MAXCOR), IOFFT1A(8,2), IOFFT1B(8,2)
C
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,
     &             NF2AA,NF2BB
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
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
C Allocate  memory for T1(ALPHA) and T1(BETA) amplitudes and get
C them from the disk.
C     
      LISTA = IAPRT1AA
      LISTB = IBPRT1AA
C     
      CALL GETPERT1(ICORE, MAXCOR, INTCOR, IUHF, IRREPX, LISTA,
     &              IOFFSET, IOFFT1A)
C     
      CALL GETPERT1(ICORE, INTCOR, MXCOR, IUHF, IRREPX, LISTB,
     &              IOFFSET, IOFFT1B)
C
C AAAA and BBBB spin cases (UHF only)
C     
      IF (IUHF .EQ. 1) THEN
C
         DO 10 ISPIN = 1, 2
C             
            LISTQ = 60 + ISPIN
C
            IF (ITYPE .EQ. 1) THEN
               LISTH = 50 + ISPIN
            ELSE
               LISTH = 230 + ISPIN
            ENDIF
C     
C Loop over irreps
C
            DO 20 IRREPMN = 1, NIRREP
C     
               IRREPHMN = IRREPMN
               IRREPHIJ = IRREPHMN
               IRREPTMN = IRREPHMN
               IRREPTAB = IRREPTMN
               IRREPQAB = IRREPTAB
               IRREPQIJ = IRREPQAB
C
               DISSYH = IRPDPD(IRREPHMN, ISYTYP(1, LISTH))
               NUMSYH = IRPDPD(IRREPHIJ, ISYTYP(2, LISTH))
               DISSYQ = IRPDPD(IRREPQAB, ISYTYP(1, LISTQ))
               NUMSYQ = IRPDPD(IRREPQIJ, ISYTYP(2, LISTQ))
C 
               IF (ITYPE .EQ. 1) THEN
                  DISTMP = DISSYQ
                  NUMTMP = DISSYH
               ELSE IF (ITYPE .EQ. 6) THEN
                  DISTMP = NUMSYH
                  NUMTMP = NUMSYQ
               ENDIF
C
               I000 = 1
               I010 = I000 + MAX(DISTMP*NUMTMP, DISSYQ*NUMSYQ)*IINTFP
               I020 = I010 + DISSYQ*NUMSYQ*IINTFP
C
               MAXSIZE = (MXCOR - I020)/IINTFP
C
               IF (MIN(DISSYQ, NUMSYQ, DISSYH, NUMSYH) .NE. 0) THEN
C
                  IF (MAXSIZE .GT. DISSYH) THEN
C     
                     CALL MKDBLAA5(ICORE(I000), ICORE(I010),
     &                             ICORE(I020), 
     &                             ICORE(IOFFT1A(IRREPX, ISPIN)),
     &                             ICORE(IOFFT1B(IRREPX, ISPIN)),
     &                             POP(1,ISPIN), POP(1, ISPIN),
     &                             VRT(1, ISPIN), VRT(1, ISPIN),
     &                             DISSYH, NUMSYH, DISSYQ, NUMSYQ,
     &                             DISTMP, NUMTMP, LISTH, LISTQ, 
     &                             IRREPHIJ, IRREPTMN, IRREPQIJ, 
     &                             IRREPX, MAXSIZE, ITYPE, ISPIN, IUHF)
                  ELSE
                     CALL INSMEM('T1T1IND2', I020, MXCOR)
                  ENDIF
               ENDIF
 20         CONTINUE
 10      CONTINUE
C
      ENDIF
C
C ABAB spin case. Have contributions in both UHF and RHF cases.
C 
      LISTQ = 63 
      ISPIN = 3
C
      IF (ITYPE .EQ. 1) THEN
         LISTH = 53
      ELSE
         LISTH = 233
      ENDIF
C
      DO 30 IRREPMN = 1, NIRREP
C     
         IRREPHMN = IRREPMN
         IRREPHIJ = IRREPHMN
         IRREPTMN = IRREPHMN
         IRREPTAB = IRREPTMN
         IRREPQAB = IRREPTAB
         IRREPQIJ = IRREPQAB
C
         DISSYH = IRPDPD(IRREPHMN, ISYTYP(1, LISTH))
         NUMSYH = IRPDPD(IRREPHIJ, ISYTYP(2, LISTH))
         DISSYQ = IRPDPD(IRREPQAB, ISYTYP(1, LISTQ))
         NUMSYQ = IRPDPD(IRREPQIJ, ISYTYP(2, LISTQ))
C
         IF (ITYPE .EQ. 1) THEN
            DISTMP = DISSYQ
            NUMTMP = DISSYH
         ELSE IF (ITYPE .EQ. 6) THEN
            DISTMP = NUMSYH
            NUMTMP = NUMSYQ
         ENDIF
C
         I000 = 1
         I010 = I000 + MAX(DISTMP*NUMTMP, DISSYQ*NUMSYQ)*IINTFP
         I020 = I010 + DISSYQ*NUMSYQ*IINTFP
C
         MAXSIZE = (MXCOR - I020)/IINTFP
         
         IF (MIN(DISSYQ, NUMSYQ, DISSYH, NUMSYH) .NE. 0) THEN
C
            IF (MAXSIZE .GT. DISSYH) THEN
C       
               CALL MKDBLAB5(ICORE(I000), ICORE(I010), ICORE(I020),
     &                       ICORE(IOFFT1A(IRREPX, 1)),
     &                       ICORE(IOFFT1B(IRREPX, 2)), POP(1, 1),
     &                       POP(1, 2), VRT(1, 1), VRT(1, 2), DISSYH,
     &                       NUMSYH, DISSYQ, NUMSYQ, DISTMP, NUMTMP, 
     &                       LISTH, LISTQ, IRREPHIJ, IRREPTMN, 
     &                       IRREPQIJ, IRREPX, MAXSIZE, ITYPE, ISPIN,
     &                       IUHF)
            ELSE
               CALL INSMEM('T1T1IND2', I020, MXCOR)
            ENDIF
         ENDIF
C     
 30   CONTINUE
C      
      RETURN
      END
