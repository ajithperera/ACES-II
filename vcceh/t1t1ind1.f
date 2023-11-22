C
      SUBROUTINE T1T1IND1(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF)
C
C This routine calculates the following term required in the 
C quadratic term
C
C  - P_(IJ)P_(AB) SUM M E  T(E,I) T(A,M) Hbar(MB,EJ)
C
C In RHF, there are two terms
C
C  - SUM m E T(E,I) T(b,m) Hbar(mA,jE)
C
C  - SUM m e T(e,j) T(b,m) Hbar(Am,Ie)   
C
C In UHF, there are additional four terms
C
C  - SUM M E T(E,J) T(A,M) Hbar(MB,IE)
C
C  - SUM m e T(e,j) T(a,m) Hbar(mb,ie)  
C
C  - SUM M E T(E,I) T(A,M) <Eb//mJ>  
C
C  - SUM M e T(e,j) T(A,M) <eI//bM>  
C
C The multiplication is carried out in (The logic is the same for all 
C six cases)
C  
C The transpositions within RHF is carried out in the next call
C to the subroutine GMBEJIND, which evaluate the ring contribution 
C to the doubles equation. Also reordering of the amplitudes etc.
C is done there. Note also that the T2 contributions are multiplied 
C by minus one in the GMBEJIND code so that here the negative contribution
C will be calculated
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER DIRPRD,DISSYW,DISSYZ,DISTMP,VRT,POP
      CHARACTER*8 SPCASE(2)
      DIMENSION ICORE(MAXCOR),IOFFT1A(8,2), IOFFT1B(8,2)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWO
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,
     &             NF2AA,NF2BB
      COMMON/FILES/LUOUT, MOINTS
      COMMON/FLAGS/ IFLAGS(100)
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
      DATA AZERO,ONE,ONEM /0.D0,1.D0,-1.D0/
      DATA SPCASE /'AAAA', 'BBBB'/
C
C Allocate memory for T1(ALPHA) and T1(BETA) amplitudes and get them 
C from the disk. 

      LISTA  = IAPRT1AA
      LISTB  = IBPRT1AA
C     
      CALL GETPERT1(ICORE, MAXCOR, INTCOR, IUHF, IRREPX, LISTA, 
     &              IOFFSET, IOFFT1A)
C
      CALL GETPERT1(ICORE, INTCOR, MXCOR, IUHF, IRREPX, LISTB, 
     &              IOFFSET, IOFFT1B)
C
C AAAA and BBBB spin cases (Only neccesary for UHF)
C 
C  - SUM M E T(E,I) T(B,M) Hbar(MA,JE)   
C
C  - SUM m e T(e,i) T(b,m) Hbar(ma,je)     
C
C Calculated quantity is the negative of the above two equations
C so the total sign is positive.
C
      IF (IUHF .EQ. 1) THEN
C     
         DO 100 ISPIN = 1, 2
C     
C  ITRANS = 1 Transpose of result is necessary 
C  INCREM = 0   Initialize target list, do not update
C
            ITRANS = 0
            INCREM = 0
C 
C Prefactor and integral lists of the corresponding contribution
C
            FACT  = ONE
            LISTW = 53 + ISPIN
            LISTZ = 39 + ISPIN
C
C Loop over irrreps
C
            DO 50 IRREP = 1, NIRREP
C
               IRREPWAJ = IRREP
               IRREPWEM = IRREPWAJ
               IRREPQAJ = IRREPWAJ
               IRREPQBI = IRREPQAJ
               IRREPIMI = DIRPRD(IRREPQAJ, IRREPX)
C
               NOCCSQ = 0
               DO 45 IRREPM = 1, NIRREP
                  IRREPI = DIRPRD(IRREPM, IRREPIMI)
                  NOCCSQ = NOCCSQ + POP(IRREPM, ISPIN)*
     &                     POP(IRREPI, ISPIN)
 45            CONTINUE
C     
               DISSYW = IRPDPD(IRREPWEM, ISYTYP(1, LISTW))
               NUMSYW = IRPDPD(IRREPWAJ, ISYTYP(2, LISTW))
               DISSYZ = IRPDPD(IRREPQAJ, ISYTYP(1, LISTZ))
               NUMSYZ = IRPDPD(IRREPQBI, ISYTYP(2, LISTZ))
               DISTMP = DISSYZ
               NUMTMP = NUMSYZ
C
               MSIZE = MAX(NUMSYW*DISSYW, NUMSYZ*DISSYZ, DISSYZ*NOCCSQ)
C     &                     DISSYW*NOCCSQ)
C
               I001 = 1
               I002 = I001 + IINTFP*MSIZE
               I003 = I002 + IINTFP*MSIZE
               IF (I003 .LT. MXCOR) THEN
C     
C In core multiplication 
C     
                  CALL MKT1DBL(ICORE(I001), ICORE(I001), ICORE(I002),
     &                         ICORE(IOFFT1A(IRREPX, ISPIN)), 
     &                         ICORE(IOFFT1B(IRREPX, ISPIN)),
     &                         DISSYW, DISSYZ, DISTMP, NUMSYW, NUMSYZ,
     &                         NUMTMP, NOCCSQ, POP(1,ISPIN), 
     &                         POP(1,ISPIN), VRT(1,ISPIN), VRT(1,ISPIN),
     &                         LISTW, LISTZ, FACT, ITRANS, INCREM, 
     &                         IRREPWEM, IRREPWAJ, IRREPQAJ, IRREPQBI,
     &                         IRREPX)
C
               ELSE
C     
C Sorry, no out of core algoritm available. 
C     
                  CALL INSMEM('T1T1IND', MXCOR, I003)
               ENDIF
C
 50         CONTINUE
 100     CONTINUE
C
      ENDIF
C     
C ABAB Spin case 
C 
C  - SUM M E T(E,I) T(b,m) Hbar(mA,jE)
C
C  - SUM M e T(e,j) T(A,M) Hbar(bM,eI)
C
      DO 200 ISPIN = 1, IUHF + 1
C     
         FACT = ONE
C     
C On the first pass, initialize output list and do not transpose 
C result, on the second pass, update output list and transpose 
C     
         IF (ISPIN .EQ. 1) THEN
            I0T1A  = IOFFT1A(IRREPX, 1)
            I0T1B  = IOFFT1B(IRREPX, 2)
            IAOFF  = 1
            IBOFF  = 2
            ITRANS = 0
            INCREM = 0
         ELSE
            I0T1A  = IOFFT1A(IRREPX, 2)
            I0T1B  = IOFFT1B(IRREPX, 1)
            IAOFF  = 2
            IBOFF  = 1
            ITRANS = 1
            INCREM = 1
         ENDIF
C     
C Integral and target list
C     
         LISTW = 57 + ISPIN
         LISTZ = 41 + 2*IUHF
C     
         DO 150 IRREP = 1, NIRREP
C
            IRREPWAJ = IRREP
            IRREPWEM = IRREPWAJ
            IRREPQAJ = IRREPWAJ
            IRREPQBI = IRREPQAJ
            IRREPIMI = DIRPRD(IRREPQAJ, IRREPX)
C     
            DISSYW = IRPDPD(IRREPWEM, ISYTYP(1, LISTW))
            NUMSYW = IRPDPD(IRREPWAJ, ISYTYP(2, LISTW))
            DISSYZ = IRPDPD(IRREPQAJ, ISYTYP(1, LISTZ))
            NUMSYZ = IRPDPD(IRREPQBI, ISYTYP(2, LISTZ))
C
            IF(ITRANS .EQ. 0) THEN
               DISTMP = DISSYZ
               NUMTMP = NUMSYZ
            ELSE
               DISTMP = NUMSYZ
               NUMTMP = DISSYZ
            ENDIF
C
            NOCCSQ = IRPDPD(IRREPIMI, 14)
C
            MSIZE = MAX(NUMSYW*DISSYW, NUMSYZ*DISSYZ, DISTMP*NOCCSQ,
     &              DISSYZ*NOCCSQ, DISSYW*NOCCSQ)
C
            I001 = 1
            I002 = I001 + IINTFP*MSIZE
            I003 = I002 + IINTFP*MSIZE 
C
            IF(I003 .LT. MXCOR) THEN
C     
C In core algorithm
C     
               CALL MKT1DBL(ICORE(I001), ICORE(I001), ICORE(I002),
     &                      ICORE(I0T1A), ICORE(I0T1B), DISSYW, DISSYZ,
     &                      DISTMP, NUMSYW, NUMSYZ, NUMTMP, NOCCSQ, 
     &                      POP(1,3-ISPIN), POP(1,ISPIN),
     &                      VRT(1,3-ISPIN), VRT(1,ISPIN), LISTW, LISTZ,
     &                      FACT, ITRANS, INCREM, IRREPWEM, IRREPWAJ,
     &                      IRREPQAJ, IRREPQBI, IRREPX)
C
            ELSE
C     
C Soory, no out of core algorithm available
C     
               CALL INSMEM('T1T1IND', MXCOR, I003)               
C
            ENDIF
C
 150     CONTINUE
 200  CONTINUE
C     
C  ABAB Spin case second contribution 
C     
C  - SUM m e T(e,j) T(b,m) Hbar(mA,eI)
C     
C  - SUM M E T(E,I) T(A,M) Hbar(Mb,Ej)  (UHF only)
C     
      DO 300 ISPIN = 1, IUHF + 1
C     
         FACT = ONE
C     
C On the first pass, initialize target list and do not transpose 
C calculated array. In the second pass, update the target list
C and transpose the result array.
C     
         IF (ISPIN .EQ. 1) THEN
            I0T1A  = IOFFT1A(IRREPX, 2)
            I0T1B  = IOFFT1B(IRREPX, 2)           
            ITRANS = 0
            INCREM = 0
         ELSE
            I0T1A  = IOFFT1A(IRREPX, 1)
            I0T1B  = IOFFT1B(IRREPX, 1)           
            ITRANS = 1
            INCREM = 1
         ENDIF
C     
C Integral and target list
C     
         LISTW = 58 - ISPIN
         IF (IUHF .EQ. 0) LISTW = 56
         LISTZ = 42
C     
         DO 250 IRREP = 1, NIRREP
C
            IRREPWBJ = IRREP
            IRREPWEM = IRREPWBJ
            IRREPWAI = IRREPWBJ
            IRREPQAI = IRREPWAI
            IRREPQBJ = IRREPQAI
            IRREPIMI = DIRPRD(IRREPQAI, IRREPX)
C     
            NOCCSQ = 0
            DO 245 IRREPJ = 1, NIRREP
               IRREPM =DIRPRD(IRREPJ, IRREPIMI)
               NOCCSQ = NOCCSQ + POP(IRREPJ, 3-ISPIN)*
     &                  POP(IRREPM, 3-ISPIN)
 245        CONTINUE
C     
            DISSYW = IRPDPD(IRREPWEM, ISYTYP(1, LISTW))
            NUMSYW = IRPDPD(IRREPWAI, ISYTYP(2, LISTW))
            DISSYZ = IRPDPD(IRREPQAI, ISYTYP(1, LISTZ))
            NUMSYZ = IRPDPD(IRREPQBJ, ISYTYP(2, LISTZ))
C
            IF (ITRANS .EQ. 0) THEN
               DISTMP = DISSYZ
               NUMTMP = NUMSYZ
            ELSE 
               DISTMP = NUMSYZ
               NUMTMP = DISSYZ
            ENDIF
C
            MSIZE = MAX(NUMSYW*DISSYW, NUMSYZ*DISSYZ, DISTMP*NOCCSQ,
     &              DISSYZ*NOCCSQ, NUMTMP*NOCCSQ)
            I001 = 1
            I002 = I001 + IINTFP*MSIZE
            I003 = I002 + IINTFP*MSIZE
C
            IF(I003 .LT. MXCOR) THEN
C     
C In core algorithm
C     
               CALL MKT1DBL(ICORE(I001), ICORE(I001), ICORE(I002),
     &                      ICORE(I0T1A), ICORE(I0T1B), DISSYW, DISSYZ,
     &                      DISTMP, NUMSYW, NUMSYZ, NUMTMP, NOCCSQ, 
     &                      POP(1,3-ISPIN), POP(1,3-ISPIN), 
     &                      VRT(1,3-ISPIN), VRT(1,3-ISPIN), LISTW, 
     &                      LISTZ, FACT, ITRANS, INCREM, IRREPWEM,
     &                      IRREPWAI, IRREPQAI, IRREPQBJ, IRREPX)
C
            ELSE 
C     
C Sorry, no out of core algorithm
C     
               CALL INSMEM('T1T1IND', MXCOR, I003)               
            ENDIF
C
 250     CONTINUE
 300  CONTINUE
C     
C All done, return
C     
      RETURN
      END
