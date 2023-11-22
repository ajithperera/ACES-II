C 
      SUBROUTINE MKT2IMIAA(T, W, Z, ISPIN, POP, VRT, DISSYW, DISSYT,
     &                     NUMSYW, NUMSYT, NOCCSQ, NOCC2SQ1, NOCC2SQ2, 
     &                     IRREPWL, IRREPWR, IRREPTL, IRREPTR, IRREPX,
     &                     TMP, IOFFSET, IUHF, FACT)
C
C Calculate the following contributions to the F(MI) intemediate
C depending on the value of ISPIN. 
C
C  Z(MI) = SUM N SUM E<F {T(IN,EF)*Hbar(MN,EF)} [ISPIN = 1]
C
C  Z(mi) = SUM n SUM e<f {T(in,ef)*Hbar(mn,ef)} [ISPIN = 2]
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER DISSYT,DISSYW,DIRPRD,POP,VRT
      DIMENSION W(DISSYW,NOCC2SQ1),T(DISSYT,NOCC2SQ2),Z(NOCCSQ)
      DIMENSION TMP(1), IWOFF(8) 
      DIMENSION POP(8),VRT(8)
C
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
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
      DATA ONEM,ONE,HALF,TWO /-1.0D0, 1.D0, 0.5D0, 2.D0/
C     
      IND(I,J) = ((I-2)*(I-1))/2+J
C
C Pick up first the relevent T2 and W pieces
C     
      LISTT = ISPIN + (IAPRT2AA2 - 1)       
      LISTW = ISPIN + 13
C
      CALL GETLST(W, 1, NUMSYW, 2, IRREPWR, LISTW)
      CALL GETLST(T, 1, NUMSYT, 1, IRREPTR, LISTT)
C     
      CALL SYMEXP(IRREPWR, POP, DISSYW, W)
      CALL SYMEXP(IRREPTR, POP, DISSYT, T)
C     
      JOFF = 1
      KOFF = 1
      IWOFF(1) = 1 
C
      DO 5000 IRREPM = 2, NIRREP
C
         IRREP = IRREPM - 1
         IRREPN = DIRPRD(IRREP, IRREPWR)
C
         IWOFF(IRREPM) = IWOFF(IRREPM - 1) + POP(IRREP)*POP(IRREPN)
C
 5000 CONTINUE
C     
      DO 90 IRREPI = 1, NIRREP
C
         IRREPN = DIRPRD(IRREPI, IRREPTR)
         IRREPM = DIRPRD(IRREPN, IRREPWR)
C     
C Determine irrepn and irrepi.
C
         NOCCM = POP(IRREPM)
         NOCCI = POP(IRREPI)
         NOCCN = POP(IRREPN)
C
         IOFF = IWOFF(IRREPM)
C     
C If zero, nothing to compute
C
         IF (NOCCM .GT. 0 .AND. NOCCN .GT. 0) THEN

            CALL XGEMM('T', 'N', NOCCM, NOCCI, DISSYW*NOCCN, FACT,
     &                  W(1,IOFF), NOCCN*DISSYW, T(1,JOFF),
     &                  NOCCN*DISSYT, ONE, Z(KOFF), NOCCM)
         ENDIF
C
C Update the intermediate for Z
C
         JOFF = JOFF + NOCCI*NOCCN
         KOFF = KOFF + NOCCM*NOCCI
C
 90   CONTINUE
C     
      RETURN
      END
