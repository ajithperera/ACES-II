C
      SUBROUTINE MKT2IMIAB(W, T, Z, ISPIN, POP1, POP2, VRT1, VRT2,
     &                     DISSYW, DISSYT, NUMSYW, NUMSYT, NOCCSQ, 
     &                     IRREPWL, IRREPTL, IRREPWR, IRREPTR, IRREPX,
     &                     TMP, IOFFSET, IUHF, FACT)
C     
C Calculate the following contributions to the F(MI) intermediate
C depending the value of ISPIN. Also handels RHF spin adapted code.
C
C  Z(MI) =  SUM n SUM (E,f) {T(In,Ef)*Hbar(Mn,Ef)} [ISPIN = 1]
C
C  Z(mi) =  SUM N SUM (e,F) {T(iN,eF)*Hbar(mN,eF)} [ISPIN = 2]
C 
C RHF spin adapted expression
C
C  Z(MI) = SUM n SUM E,f {2*Hbar(Mn,Ef) - Hbar(Mn,Fe)}T(In,Ef)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER DISSYT, DISSYW, DIRPRD,POP1,POP2,VRT1,VRT2
      DIMENSION W(DISSYW,NUMSYW),T(DISSYT,NUMSYT),Z(NOCCSQ)
      DIMENSION TMP(1), POP1(8), POP2(8), VRT1(8), VRT2(8),
     &          IWOFF(8)
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
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
C     
      DATA ONE,ONEM,TWO /1.0D0,-1.D0,2.D0/,HALF /0.5D0/
C     
C Pick up first the relevent T2 and W pieces
C     
      LISTW  = 16
      LISTT  = IAPRT2AB5
C
      CALL GETLST(W, 1, NUMSYW, 2, IRREPWR, LISTW)
C     
C Spin adapted code for RHF
C     
      IF(IUHF .EQ. 0) THEN
         CALL SPINAD1(IRREPWR, POP1, DISSYW, W, TMP, TMP(1 + DISSYW))
      ENDIF
C     
      CALL GETLST(T, 1, NUMSYT, 1, IRREPTR, LISTT)
C     
      IF (IUHF .EQ. 1 .AND. ISPIN .EQ. 1) THEN 
         CALL SYMTR1(IRREPWR, POP1, POP2, DISSYW, W, TMP,
     &               TMP(1 + DISSYW), TMP(1 + 2*DISSYW))
         CALL SYMTR1(IRREPTR, POP1, POP2, DISSYT, T, TMP,
     &               TMP(1 + DISSYT), TMP(1 + 2*DISSYT))
      ENDIF
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
         IWOFF(IRREPM) = IWOFF(IRREPM - 1) + POP1(IRREP)*POP2(IRREPN)
C
 5000 CONTINUE
C
      DO 90 IRREPI = 1, NIRREP
C     
C Determine IRREPN and IRREPM
C
         IRREPN = DIRPRD(IRREPI, IRREPTR)
         IRREPM = DIRPRD(IRREPN, IRREPWR)
C         
         NOCCI = POP1(IRREPI)
         NOCCN = POP2(IRREPN)
         NOCCM = POP1(IRREPM)
C
         IOFF = IWOFF(IRREPM)
C
C If zero nothing to compute
C     
         IF (NOCCN .GT. 0 .AND. NOCCM .GT. 0) THEN
C     
            CALL XGEMM('T', 'N', NOCCM, NOCCI, DISSYW*NOCCN, FACT,
     &                  W(1,IOFF), NOCCN*DISSYW, T(1,JOFF),      
     &                  NOCCN*DISSYT, ONE, Z(KOFF), NOCCM)
         ENDIF
C
C Update the offset for the intermediate Z
C
         JOFF = JOFF + NOCCN*NOCCI
         KOFF = KOFF + POP1(IRREPM)*POP1(IRREPI)
C
 90   CONTINUE
C     
      RETURN
      END
