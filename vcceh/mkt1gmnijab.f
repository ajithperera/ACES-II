C
      SUBROUTINE MKT1GMNIJAB(Z, W, TA, TB, POP1, POP2, VRT1, VRT2,
     &                       DISSYZ, DISSYWA, DISSYWB, NUMSYZ, NUMSYWA,
     &                       NUMSYWB, NTAA, NTBB, LISTWA, LISTWB, LISTZ,
     &                       IRREPWAR, IRREPWBR, IRREPGL, IRREPGR, 
     &                       IRREPX, IUHF, TMP)
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER DISSYZ,DISSYWA,DISSYWB,DIRPRD,POP1,POP2,VRT1,VRT2
      DIMENSION Z(DISSYZ,NUMSYZ),W(DISSYWA,1),TA(NTAA),TB(NTBB)
      DIMENSION POP1(8),POP2(8),VRT1(8),VRT2(8),IWOFF(8)
      DIMENSION TMP(1)
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
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
C
C Common blocks used in the quadratic term
C
      DATA AZERO,ONE /0.0D0,1.0D0/
C
      CALL ZERO(Z, DISSYZ*NUMSYZ)
C
C Get Hbar(Mn,Ie) from (LISTWA, LISTWB)
C       
      IF (MIN (NUMSYWA, DISSYWA) .NE. 0) THEN
C     
         CALL GETLST(W, 1, NUMSYWA, 2, IRREPWAR, LISTWA)
C     
C Perform multiplication 
C     
         JOFFZ = 1
         IOFF  = 1
         IWOFF(1) = 1
C
         DO 5000 IRREPE = 2, NIRREP
C
            IRREP  = IRREPE - 1
            IRREPI = DIRPRD(IRREP, IRREPWAR)
C
            IWOFF(IRREPE) = IWOFF(IRREPE - 1) + POP1(IRREPI)*
     &                      VRT2(IRREP)
C
 5000    CONTINUE
C     
         DO 90 IRREPJ = 1, NIRREP
C
            IRREPE = DIRPRD(IRREPJ, IRREPX)
            IRREPI = DIRPRD(IRREPE, IRREPWAR)
C
            NOCCI = POP1(IRREPI)         
            NOCCJ = POP2(IRREPJ)
            NVRTE = VRT2(IRREPE)
C
            JOFFW = IWOFF(IRREPE)
C     
            IF(NVRTE .EQ. 0 .OR. NOCCJ .EQ. 0 .OR. NOCCI .EQ. 0) 
     &      GO TO 80
C     
            CALL XGEMM('N', 'N', DISSYWA*NOCCI, NOCCJ, NVRTE, ONE, 
     &                  W(1,JOFFW), DISSYWA*NOCCI, TB(IOFF), NVRTE,
     &                  AZERO, Z(1,JOFFZ), DISSYZ*NOCCI)
 80        CONTINUE
C
           JOFFZ = JOFFZ + NOCCJ*NOCCI
           IOFF  = IOFF +  NOCCJ*NVRTE
C
 90      CONTINUE
C
      ENDIF
C     
      IF (MIN (NUMSYWA, DISSYWA) .NE. 0) THEN
C     
C In RHF, the second term is transposed of the first
C
         IF(IUHF .EQ. 0) THEN
            CALL QSYMRHF(IRREPGL, IRREPGR, POP1, POP1, DISSYZ, Z, TMP,
     &                    TMP(1 + DISSYZ), TMP(1 + 2*DISSYZ))
         ENDIF
C
C The product is G(Mn,Ij) and update (MN,IJ) lists.
C
         CALL GETLST(W, 1, NUMSYZ, 1, IRREPGR, LISTZ) 
         CALL VADD(W, W, Z, NUMSYZ*DISSYZ, ONE)
         CALL PUTLST(W, 1, NUMSYZ, 1, IRREPGR, LISTZ)
      ENDIF
C
      IF (IUHF .EQ. 1) THEN
C     
C For UHF, the second term must be calculated explicitly.
C     
         IF (MIN(NUMSYWB, DISSYWB) .NE. 0) THEN
C     
            CALL ZERO(Z, DISSYZ*NUMSYZ)
C     
C Get Hbar(Mn,Aj) form (LISTWB)
C     
            CALL GETLST(W, 1, NUMSYWB, 2, IRREPWBR, LISTWB)
            CALL SYMTR1(IRREPWBR, VRT1, POP2, DISSYWB, W, TMP,
     &                  TMP(1 + DISSYWB), TMP(1 + 2*DISSYWB))
C     
            JOFFZ = 1
            IOFF  = 1
            IWOFF(1) = 1
C
            DO 5100 IRREPE = 2, NIRREP
C
               IRREP  = IRREPE - 1
               IRREPI = DIRPRD(IRREP, IRREPWBR)
C
               IWOFF(IRREPE) = IWOFF(IRREPE - 1) + POP2(IRREPI)*
     &                         VRT1(IRREP)
C
 5100       CONTINUE
C     
            DO 190 IRREPI = 1, NIRREP

               IRREPE = DIRPRD(IRREPI, IRREPX)
               IRREPJ = DIRPRD(IRREPE, IRREPWBR)
C
               NOCCI = POP1(IRREPI)
               NOCCJ = POP2(IRREPJ)
               NVRTE = VRT1(IRREPE)
C
               JOFFW = IWOFF(IRREPE)
C     
               IF (NOCCI .EQ. 0 .OR. NOCCJ .EQ. 0 .OR. NVRTE .EQ. 0)
     &         GO TO 180
C     
               CALL XGEMM('N', 'N', DISSYWB*NOCCJ, NOCCI, NVRTE, ONE,
     &                     W(1,JOFFW), DISSYWB*NOCCJ, TA(IOFF), NVRTE,
     &                     ONE, Z(1, JOFFZ), DISSYZ*NOCCJ)
C     
 180           CONTINUE
C
               JOFFZ = JOFFZ + NOCCJ*NOCCI
               IOFF  = IOFF  + NOCCI*NVRTE
C
 190        CONTINUE
C
            CALL SYMTR1(IRREPGR, POP2, POP1, DISSYZ, Z, TMP,
     &                  TMP(1 + DISSYZ), TMP(1 + 2*DISSYZ))
C
C The product is G(Mn,Ij) and update (MN,IJ) lists.
C
            CALL GETLST(W, 1, NUMSYZ, 1, IRREPGR, LISTZ)
            CALL VADD(W, W, Z, NUMSYZ*DISSYZ, ONE)
            CALL PUTLST(W, 1, NUMSYZ, 1, IRREPGR, LISTZ)
         ENDIF
C
      ENDIF  
C
      RETURN
      END
