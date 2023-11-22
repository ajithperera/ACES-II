C
      SUBROUTINE MKT1GMNIJAA(W, Z, T, POP, VRT, DISSYZ, DISSYW, 
     &                       NUMSYZ, NUMSYW, NOCCSQ, NFT, LISTW,
     &                       LISTZ, IRREPWR, IRREPGR, IRREPX, ISPIN, 
     &                       TMP)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DISSYW,DISSYZ,DIRPRD,POP,VRT
      CHARACTER*8 SPCASE(2)
C
      DIMENSION Z(DISSYZ,NOCCSQ),W(DISSYW,NUMSYW),T(NFT),
     &          POP(8),VRT(8),TMP(1), IWOFF(8)
C
      COMMON /SYMINF/ NSTART,NIRREP,IRREPBA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /FILES/ LUOUT, MOINTS
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
      DATA SPCASE /'AAAA =  ', 'BBBB =  '/
      DATA AZERO, ONE, ONEM /0.0D0, 1.0D0, -1.D0/
C
      CALL ZERO(Z, DISSYZ*NOCCSQ)
C
C Get integrals Hbar(MN,IE) from (LISTW1, LISTW2)
C
      IF (MIN (NUMSYW, DISSYW) .NE. 0) THEN
C
         CALL GETLST(W, 1, NUMSYW, 2, IRREPWR, LISTW)
C
C Set up the offset arrays for T1 amplitudes and Hbar elements
C and Perform multiplication,
C     
         JOFFZ = 1
         IOFF  = 1 
         IWOFF(1) = 1
C
         DO 5000 IRREPE = 2, NIRREP
C
            IRREP  = IRREPE - 1
            IRREPI = DIRPRD(IRREP, IRREPWR)
C
            IWOFF(IRREPE) = IWOFF(IRREPE - 1) + POP(IRREPI)*VRT(IRREP)
C
 5000    CONTINUE
C
         DO 90 IRREPJ = 1, NIRREP
C
            IRREPE = DIRPRD(IRREPJ, IRREPX)
            IRREPI = DIRPRD(IRREPE, IRREPWR)
C          
            NOCCJ = POP(IRREPJ)
            NVRTE = VRT(IRREPE)
            NOCCI = POP(IRREPI)
C
            JOFFW = IWOFF(IRREPE)
C     
            IF(NVRTE .EQ. 0 .OR. NOCCJ .EQ. 0 .OR. NOCCI .EQ. 0)
     &      GO TO 80
C     
            CALL XGEMM('N', 'N', DISSYW*NOCCI, NOCCJ, NVRTE, ONE, 
     &                  W(1,JOFFW), DISSYW*NOCCI, T(IOFF), NVRTE, 
     &                  AZERO, Z(1,JOFFZ), DISSYZ*NOCCI)
 80         CONTINUE
C
            JOFFZ = JOFFZ + NOCCJ*NOCCI
            IOFF  = IOFF  + NOCCJ*NVRTE
C
 90      CONTINUE
C
      CALL ASSYM(IRREPGR, POP, DISSYZ, DISSYZ, W, Z)
C
C The resulting product is G(MN,IJ) in the fist call and G(mn,ij)
C in the second call. Now update the (MN,IJ) list in the disk.
C
      CALL GETLST(Z, 1, NUMSYZ, 1, IRREPGR, LISTZ)
      CALL VADD(Z, Z, W, NUMSYZ*DISSYZ, ONE)
      CALL PUTLST(Z, 1, NUMSYZ, 1, IRREPGR, LISTZ)
C     
C Update the doubles contribution in the disk
C     
         IF (IFLAGS(1) .GE. 20) THEN
C
            NSIZE = NUMSYZ*DISSYZ
            CALL HEADER('Checksum @-T1QGMNIJ per sym. block', 0, LUOUT)
C            
            WRITE(LUOUT, *) SPCASE(ISPIN), SDOT(NSIZE, Z, 1, Z, 1)
         ENDIF
C
      ENDIF
C
      RETURN
      END 
