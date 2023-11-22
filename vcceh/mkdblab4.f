C
      SUBROUTINE MKDBLAB4(T2, Z2, FMIAA, FMIBB, POP1, POP2, VRT1, 
     &                    VRT2, DISSYT, DISSYZ, NUMSYT, NUMSYZ,
     &                    NFSIZA, NFSIZB, LISTT, LISTZ, IRREPTL,
     &                    IRREPTR, IRREPQL, IRREPQR, IRREPX,
     &                    IUHF, TMP)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DISSYT,DISSYZ,DIRPRD,POP1,POP2,VRT1,VRT2
      DIMENSION T2(DISSYT,NUMSYT),Z2(DISSYZ,NUMSYZ),FMIAA(NFSIZA),
     &          FMIBB(NFSIZB)
      DIMENSION TMP(1),POP1(8),POP2(8),VRT1(8),VRT2(8),ITOFF1(8),
     &          ITOFF2(8)
C
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
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
      DATA AZERO,ONE,ONEM /0.0D0,1.0D0,-1.0D0/
C
C Get T2(beta) amplitudes
C
      CALL GETLST(T2, 1, NUMSYT, 2, IRREPTR, LISTT)
      CALL ZERO(Z2, NUMSYZ*DISSYZ)
C
C Perform multiplication
C
      JOFF = 1     
      KOFF = 1
      ITOFF1(1) = 1 
C
      DO 5000 IRREPM = 2, NIRREP
C
         IRREP = IRREPM - 1
         IRREPI = DIRPRD(IRREP, IRREPTR)
         ITOFF1(IRREPM) = ITOFF1(IRREPM - 1) + POP2(IRREP)*POP1(IRREPI)
C
 5000 CONTINUE
C
      DO 90 IRREPJ = 1, NIRREP
C
         IRREPIJ = IRREPJ
         IRREPIM = DIRPRD(IRREPIJ, IRREPX)
         IRREPTM = IRREPIM
         IRREPTI = DIRPRD(IRREPTM, IRREPTR)
C
         NOCCIM = POP2(IRREPIM)
         NOCCIJ = POP2(IRREPIJ)
         NOCCTM = POP2(IRREPTM)
         NOCCTI = POP1(IRREPTI)
C
         IOFF = ITOFF1(IRREPTM)
C
         IF (NOCCIM .GT.0 .AND. NOCCTI .GT. 0) THEN
C
            CALL XGEMM('N', 'N', DISSYT*NOCCTI, NOCCIJ, NOCCIM, ONEM,
     &                  T2(1,IOFF), NOCCTI*DISSYT, FMIBB(JOFF),
     &                  NOCCIM, AZERO, Z2(1,KOFF), NOCCTI*DISSYZ)
         ENDIF
C
         KOFF = KOFF + NOCCIJ*NOCCTI
         JOFF = JOFF + NOCCIJ*NOCCIM
C
 90   CONTINUE
C
      IF (IUHF .EQ. 0) THEN
C
C In RHF this is simply a transposition, update the doubles contribution
C in the disk.
C
         CALL GETLST(T2, 1, NUMSYZ, 1, IRREPQR, LISTZ)
         CALL SYMRHF(IRREPQR, VRT1, POP1, DISSYZ, Z2, TMP,
     &               TMP(1 + DISSYZ), TMP(1 + 2*DISSYZ))
         CALL VADD(T2, T2, Z2, NUMSYZ*DISSYZ, ONE)
         CALL PUTLST(T2, 1, NUMSYZ, 1, IRREPQR, LISTZ)
C    
      ELSE
C
C If UHF we have to calculate the second part of the total contribution.
C
         CALL SYMTR1(IRREPTR, POP1, POP2, DISSYT, T2, TMP,
     &               TMP(1 + DISSYT), TMP(1 + 2*DISSYT))
         CALL SYMTR1(IRREPQR, POP1, POP2, DISSYZ, Z2, TMP,
     &               TMP(1 + DISSYZ), TMP(1 + 2*DISSYZ))
C
         JOFF = 1     
         KOFF = 1
         ITOFF2(1) = 1
C
         DO 5100 IRREPM = 2, NIRREP
C
            IRREP = IRREPM - 1
            IRREPJ = DIRPRD(IRREP, IRREPTR)
            ITOFF2(IRREPM) = ITOFF2(IRREPM - 1) + POP1(IRREP)*
     &                       POP2(IRREPJ)
C
 5100    CONTINUE
C
         DO 190 IRREPI = 1, NIRREP
C
            IRREPII = IRREPI
            IRREPIM = DIRPRD(IRREPII, IRREPX)
            IRREPTM = IRREPIM
            IRREPTJ = DIRPRD(IRREPTM, IRREPTR)
C
            NOCCIM = POP1(IRREPIM)
            NOCCII = POP1(IRREPII)
            NOCCTJ = POP2(IRREPTJ)
            NOCCTM = POP1(IRREPTM)
C
            IOFF = ITOFF2(IRREPTM)
C     
            IF (NOCCIM .GT.0 .AND. NOCCTJ .GT. 0) THEN
C
               CALL XGEMM('N', 'N', DISSYT*NOCCTJ, NOCCII, NOCCIM, 
     &                     ONEM, T2(1,IOFF), NOCCTJ*DISSYT,
     &                     FMIAA(JOFF), NOCCIM, ONE, Z2(1,KOFF),
     &                     NOCCTJ*DISSYZ)
            ENDIF
C
            KOFF = KOFF + NOCCTJ*NOCCII
            JOFF = JOFF + NOCCII*NOCCIM
C
 190     CONTINUE
C
C Update the doubles list in the disk
C
         IF (IFLAGS(1) .GE. 20) THEN
            NSIZE = DISSYZ*NUMSYZ
            CALL HEADER('Checksum @-IMIINSD Doubles', 0, LUOUT)
            WRITE(LUOUT, *) 'ABAB =  ', SDOT(NSIZE, Z2, 1,
     &                      Z2, 1)
         ENDIF
C
         CALL SYMTR1(IRREPQR, POP2, POP1, DISSYZ, Z2, TMP, 
     &               TMP(1 + DISSYZ), TMP(1 + 2*DISSYZ))
C
         CALL GETLST(T2, 1, NUMSYZ, 1, IRREPQR, LISTZ)
         CALL VADD(Z2, Z2, T2, NUMSYZ*DISSYZ, ONE)
         CALL PUTLST(Z2, 1, NUMSYZ, 1, IRREPQR, LISTZ)
C
      ENDIF
C
      RETURN
      END
