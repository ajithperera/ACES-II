C
      SUBROUTINE MKDBLAB3(T2, Z2, T, Z, FEAAA, FEABB, POP1, POP2, VRT1,
     &                    VRT2, DISSYT, DISSYZ, NUMSYT, NUMSYZ, NFSIZA,
     &                    NFSIZB, LISTT, LISTZ, IRREPTL, IRREPTR, 
     &                    IRREPQL, IRREPQR, IRREPX, IUHF, TMP)
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER DISSYT,DISSYZ,DIRPRD,POP1,POP2,VRT1,VRT2
      DIMENSION T2(DISSYT,NUMSYT),Z2(DISSYZ,NUMSYZ),FEAAA(NFSIZA),
     &           FEABB(NFSIZB),T(NUMSYT,DISSYT),Z(NUMSYZ,DISSYZ)
      DIMENSION TMP(1),POP1(8),POP2(8),VRT1(8),VRT2(8), ITOFF1(8),
     &           ITOFF2(8)
C
      COMMON/SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON/FLAGS/ IFLAGS(100)
      COMMON/FILES/ LUOUT, MOINTS
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
C Get T2 amplitudes
C
      CALL GETLST(Z2, 1, NUMSYT, 2, IRREPTR, LISTT)
      CALL TRANSP(Z2, T, NUMSYT, DISSYT)
      CALL ZERO(Z, NUMSYZ*DISSYZ)
C
C Perform multiplication
C
      JOFF = 1
      KOFF = 1
      ITOFF1(1) = 1
C
      DO 5000 IRREPE = 2, NIRREP
C
         IRREP = IRREPE - 1
         IRREPAO = DIRPRD(IRREP, IRREPTL)
         ITOFF1(IRREPE) = ITOFF1(IRREPE - 1) + VRT2(IRREP)*
     &                     VRT1(IRREPAO)
C
 5000 CONTINUE
C
      DO 90 IRREPBO = 1, NIRREP
C     
         IRREPIB = IRREPBO
         IRREPIE = DIRPRD(IRREPIB, IRREPX)
         IRREPTE = IRREPIE
         IRREPTA = DIRPRD(IRREPTE, IRREPTL)
C     
         NVRTTE = VRT2(IRREPTE)
         NVRTTA = VRT1(IRREPTA)
         NVRTIE = VRT2(IRREPIE)
         NVRTIB = VRT2(IRREPIB)
C
         IOFF = ITOFF1(IRREPTE)
C     
         IF (NVRTIE .GT. 0 .AND. NVRTTA .GT. 0) THEN
C     
            CALL XGEMM('N', 'N', NUMSYT*NVRTTA, NVRTIB, NVRTIE,
     &                  ONE, T(1,IOFF), NVRTTA*NUMSYT, FEABB(JOFF),
     &                  NVRTIE, AZERO, Z(1, KOFF), NVRTTA*NUMSYZ)
         ENDIF
C     
         KOFF = KOFF + NVRTTA*NVRTIB
         JOFF = JOFF + NVRTIE*NVRTIB
C
 90   CONTINUE
C     
C RHF spin adapted code
C
      IF(IUHF .EQ. 0) THEN
C
         CALL SCOPY(DISSYZ*NUMSYZ, Z, 1, T2, 1)
         CALL TRANSP(T2, Z2, DISSYZ, NUMSYZ)
         CALL GETLST(T2, 1, NUMSYZ, 1, IRREPQR, LISTZ)
         CALL SYMRHF(IRREPQL, VRT1, POP1, DISSYZ, Z2, TMP,
     &               TMP(1 + DISSYZ), TMP(1 + 2*DISSYZ))
C
C Update the existing list in the disk
C
         CALL VADD(Z2, Z2, T2, NUMSYZ*DISSYZ, ONE)
         CALL PUTLST(Z2, 1, NUMSYZ, 1, IRREPQR, LISTZ)
C
      ELSE
C
         CALL SYMTR1(IRREPTL, VRT1, VRT2, NUMSYT, T, TMP, 
     &               TMP(1 + NUMSYT), TMP(1 + 2*NUMSYT))
         CALL SYMTR1(IRREPQL, VRT1, VRT2, NUMSYZ, Z, TMP,
     &               TMP(1 + NUMSYZ), TMP(1 + 2*NUMSYZ))
C     
C Perform second multiplication
C     
         JOFF = 1 
         KOFF = 1
         ITOFF2(1) = 1 
C
         DO 5100 IRREPE = 2, NIRREP
C
            IRREP = IRREPE - 1
            IRREPBO = DIRPRD(IRREP, IRREPTL)
            ITOFF2(IRREPE) = ITOFF2(IRREPE - 1) + VRT1(IRREP)*
     &                       VRT2(IRREPBO)
C
 5100    CONTINUE
C     
         DO 190 IRREPAO = 1, NIRREP
C
            IRREPIA = IRREPAO
            IRREPIE = DIRPRD(IRREPIA, IRREPX)
            IRREPTE = IRREPIE
            IRREPTB = DIRPRD(IRREPTE, IRREPTL)
C     
            NVRTTE = VRT1(IRREPTE)
            NVRTTB = VRT2(IRREPTB)
            NVRTIE = VRT1(IRREPIE)
            NVRTIA = VRT1(IRREPIA)
C
            IOFF = ITOFF2(IRREPTE)
C     
            IF (NVRTTB .GT. 0 .AND. NVRTIE .GT. 0) THEN
C     
               CALL XGEMM('N', 'N', NUMSYT*NVRTTB, NVRTIA, NVRTIE,
     &                     ONE, T(1,IOFF), NVRTTB*NUMSYT, FEAAA(JOFF),
     &                     NVRTIE, ONE, Z(1, KOFF), NVRTTB*NUMSYZ)
            ENDIF
C     
            KOFF = KOFF + NVRTTB*NVRTIA
            JOFF = JOFF + NVRTIE*NVRTIA
C
 190     CONTINUE
C     
         CALL SYMTR1(IRREPQL, VRT2, VRT1, NUMSYZ, Z, TMP,
     &               TMP(1 + NUMSYZ), TMP(1 + 2*NUMSYZ))
C
         CALL SCOPY(DISSYZ*NUMSYZ, Z, 1, T2, 1)
         CALL TRANSP(T2, Z2, DISSYZ, NUMSYZ)
C
         IF (IFLAGS(1) .GE. 20) THEN
            NSIZE = DISSYZ*NUMSYZ
            CALL HEADER('Checksum @-IAEINSD Doubles', 0, LUOUT)
            WRITE(LUOUT, *) 'ABAB =  ', SDOT(NSIZE, Z2, 1, Z2, 1)
         ENDIF
C
C Update the lists
C
         CALL GETLST(T, 1, NUMSYZ, 1, IRREPQR, LISTZ)
         CALL VADD(Z2, Z2, T, NUMSYZ*DISSYZ, ONE)
         CALL PUTLST(Z2, 1, NUMSYZ, 1, IRREPQR, LISTZ)
C
      ENDIF
C
      RETURN
      END
