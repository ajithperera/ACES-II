      SUBROUTINE ANALEXCI(IUHF)
C
C  THE EXCITED STATE ENERGY SPECTRUM OF THE ELECTRO-ATTACHED SYSTEM
C  IS CALCULATED AND THE ENERGY DIFFERENCES ARE PRINTED IN A 
C  VARIETY OF UNITS
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER *5 SPIN(2)
      LOGICAL SPECIAL
C
      COMMON/ROOTS/EIGVAL(100,8,3), OSCSTR(100,8,3)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/EAINFO/NUMROOT(8,3)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA SPIN /'ALPHA', 'BETA '/
      DATA FACTEV, FACTCM /27.2113957D0,  2.19474625D5/
C
C  FIRST FIND THE MINIMUM ELECTO-AFFINITY
C
C      irrepG = 2
C      write(6,*) 'The ground state is searched for in irrep', irrepg
      EMIN = 10.0
      DO 10 ISPIN = 1 , 1+IUHF
         DO 20 IRREP = 1, NIRREP
            DO 30 IROOT = 1, NUMROOT(IRREP, ISPIN)
            IF (EIGVAL(IROOT, IRREP, ISPIN) . LT. EMIN) THEN
               EMIN = EIGVAL(IROOT, IRREP, ISPIN)
               MINROOT = IROOT
               MINREP = IRREP
               MINSPIN = ISPIN
            ENDIF
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE
C
      WRITE(6, *) '   ',('*', i=1,70)
      WRITE(6,*) '          GROUND STATE OF ELECTRO-ATTACHED SYSTEM '
      WRITE(6,1000) SPIN(MINSPIN), MINREP, MINROOT
 1000 FORMAT(3X, ' SPIN = ',A5,',', 5X, ' SYMMETRY BLOCK = ',I2,',',
     $   5X,'ROOT ', I3)
      WRITE(6,1010) EMIN, EMIN*FACTEV, EMIN*FACTCM
 1010  FORMAT(3X, ' ELECTRO-AFFINITY :', F10.6, ' AU', 4X,
     $   F10.6, ' eV',4x, F12.3, ' cm(-1)')
C
       CALL GETREC(20, 'JOBARC', 'TOTENERG', IINTFP, ECC)
       WRITE(6,*)'   TOTAL ELECTRONIC ENERGY ', ECC + EMIN
C
       WRITE(6,*)'       ----------------------------------------------'
 1300  FORMAT(/)
       WRITE(6,*)'       EXCITATION ENERGIES OF ELECTRO-ATTACHED SYSTEM'
       WRITE(6, 1015)
 1015  FORMAT(5X, 'ROOT',12X, 'A.U.',16X, 'eV',17x, 'cm(-1)')
       DO 110 ISPIN = 1, 1 + IUHF
          DO 120 IRREP = 1, NIRREP
             SPECIAL = ((ISPIN .EQ. MINSPIN) .AND. (IRREP.EQ.MINREP))
             NMIN = 0
             IF (SPECIAL) NMIN = 1
             IF (NUMROOT(IRREP,ISPIN) .GT. NMIN) THEN
                WRITE(6,1300)
                WRITE(6,1020) SPIN(ISPIN), IRREP
 1020           FORMAT(3X, A5, ' SPIN,   SYMMETRY BLOCK: ', I2)
                WRITE(6,*) '  ---------------------------------'
             ENDIF
             DO 130 IROOT = 1, NUMROOT(IRREP, ISPIN)
                IF (.NOT. (SPECIAL .AND. (IROOT.EQ.MINROOT))) THEN
                   XAU = EIGVAL(IROOT, IRREP, ISPIN) - EMIN
                   XEV = XAU * FACTEV
                   XCM = XAU * FACTCM
                   WRITE(6,1030) IROOT, XAU, XEV, XCM
                ENDIF
 130         CONTINUE
 120      CONTINUE
 110   CONTINUE
C
 1030  FORMAT(3X, I4, 5X, F15.8, 5X, F15.8, 5X, F16.4)
       WRITE(6,1300)
      WRITE(6, *) '   ',('*', i=1,70)
       WRITE(6,1300)
C
       RETURN
       END
