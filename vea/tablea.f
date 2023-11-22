      SUBROUTINE TABLEA(IUHF)
C
C  THE CALCULATION OF ELECTRON AFFINITIES IS SUMMARIZED
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER *5 SPIN(2)
      LOGICAL  LEFTHAND, EXCICORE, SINGONLY, DROPCORE
C
      COMMON/ROOTS/EIGVAL(100,8,3), OSCSTR(100,8,3)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/EAINFO/NUMROOT(8,3)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/EACALC/LEFTHAND, EXCICORE, SINGONLY, DROPCORE
C
      DATA SPIN /'ALPHA', 'BETA '/
      DATA FACTEV, FACTCM /27.2113957D0,  2.19474625D5/
C
C
      IF (EXCICORE) THEN
         RETURN
      ENDIF
         WRITE(6,1300)
         WRITE(6,*)'    SUMMARY OF ELECTRON-AFFINITY EOM-CC CALCULATION'
         WRITE(6,1300)
         WRITE(6, 2015)
 2015    FORMAT(4X, 'SPIN ','  SYM. BLOCK',4X,' ROOT #',9X, 'EA(eV)',
     $      10X, 'TOTAL ENERGY')
         WRITE(6,*)'  --------------------------------------',
     $      '------------------------------'
         WRITE(6,*)
         CALL GETREC(20,'JOBARC','TOTENERG',IINTFP,ECC)
         DO ISPIN = 1, 1+IUHF
            DO IRREP = 1, NIRREP
               DO  IROOT = 1, NUMROOT(IRREP, ISPIN)
                  EA = EIGVAL(IROOT, IRREP, ISPIN) * FACTEV
                  ET = EIGVAL(IROOT, IRREP, ISPIN) + ECC
                  WRITE(6,2030) SPIN(ISPIN), IRREP,IROOT, EA, ET
               ENDDO
            ENDDO
         ENDDO
         WRITE(6,*)
         WRITE(6,*)'  --------------------------------------',
     $      '------------------------------'
         WRITE(6,1300)
C
 2030    FORMAT(3X,A5,7X,I3,7X, I4, 5X, F15.8, 5X, F15.8)
C
         RETURN
 1300    FORMAT(/)
         END
