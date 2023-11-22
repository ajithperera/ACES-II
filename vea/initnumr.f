      SUBROUTINE INITNUMR(SCR, MAXCOR, IUHF)
C
C     IF THE NUMBER OF ROOTS TO BE DETERMINED IS NOT GIVEN (I.E. ALL ARE ZERO) THEN
C     THE NUMBER OF ROOTS TO BE DETERMINED IN THE EA_EOM CALCULATION IS DETERMINED 
C     FROM INSPECTION OF THE SCF ENERGY EIGENVALUES. POSSIBLE DEGENERACIES BETWEEN
C     BLOCKS IS ACCOUNTED FOR.
C
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION SCR,  EA(8,2), EAMIN, DIFFE
      DIMENSION SCR(MAXCOR), IOFFVRT(8,2), IOFFPOP(8,2)
      logical print, EXCICORE
C
      COMMON/SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/EACALC/LEFTHAND, EXCICORE, SINGONLY, DROPCORE
      COMMON/EAINFO/NUMROOT(8,3)
C
      IF (EXCICORE) THEN
         IROOT = 0
         DO IRREP = 1, NIRREP
            DO ISPIN = 1, 3
               IROOT = IROOT + NUMROOT(IRREP, ISPIN)
            ENDDO
         ENDDO
         IF (IROOT .EQ. 0) THEN 
            WRITE(6,*) ' IN CORE-EXCITATION CALCULATIONS NUMROOT',
     $          ' HAS TO BE SPECIFIED'
            CALL ERREX
         ENDIF
         RETURN
      ENDIF
C
      IROOT = 0
      DO IRREP = 1, NIRREP
         DO ISPIN = 1, 1+IUHF
            IROOT = IROOT + NUMROOT(IRREP, ISPIN)
         ENDDO
      ENDDO
      IF (IROOT .GT. 0) THEN
         RETURN
      ELSE
         WRITE(6,*)
         WRITE(6,*) ' NUMBER OF DESIRED ROOTS',
     $      ' IS ESTIMATED FROM ORBITAL EIGENVALUES'
         WRITE(6,*)
      ENDIF
C
      print =.false.
      NBASA = NOCCO(1) + NVRTO(1)
      IF (IUHF.NE.0) THEN
         NBASB = NOCCO(2)  + NVRTO(2)
      ELSE
         NBASB = 0
      ENDIF
      I000 = 1
      I010 = I000 + NBASA
      I020 = I010 + NBASB
      CALL GETREC(20, 'JOBARC', 'SCFEVALA', IINTFP*NBASA, SCR(I000))
      if (print) then
         write(6,*) ' Hartree Fock orbital energies : alfa'
         call output(scr(i000), 1, 1, 1, nbasa, 1, nbasa, 1)
      endif
      IF (IUHF. NE. 0) THEN
         CALL GETREC(20, 'JOBARC', 'SCFEVALB', IINTFP*NBASB, SCR(I010))
      if (print) then
         write(6,*) ' Hartree Fock orbital energies : beta'
         call output(scr(i010), 1, 1, 1, nbasb, 1, nbasb, 1)
      endif
      ENDIF
C
C  CALCULATE OFFSETS IN ORBITAL ENERGY ARRAYS
C
      IOFFPOP(1, 1) = I000
      IOFFPOP(1, 2) = I010
      DO IRREP = 2, NIRREP
         DO ISPIN = 1, 1 + IUHF
            IOFFPOP(IRREP, ISPIN) = IOFFPOP(IRREP-1, ISPIN) +
     $         POP(IRREP-1, ISPIN)
         ENDDO 
      ENDDO 
      IOFFVRT(1, 1) = IOFFPOP(NIRREP, 1) + POP(NIRREP,1)
      IF (IUHF .NE. 0) THEN
         IOFFVRT(1, 2) = IOFFPOP(NIRREP, 2) + POP(NIRREP,2)
      ENDIF
      DO IRREP = 2, NIRREP
         DO ISPIN = 1, 1 + IUHF
            IOFFVRT(IRREP, ISPIN) = IOFFVRT(IRREP-1, ISPIN) +
     $         VRT(IRREP-1, ISPIN)
         ENDDO 
      ENDDO 
C
C   DETERMINE MINIMUM SCF EA PER IRREP
C
      DO IRREP = 1, NIRREP
         DO ISPIN = 1, 1 + IUHF
            EA(IRREP,ISPIN) = SCR(IOFFVRT(IRREP,ISPIN))
         ENDDO 
      ENDDO 
C
      print = .false.
      IF (PRINT) THEN
         WRITE(6,*) ' MINIMUM SCF ELECTRO-AFFINITIES '
         CALL OUTPUT(EA,1,nirrep,1,2,8,2,1)
      ENDIF
C
C      DETERMINE MINIMUM ELECTRO-AFFINITY
C
      EAMIN = 1.d100
      DO  ISPIN = 1, 1 + IUHF
         DO  IRREP = 1, NIRREP
            IF (VRT(IRREP,ISPIN).GT.0) EAMIN=MIN(EAMIN,EA(IRREP,ISPIN))
         ENDDO 
      ENDDO
C      
      if (print) write(6,*) ' EAMIN: ', EAMIN
C
      DO  ISPIN =1, 1+IUHF
         DO  IRREP = 1, NIRREP
            NUMROOT(IRREP,ISPIN) = 0
            IF (VRT(IRREP,ISPIN).GT.0) THEN
            IF (ABS(EAMIN-EA(IRREP,ISPIN)).LT.0.05) THEN
               NUMROOT(IRREP, ISPIN) = 1
            ENDIF
            ENDIF
         ENDDO
      ENDDO
C
C  CHECK FOR DEGENERACIES
C
      DO ISPIN = 1, 1+IUHF
         DO IRREP = 1, NIRREP - 1
            IF (NUMROOT(IRREP, ISPIN) .NE. 0) THEN
               DO JRREP = IRREP + 1, NIRREP
                  IF ((NUMROOT(JRREP, ISPIN) .NE. 0) .AND.
     $               (POP(IRREP, ISPIN) .EQ. POP(JRREP, ISPIN)) .AND.
     $               (VRT(IRREP, ISPIN) .EQ. VRT(JRREP,ISPIN))) THEN
                     DIFFE = 0.0D0
                     DO I = 0, POP(IRREP,ISPIN) - 1
                  DIFFE = DIFFE + ABS(SCR(IOFFPOP(IRREP,ISPIN) + I) -
     $                     SCR(IOFFPOP(JRREP,ISPIN) + I))
                     ENDDO
                     DO I = 0, VRT(IRREP,ISPIN) - 1
                  DIFFE = DIFFE + ABS(SCR(IOFFVRT(IRREP,ISPIN) + I) -
     $                     SCR(IOFFVRT(JRREP,ISPIN) + I))
                     ENDDO
                  IF (ABS(DIFFE) .LT. 0.0001) NUMROOT(JRREP,ISPIN) = 0
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO
C
      RETURN
      END
