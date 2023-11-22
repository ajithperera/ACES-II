      SUBROUTINE CALCEXCP(IUHF, SCR, MAXCOR, IRREPX, 
     &   DOUBLE, NONSTD, NSIZEC, IOPT)
C
C  THIS SUBROUTINE CALCULATES EXCITATION PATTERNS OR A 'MASK'
C
C     IPATTERN INDICATES THE TYPE OF EXCITATION PATTERN. 
C     THESE PATTERNS CAN BE EXTENDED OVER TIME. CURRENTLY IMPLEMENTED ARE
C
C         0   NO PATTERN
C         1   CORE:  EXCITATIONS INCLUDE AT LEAST ONE CORE-ORBITAL
C         2   LUMO:  EXCITATIONS INTO THE LOWEST UNOCCUPIED ORBITAL
C         3   HOMO:  EXCITATIONS FROM THE HIGHEST OCCUPIED ORBITAL
C
C   IOPT YIELDS DIFFERENT OPTIONS. CURRENTLY ONLY THE CORE VARIANT HAS 
C               MORE OPTIONS
C    IOPT = 1:  DETERMINE SPECIFIC STATES OF INTEREST; MORE STRICT PROJECTION.
C    IOPT = 2:  DETERMINE A 'LOOSE' PATTERN, ALL STATES THAT MAY CONTRIBUTE 
C               SIGNIFICANTLY
C               THE LATTER OPTION MAY EQUAL THE ACTUAL PROJECTION.
C
C IOPT=1,2 combination can be used to built a strict projection, but
C it is not used. Most likely Marcel is experimenting.

      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION SCR(MAXCOR), EMIN, EMAX, FACTOR,
     &   THRESH, ONE, EMINPREV
      DIMENSION IOFFPOP(8,2), IOFFVRT(8,2), INDEX(2)
      LOGICAL DOUBLE,NONSTD
      logical print
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/PROJECT/IPROJECT, IPATTERN, NCALC, ICALC, IWINDOW(8)
      COMMON/FLAGS/IFLAGS(100)
C
      ONE = 1.0D0
      CORE_WINDOW = 0
      DO IRREP = 1, NIRREP
         CORE_WINDOW=CORE_WINDOW + IWINDOW(IRREP)
      ENDDO
      I000 = 1
C
      IF (IPATTERN .EQ. 0) THEN
C
C ALL EXCITATIONS ARE TO BE INCLUDED
C
        DO I = 1, NSIZEC
          SCR(I) = ONE
        ENDDO
        CALL UPDATES(IRREPX,SCR(I000),444,0,490,IUHF)
      ELSE
C
C  CONSIDER VARIOUS CASES.
C
      CALL ZERO(SCR(I000), NSIZEC)
      CALL UPDATES(IRREPX,SCR(I000),444,0,490,IUHF)
C      
C  GET ORBITAL EIGENVALUES TO DETERMINE RELEVANT ORBITALS
C
      print = .false.
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
      FACTOR = 1.0D0
C
      IF (IPATTERN .EQ. 1) THEN

        IF (CORE_WINDOW .EQ. 0 .AND. .NOT. NONSTD) THEN 
C
C  CORE-PATTERN: DETERMINE SPECIFIC CORE-ENERGY CORRESPONDING TO ICALC
C  (SEE ALSO DETNCALC)
C
        THRESH = 2.0D0 / 27.2113957D0
C
C  CORE ENERGIES CORRESPOND TO SAME 'SITE' IF THEY LY WITHIN THRESH (2 EV)
C
        DO ISPIN = 1, 1 + IUHF
        EMINPREV = - 1.0D8
        JCALC = 0
        IF (ISPIN .EQ. 1) IOFF = I000-1
        IF (ISPIN .EQ. 2) IOFF = I010-1

  100   EMIN = 100.0D0
        DO I = 1, NOCCO(ISPIN)
          IF (SCR(IOFF+I) .GT. EMINPREV .AND. SCR(IOFF+I) .LT. EMIN)
     &       EMIN = SCR(IOFF+I)
        ENDDO
        write(6,"(2I2,1x,F15.10)") jcalc,icalc, Emin
        IF (EMIN .LT. -4.0) THEN
C
C   CORE-LEVEL ENERGY IS FOUND
C
          JCALC = JCALC + 1

          IF (JCALC .LT. ICALC) THEN
            EMINPREV = EMIN + THRESH
          ENDIF 
 
        ENDIF
C
C  FIND NEXT CORE-SITE
C
        write(6,"(2I2,1x,F15.10)") jcalc,icalc, Emin
        IF (JCALC .LT. ICALC) GOTO 100
C
C  OCCUPIED ORBITALS NEAR EMIN ARE OF 'INTEREST', THIS DEPENDS ON IOPT
C
        IF (IOPT .EQ. 2) THEN
          THRESH = 2.0D0 / 27.2113957D0
        ELSE
          THRESH = 2.0D0 / 27.2113957D0
        ENDIF
C          
        IF (IOPT .EQ. 1 .AND. ISPIN .EQ. 1) THEN
        write(6,*)
        write(6,"(a,a)") '  Core excited state are searched using',
     &                   ' excitation mask' 
        WRITE(6,*) ' Core-orbitals included in excitation mask '
        ENDIF 
        Write(6,*)

C Moved to above so different Alpha and Beta blocks can be handled 
C without assuming that they are identical.

CSS        DO ISPIN = 1, 1 + IUHF

          IF (ISPIN .EQ. 1 .AND. IOPT .EQ. 1) Write(6,"(a)") 
     &                     " Alpha block"
          IF (ISPIN .EQ. 2 .AND. IOPT .EQ. 1) Write(6,"(a)") 
     &                     " Beta  block"
          DO IRREP = 1, NIRREP
            DO I = 1, POP(IRREP, ISPIN)
              IF ((SCR(IOFFPOP(IRREP,ISPIN)+I-1) - EMIN)
     &           .LT. THRESH) THEN
                INDEX(ISPIN) = I
                IF (IUHF .EQ. 0) THEN
                  INDEX(3-ISPIN) = I
                ELSE
                  INDEX(3-ISPIN) = 0
                ENDIF
                IF (IOPT .EQ. 1) THEN
                write(6,1000) I, IRREP, SCR(IOFFPOP(IRREP,ISPIN)+I-1) 
1000            Format(2x, I4, " [",I1,"]",F12.4)
                ENDIF
                CALL UPDEXCP(IUHF, SCR(I020), MAXCOR-I020+1,IRREPX,
     &             IRREP, INDEX, 'HOLE', FACTOR, DOUBLE)
              ENDIF
            ENDDO
          ENDDO
        ENDDO

        ELSEIF (CORE_WINDOW .GT. 0 .AND. .NOT. NONSTD) THEN

        CALL DRIVE_CORE_TDA_WINDOW(SCR,MAXCOR,IRREPX,IUHF)
        IF (IOPT .EQ. 1) THEN
        write(6,*)
        write(6,"(a,a)") '  Core excited state are searched using',
     &                   ' the core-window defined by the user'
        Write(6,"(a)")   '  input'
        ENDIF 
        ENDIF 
C
      ELSEIF (IPATTERN .EQ. 2) THEN
C
C  EXCITATIONS INTO THE LUMO. FIRST FIND THE LUMO-ENERGY
C
        EMIN = 100.0D0
        DO ISPIN = 1, 1 + IUHF
          DO IRREP = 1, NIRREP
            DO A = 1, VRT(IRREP, ISPIN)
              IF ( SCR(IOFFVRT(IRREP,ISPIN) + A-1) .LT. EMIN) THEN
                EMIN = SCR(IOFFVRT(IRREP,ISPIN) + A-1) 
              ENDIF
            ENDDO
          ENDDO
        ENDDO
C
C  NOW INCLUDE ALL ORBITALS THAT HAVE THIS ENERGY INTO EXCITATION PATTERN        
C
        THRESH = 0.01D0
C
        DO ISPIN = 1, 1 + IUHF
          DO IRREP = 1, NIRREP
            DO A = 1, VRT(IRREP, ISPIN)
              IF (ABS( SCR(IOFFVRT(IRREP,ISPIN) + A-1) - EMIN) .LT.
     &           THRESH) THEN
                INDEX(ISPIN) = A
                IF (IUHF .EQ. 0) THEN
                  INDEX(3-ISPIN) = A
                ELSE
                  INDEX(3-ISPIN) = 0
                ENDIF
                CALL UPDEXCP(IUHF, SCR(I020), MAXCOR-I020+1,IRREPX,
     &             IRREP, INDEX, 'PART',FACTOR,DOUBLE)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
C
      ELSEIF (IPATTERN .EQ. 3) THEN
C
C  EXCITATIONS FROM THE HOMO. FIRST FIND THE HOMO-ENERGY
C
        EMAX = -100.0D0
        DO ISPIN = 1, 1 + IUHF
          DO IRREP = 1, NIRREP
            DO I = 1, POP(IRREP, ISPIN)
              IF ( SCR(IOFFPOP(IRREP,ISPIN) + I-1) .GT. EMAX) THEN
                EMAX = SCR(IOFFPOP(IRREP,ISPIN) + I-1) 
              ENDIF
            ENDDO
          ENDDO
        ENDDO
C
C  NOW INCLUDE ALL ORBITALS THAT HAVE THIS ENERGY INTO EXCITATION PATTERN        
C
        THRESH = 0.01D0
C
        DO ISPIN = 1, 1 + IUHF
          DO IRREP = 1, NIRREP
            DO I = 1, POP(IRREP, ISPIN)
              IF (ABS( SCR(IOFFPOP(IRREP,ISPIN) + I-1)-EMAX) .LT.
     &           THRESH) THEN
                INDEX(ISPIN) = I
                IF (IUHF .EQ. 0) THEN
                  INDEX(3-ISPIN) = I
                ELSE
                  INDEX(3-ISPIN) = 0
                ENDIF
                CALL UPDEXCP(IUHF, SCR(I020), MAXCOR-I020+1,IRREPX,
     &             IRREP, INDEX, 'HOLE', FACTOR, DOUBLE)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
C
C  END IF'S OVER DIFFERENT EXCITATION PATTERNS
C
      ENDIF
C
      ENDIF
C
C  EXCITATION PATTERNS ARE DETERMINED AND PUT ON LISTS 444-446, 490
C  NOW NORMALIZE THEM SUCH THAT MASK IS EITHER ONE OR ZERO
C
      CALL LOADVEC1(IRREPX,SCR,MAXCOR,IUHF,490,0,443,NSIZEC,
     &              .FALSE.)
      DO I = 1, NSIZEC
        IF (SCR(I) .GT. 0.75) THEN
          SCR(I) = 1.0D0
        ELSE
          SCR(I) = 0.0D0
        ENDIF
      ENDDO
C
C  PUT NORMALIZED PROJECTION VECTORS BACK ON LIST
C
      CALL UPDATES(IRREPX,SCR(I000),444,0,490,IUHF)
C
      print = (iflags(1) .ge. 50)
      if (print) then
       CALL PRMAINX(IUHF, SCR, MAXCOR, IRREPX,1)
      endif
C
      RETURN
      END
