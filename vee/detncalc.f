      SUBROUTINE DETNCALC(SCR,IUHF)
C
C THE NUMBER OF DIFFERENT EXCITATION SPECTRA TO BE CALCULATED
C IS DETERMINED. AT PRESENT NCALC = 1, UNLESS FOR CORE-EXCITATION
C SPECTRA. THEN IT EQUALS THE NUMBER OF DIFFERENT 'CORE-SITES'
C IN THE SYSTEM
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION SCR(*), EMIN, EMINPREV, THRESH
C
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/PROJECT/IPROJECT, IPATTERN, NCALC, ICALC, IWINDOW(8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
C
      CORE_WINDOW = 0
      DO IRREP = 1, NIRREP
         CORE_WINDOW=CORE_WINDOW + IWINDOW(IRREP)
      ENDDO

      IF (IPATTERN .NE. 1) THEN
        NCALC = 1
      ELSE
C      
C  GET ORBITAL EIGENVALUES TO DETERMINE RELEVANT ORBITALS
C
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
        IF (IUHF. NE. 0) THEN
          CALL GETREC(20, 'JOBARC', 'SCFEVALB', IINTFP*NBASB, SCR(I010))
        ENDIF

        IF (CORE_WINDOW .NE. 0) THEN
           NCALC=1
        ELSE
C  
C  DETERMINE NUMBER OF DIFFERENT CORE-SITES IN MOLECULE
C
        THRESH = 2.0D0 / 27.2113957D0
C
C  CORE ENERGIES CORRESPOND TO SAME 'SITE' IF THEY LY WITHIN THRESH (2 EV)
C
        DO ISPIN = 1, IUHF+1
        IF (ISPIN .EQ. 1) IOFF = I000 - 1
        IF (ISPIN .EQ. 2) IOFF = I010 - 1

        NCALC = 0
        EMINPREV = - 1.0D8

  100   EMIN = 100.0D0
        DO I = 1, NOCCO(ISPIN)
          IF (SCR(IOFF+I) .GT. EMINPREV .AND. SCR(IOFF+I) .LT. EMIN)
     &       EMIN = SCR(IOFF+I)
        ENDDO
        IF (EMIN .LT. -4.0) THEN
C
C   CORE-LEVEL ENERGY IS FOUND
C
          NCALC = NCALC + 1
          EMINPREV = EMIN + THRESH
C
        ENDIF
C
C  FIND NEXT CORE-SITE
C
        IF (EMIN .LT. -4.0) GOTO 100
C
        IF (ISPIN .EQ. 2) THEN
        WRITE(6,"(a,I3)")'  Number of different core-sites: ', NCALC
        Write(6,"(a,a)") "  Number of core-sites are determined by",
     +                   " counting the number of consecutive"
        Write(6,"(a,a)") "  eigenvalues that are within 2 eV. Note", 
     +                   " that counting start start from the"
        Write(6,"(a)")   "  lowest eigenvalue."
        Write(6,*)
        ENDIF 
        ENDDO 
        ENDIF
C
      ENDIF
C
      RETURN
      END
