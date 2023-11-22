      SUBROUTINE INITNUMR(SCR, IUHF, NUMROOT)
C
C     IF THE NUMBER OF ROOTS TO BE DETERMINED IS NOT GIVEN (I.E. ALL ARE ZERO)
C     THEN THE NUMBER OF ROOTS TO BE DETERMINED IN THE IP_EOM CALCULATION IS
C     DETERMINED  FROM INSPECTION OF THE SCF ENERGY EIGENVALUES. POSSIBLE
C     DEGENERACIES BETWEEN BLOCKS IS ACCOUNTED FOR.
C
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION SCR,  LOWIP(8,2), IPMIN, DIFFE, EMIN
      DIMENSION SCR(*), IOFFPOP(8,2), IOFFVRT(8,2), NUMROOT(NIRREP)
      CHARACTER*10 SEARCH
      logical print
C
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON/SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/INFO/NOCCO(2),NVRTO(2)
C
      IROOT = 0
      DO IRREP = 1, NIRREP
        IROOT = IROOT + NUMROOT(IRREP)
      ENDDO
      IF (IROOT .GT. 0) THEN
C
C  Check if NUMROOT is consistent with IP-SEARCH VARIABLE
C
        IF (IFLAGS2(116) .NE. 1) THEN
          DO IRREP = 1, NIRREP
            NUMROOT(IRREP) =
     &           MIN(NUMROOT(IRREP),POP(IRREP,1))
          ENDDO
        ENDIF
        RETURN
      ELSE
        WRITE(6,*)
        WRITE(6,*) ' NUMBER OF DESIRED ROOTS',
     $     ' IS ESTIMATED FROM HARTREE FOCK ORBITAL EIGENVALUES'
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
C
C THE NUMBER OF ROOTS PER IRREP DEPENDS ON INFORMATION IN FLAG 216 'IP_SEARCH'
C
C     0 : VALENCE IP'S ONLY (IP'S SMALLER THAN 1 A.U. = 27 EV)
C     1 : FIRST IP
C     2 : CORE IP'S  (LARGER THAN 4 A.U.)
C     3 : SHAKEUPS (NUMROOT IS REQUIRED FROM JOBARC)
C     4 : KOOPMANS (ALL PRINCIPAL IP'S)
C
C      IPTYPE = IFLAGS2(116)
C
C ONLY IPTYPE = 0 WORKS AT PRESENT
C
      IPTYPE = 0
      IF (IPTYPE .EQ. 0) THEN
        SEARCH = 'VALENCE   '
      ELSEIF (IPTYPE .EQ. 1) THEN
        SEARCH = 'FIRST IP  '
      ELSEIF (IPTYPE .EQ. 2) THEN
        SEARCH = 'CORE IPS  '
      ELSEIF (IPTYPE .EQ. 3) THEN
        SEARCH = 'SHAKE UP  '
      ELSEIF (IPTYPE .EQ. 4) THEN
        SEARCH = 'PRINCIPALS'
      ENDIF
      write(6,*) ' IP_SEARCH: ', SEARCH
C
      IF (IPTYPE .EQ. 3) THEN
        WRITE(6,*) ' AUTOMATIC CONSTRUCTION OF SHAKEUPS NOT IMPLEMENTED'
        CALL ERREX
      ENDIF
C
C  DETERMINE IRREP OF LOWEST VIRTUAL
C
      EMIN = 100.0D0
      ISPIN = 1
      DO IRREP = 1, NIRREP
        DO A = 1, VRT(IRREP, ISPIN)
          IF ( SCR(IOFFVRT(IRREP,ISPIN) + A-1) .LT. EMIN) THEN
            EMIN = SCR(IOFFVRT(IRREP,ISPIN) + A-1) 
            IRREPEA = IRREP
          ENDIF
        ENDDO
      ENDDO
C
      DO IRREP = 1, NIRREP
        ISPIN = 1
        IROOT = 0
        IOFF = IOFFPOP(IRREP,ISPIN)
        IF (IPTYPE .EQ. 0) THEN            
          DO I = IOFF, IOFF + POP(IRREP,ISPIN) - 1
            IF (ABS(SCR(I)) .LT. 1.50D0) IROOT = IROOT + 1
          ENDDO
        ELSEIF (IPTYPE .EQ. 2) THEN
          DO I = IOFF, IOFF + POP(IRREP,ISPIN) - 1
            IF (ABS(SCR(I)) .GT. 4.0D0) IROOT = IROOT + 1
          ENDDO
        ELSEIF(IPTYPE.EQ.4) THEN
          IROOT = POP(IRREP,ISPIN)
        ENDIF
        IRREPEE =DIRPRD(IRREP,IRREPEA)
        NUMROOT(IRREPEE) = IROOT
      ENDDO
C
      IF (IPTYPE.EQ.1) THEN
        
C   DETERMINE MINIMUM SCF IP PER IRREP
C
        DO IRREP = 1, NIRREP
          DO ISPIN = 1, 1 + IUHF
            LOWIP(IRREP,ISPIN) =
     &         SCR(IOFFPOP(IRREP,ISPIN)+POP(IRREP,ISPIN)-1)
          ENDDO 
        ENDDO 
C
        print = .false.
        IF (PRINT) THEN
          WRITE(6,*) ' MINIMUM SCF IONIZATION POTENTIALS '
          CALL OUTPUT(LOWIP,1,nirrep,1,2,8,2,1)
        ENDIF
C
C      DETERMINE MINIMUM IONIZATION POTENTIALS PER IRREP
C
        IPMIN = 1.d100
        DO  ISPIN = 1, 1 + IUHF
          DO  IRREP = 1, NIRREP
            IF (POP(IRREP,ISPIN).GT.0)
     &         IPMIN = MIN(IPMIN,ABS(LOWIP(IRREP,ISPIN)))
          ENDDO 
        ENDDO
C      
        if (print) write(6,*) ' IPMIN: ', IPMIN
C
        DO  IRREP = 1, NIRREP
          ISPIN = 1
          IRREPEE =DIRPRD(IRREP,IRREPEA)
          NUMROOT(IRREPEE) = 0
          IF (POP(IRREP,ISPIN).GT.0) THEN
          IF (ABS(IPMIN-ABS(LOWIP(IRREP,ISPIN))).LT.0.05) THEN
            NUMROOT(IRREPEE) = 1
          ENDIF
          ENDIF
        ENDDO
      ENDIF
C
C  CHECK FOR DEGENERACIES
C
      ISPIN = 1
      DO IRREP = 1, NIRREP - 1
        IF (NUMROOT(IRREP) .NE. 0) THEN
          DO JRREP = IRREP + 1, NIRREP
            IF ((NUMROOT(JRREP) .NE. 0) .AND.
     $         (POP(IRREP, ISPIN) .EQ. POP(JRREP, ISPIN)) .AND.
     $         (VRT(IRREP, ISPIN) .EQ. VRT(JRREP,ISPIN))) THEN
              DIFFE = 0.0D0
              DO I = 0, POP(IRREP,ISPIN) - 1
                DIFFE = DIFFE + ABS(SCR(IOFFPOP(IRREP,ISPIN) + I) -
     $             SCR(IOFFPOP(JRREP,ISPIN) + I))
              ENDDO
              DO I = 0, VRT(IRREP,ISPIN) - 1
                DIFFE = DIFFE + ABS(SCR(IOFFVRT(IRREP,ISPIN) + I) -
     $             SCR(IOFFVRT(JRREP,ISPIN) + I))
              ENDDO
              IF (ABS(DIFFE) .LT. 0.0001) NUMROOT(JRREP) = 0
            ENDIF
          ENDDO
        ENDIF
      ENDDO
C
      WRITE(6,*) ' NUMBER OF ROOTS PER IRREP TO BE DETERMINED'
      WRITE(6,1000) (NUMROOT(I), I=1,NIRREP)
 1000 FORMAT(8X, 8I5)
      write(6,*)
      RETURN
      END
