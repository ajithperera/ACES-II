C
      SUBROUTINE SYMADAPT(ICORE, MAXCOR, IUHF)
C
C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C & Generates symmetry adapted Fermi-Contact, Spin-dipole &
C & Paramagnetic Spin-Orbit integrals for NMR spin-spin   &
C & coupling constant calculations. Coded by Ajith        &
C & 03/94. Benifitted  from routines written by Jurgen    &
C & Gauss.                                                &
C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C     
      PARAMETER (MAXCENT = 100)
C
      INTEGER DIRPRD
C
      DIMENSION ISYTYP(3), MULNUC(MAXCENT), IPTCNT(6*MAXCENT, 0:7),
     &          IDEGEN(MAXCENT), NUCDEG(MAXCENT), NUCPRE(MAXCENT),
     &          NUCNUM(MAXCENT, 8), COOO(MAXCENT, 3),
     &          CSTRA(6*MAXCENT, 6*MAXCENT), INDEX(6*MAXCENT),
     &          SYMTRN(6*MAXCENT, 6*MAXCENT), IRREPERT(6*MAXCENT),
     &          SYMINVRT(6*MAXCENT, 6*MAXCENT), ICORE(MAXCOR)
C
      LOGICAL JFC, JPSO, JSD, NUCLEI
C
      COMMON /MACHSP/ IINTLN, IFLTLN, IINTFP, IALONE, IBITWD
      COMMON /INFO/ NOCCO(2), NVRTO(2)
      COMMON /SYMINF/ NSTART, NIRREP, IRREPA(255), IRREPB(255), 
     &                DIRPRD(8,8)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /PERT/NTPERT, NPERT(8), IPERT(8)
      COMMON /NMR/JFC, JPSO, JSD, NUCLEI
      COMMON /FILES/ LUOUT, MOINTS
C
      IONE     = 1
      ITHREE   = 3
      IHUNDRD  = 100
      IUNIT    = 30
C
      MXCOR = MAXCOR/IINTFP
C     
C Get the symmetry information from JOBARC file
C
      CALL GETREC (20, 'JOBARC', 'NBASTOT ', IONE,     NAO)
      CALL GETREC (20, 'JOBARC', 'NREALATM', IONE,     NATOM)
      CALL GETREC (20, 'JOBARC', 'MULCINP ', IHUNDRD, MULNUC)
      CALL GETREC (20, 'JOBARC', 'SYMTYPE ', ITHREE,   ISYTYP)
      CALL GETREC (20, 'JOBARC', 'NOSYMMOP', IONE,     NSYMOP)
      CALL GETREC (20, 'JOBARC', 'NOSYATOM', IONE,     NMAX)
C
      IF (IFLAGS(18) .EQ. 8) THEN
        JPSO = .TRUE.
        JFC = .FALSE.
        JSD = .FALSE.
      ENDIF
      IF (IFLAGS(18) .EQ. 9) THEN
        JFC  = .TRUE.
        JPSO = .FALSE.
        JSD = .FALSE.
      ENDIF
      IF (IFLAGS(18) .EQ. 10) THEN
        JSD  = .TRUE.
        JPSO =.FALSE.
        JFC = .FALSE.
      ENDIF
C
      IF (.NOT. JFC .AND. .NOT. JPSO .AND. .NOT. JSD) NUCLEI = .TRUE.
C     
C Make sure MULNUC and ISYTYP are consistent with the rest
C of the program
C
      DO 10 I = 1, NMAX
         MULNUC(I) = MULNUC(I) - 1
 10   CONTINUE
C
      DO 20 I = 1, 3
         ISYTYP(I) = ISYTYP(I) - 1
 20   CONTINUE
C
C Print out the a message about the program, authors
C     
      CALL HEADER('Symmetry adaptation of NMR perturbation integrals',
     &             -1, LUOUT)
C
C Do the symmetry adaptation of nuclear perturbations
C
      CALL DRVSYM (MAXCENT, NMAX, NSYMOP, ISYTYP, MULNUC, IPTCNT,
     &             IDEGEN, NUCDEG, NUCPRE, NUCNUM, COOO, INDEX, 
     &             IRREPERT, CSTRA, SYMTRN, SYMINVRT)
C
      CALL DETPERT (MAXCENT, NMAX, NSYMOP, IPTCNT, CSTRA)
C 
C Open VPOUT and IIII file
C
      OPEN (UNIT=IUNIT, FILE='VPOUT', FORM='UNFORMATTED', STATUS='OLD')
C
C Do the symmetry adaptation of perturbation integrals
C
      IF (JFC) THEN
         CALL DRVFERMI (ICORE, MXCOR, MAXCENT, NATOM, NAO, NMAX, 
     &                  NUCDEG, SYMTRN, IRREPERT, IUHF) 
      ELSE IF (JPSO) THEN
         CALL DRVPSO (ICORE, MXCOR, MAXCENT, NATOM, NAO, NMAX, 
     &                NUCDEG, SYMTRN, IRREPERT, IUHF) 
      ELSE IF (JSD) THEN
         CALL DRVSPNDIP (ICORE, MXCOR, MAXCENT, NATOM, NAO, NMAX, 
     &                   NUCDEG, SYMTRN, IRREPERT, IUHF) 
      ENDIF
C     
C Write several important information into JOBARC file
C
      CALL PUTREC (20, 'JOBARC', 'NTOTPERT', 1, NTPERT)
      CALL PUTREC (20, 'JOBARC', 'NUCIND  ', 1, NMAX)
      CALL PUTREC (20, 'JOBARC', 'NPERTIRR', 8, NPERT)
      CALL PUTREC (20, 'JOBARC', 'IRREPERT', 6*MAXCENT, IRREPERT)
C
      CALL DUMPSYMTRN (ICORE, SYMINVRT, MAXCENT)
C
      CLOSE (IUNIT)
C
      WRITE(LUOUT, *)
C
      RETURN
      END
