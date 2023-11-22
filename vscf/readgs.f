C     SUBROUTINE READGS(NOCC,SWAP,LOCK,IPRTGS,ISTOP1,ISTOP2,
C    1                  READMO,WRITMO,IUHFRHF,NIRREP,IUHF)
      SUBROUTINE READGS(NOCC,NIRREP,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*80 TITLE
      INTEGER SWAP,READMO,WRITMO
      LOGICAL YESNO
      LOGICAL GSSOPT,GSSALW,GSSALT,GSSLOK,GSSRED,GSSWRT,GSSUFR
      DIMENSION NOCC(16)
      COMMON /FILES/ LUOUT,MOINTS
C
      COMMON /GSCOMA/ GSSOPT,GSSALW,GSSALT,GSSLOK,GSSRED,GSSWRT,GSSUFR
      COMMON /GSCOMB/ SWAP(4,8,2),LOCK(8,2),IPRTGS(8,2),ISTOP1,ISTOP2,
     1                READMO,WRITMO,IUHFRHF,LUGSS
C
C     GSSOPT - .TRUE. if the GUESS file exists. If it is .TRUE. then the
C              initial guess options are read from the GUESS file.
C     GSSALW - If it is .TRUE. then initial guess parameters are always
C              read from the GUESS file. If it is .FALSE. and GSSOPT is
C              .TRUE. then initial guess parameters are read from GUESS
C              for the first SCF calculation only.
C     GSSALT - If it is .TRUE. then the initial guess orbitals will be
C              swapped according to parameters in the SWAP array.
C     GSSLOK - If it is .TRUE. then attempts will be made in SCFIT to
C              maintain the character of the starting orbitals by monitoring
C              C(OLD)T * S C(NEW). Which symmetry and spin blocks to try
C              to lock are specified by the LOCK array.
C     GSSRED - If it is .TRUE. then the initial orbitals will be read from
C              the formatted file OLDMOS.
C     GSSWRT - If it is .TRUE. then the converged orbitals will be written
C              to the formatted file NEWMOS (in exactly the same format
C              as OLDMOS).
C     GSSUFR - If it is .TRUE. then a UHF guess is generated from a single
C              set of orbitals (the first set on OLDMOS). These may originate
C              from a multitude of sources, including a closed-shell RHF run
C              and the alpha orbitals of a previous  UHF run.
C
      GSSOPT = .FALSE.
      GSSALW = .FALSE.
      GSSALT = .FALSE.
      GSSLOK = .FALSE.
      GSSRED = .FALSE.
      GSSWRT = .FALSE.
      GSSUFR = .FALSE.
C
C     MAXNPR is the maximum number of pairs which may be swapped in
C     a given symmetry block of a given spin.
C
      MAXNPR = 2
      CALL IZERO(SWAP,MAXNPR*2*16)
      CALL IZERO(LOCK,8*2)
      CALL IZERO(IPRTGS,8*2)
      READMO  = 0
      WRITMO  = 0
      ISTOP1  = 0
      ISTOP2  = 0
      IUHFRHF = 0
      IALWAYS = 0
C
      INQUIRE(FILE='GUESS',EXIST=YESNO)
C
      IF(.NOT.YESNO) RETURN
C
C     Since the GUESS file exists, we read the initial guess options
C     therefrom.
C
      GSSOPT = .TRUE.
C
      OPEN(UNIT=LUGSS,FILE='GUESS',STATUS='OLD',FORM='FORMATTED',
     1     ACCESS='SEQUENTIAL')
C
      REWIND LUGSS
C
C     First read the title line.
      READ(LUGSS,'(A)') TITLE
C
C     Read occupations.
      READ(LUGSS,1010) (NOCC(I),I=1,NIRREP)
      IF(IUHF.GT.0)THEN
         READ(LUGSS,1010) (NOCC(I+8),I=1,NIRREP)
      ELSE
         CALL ICOPY(8,NOCC(1),1,NOCC(9),1)
      ENDIF
 1010 FORMAT(8I3)
C
      WRITE(LUOUT,5000)
      WRITE(LUOUT,5001)(NOCC(I),I=1,NIRREP)
      WRITE(LUOUT,5002)(NOCC(8+I),I=1,NIRREP)
      WRITE(LUOUT,5010)
 5000 FORMAT(T3,'@READGS-I, Occupancies from guess file:',/)
 5001 FORMAT(T8,'   Alpha population by irrep: ',8(I3,2X))
 5002 FORMAT(T8,'    Beta population by irrep: ',8(I3,2X))
 5010 FORMAT(/)
C
C     Read orbital swapping information.
      DO   20 ISPIN=1,IUHF+1
      DO   10 IRREP=1,NIRREP
      READ(LUGSS,1020) (SWAP(I,IRREP,ISPIN),I=1,2*MAXNPR)
 1020 FORMAT(4I3)
   10 CONTINUE
   20 CONTINUE
C
      DO   50 ISPIN=1,IUHF+1
      DO   40 IRREP=1,NIRREP
      DO   30 IPAIR=1,MAXNPR
      IF(SWAP(2*IPAIR-1,IRREP,ISPIN).GT.0.AND.
     1   SWAP(2*IPAIR  ,IRREP,ISPIN).GT.0)THEN
      GSSALT = .TRUE.
      ENDIF
   30 CONTINUE
   40 CONTINUE
   50 CONTINUE
C
      IF(GSSALT)THEN
C
      WRITE(LUOUT,1030)
 1030 FORMAT(' @READGS-I, The following orbitals will be swapped: ',/)
      DO   80 ISPIN=1,IUHF+1
      DO   70 IRREP=1,NIRREP
      DO   60 IPAIR=1,MAXNPR
      IF(SWAP(2*IPAIR-1,IRREP,ISPIN).GT.0.AND.
     1   SWAP(2*IPAIR  ,IRREP,ISPIN).GT.0)THEN
      WRITE(LUOUT,1040) ISPIN,IRREP,SWAP(2*IPAIR-1,IRREP,ISPIN),
     1                              SWAP(2*IPAIR  ,IRREP,ISPIN)
 1040 FORMAT('        Spin ',I3,' Symmetry ',I3,' Orbitals ',2I3)
      ENDIF
   60 CONTINUE
   70 CONTINUE
   80 CONTINUE
C
      ELSE
C
      WRITE(LUOUT,1050)
 1050 FORMAT(' @READGS-I, No orbitals will be swapped. ')
      ENDIF
C
C     Read for information about locking guess.
C
      DO  100 ISPIN=1,IUHF+1
      READ(LUGSS,*) (LOCK(I,ISPIN),I=1,NIRREP)
      DO   90 IRREP=1,NIRREP
      IF(LOCK(IRREP,ISPIN).GT.0)THEN
      GSSLOK = .TRUE.
      WRITE(LUOUT,1060) ISPIN,IRREP
 1060 FORMAT(' @READGS-I, For spin ',I3,' symmetry ',I3,
     1       ' occupied orbitals will be locked. ')
      ENDIF
   90 CONTINUE
  100 CONTINUE
C
      IF(.NOT.GSSLOK)THEN
      WRITE(LUOUT,1070)
 1070 FORMAT(' @READGS-I, No symmetry blocks will be locked. ')
      ENDIF
C
C     Read for information about printing guesses.
C
      DO  120 ISPIN=1,IUHF+1
      READ(LUGSS,*) (IPRTGS(I,ISPIN),I=1,NIRREP)
      DO  110 IRREP=1,NIRREP
      IF(IPRTGS(IRREP,ISPIN).GT.0)THEN
      WRITE(LUOUT,1080) ISPIN,IRREP
 1080 FORMAT(' @READGS-I, For spin ',I3,' symmetry ',I3,
     1       ' guess orbitals will be printed. ')
      ENDIF
  110 CONTINUE
  120 CONTINUE
C
C     Read stopping parameters.
C
      READ(LUGSS,*) ISTOP1,ISTOP2
C
C     Read parameters which determine if initial guess is to be read from
C     a file of previously calculated MOs.
C
      READ(LUGSS,*) READMO,WRITMO
C
      IF(READMO.NE.0)THEN
      GSSRED = .TRUE.
      WRITE(LUOUT,2010)
 2010 FORMAT(' @READGS-I, Initial orbitals will be read from OLDMOS. ')
      ENDIF
C
      IF(WRITMO.NE.0)THEN
      GSSWRT = .TRUE.
      WRITE(LUOUT,2020)
 2020 FORMAT(' @READGS-I, Final orbitals will be written to NEWMOS. ')
      ENDIF
C
C     Read parameter which determines if a UHF guess will be generated
C     a single set of orbitals (eg a closed-shell set or the alpha set
C     of a previous UHF solution; the first set of MOs on OLDMOS will
C     be read).
C
      READ(LUGSS,*) IUHFRHF
      IF(IUHFRHF.NE.0)THEN
      GSSUFR = .TRUE.
      WRITE(LUOUT,2030)
 2030 FORMAT(' @READGS-I, UHF guess generated from one set of MOs. ')
      ELSE
      WRITE(LUOUT,2040)
 2040 FORMAT(' @READGS-I, UHF guess generated from two sets of MOs. ')
      ENDIF
C
C     Read parameter which determines if GUESS file is read every time an
C     SCF calculation is performed or just the first time an SCF calculation
C     is done. If the latter, the GUESS file is deleted just before
C     aces_fin is called and some guess options are reset appropriately.
C
      READ(LUGSS,*) IALWAYS
      IF(IALWAYS.NE.0)THEN
      GSSALW = .TRUE.
      WRITE(LUOUT,2050)
 2050 FORMAT(
     1 ' @READGS-I, Guess parameters will always be read from GUESS. ')
      ELSE
      WRITE(LUOUT,2060)
 2060 FORMAT(
     1 ' @READGS-I, Guess parameters will be read from GUESS once. ')
      ENDIF
C
C     We keep the GUESS file for the time being.
C
C     If the GUESS file is always to be read, we keep it. Otherwise, we
C     remove it.
C
      IF(     GSSOPT.AND.GSSALW)THEN
      CLOSE(LUGSS,STATUS='KEEP')
      ELSE
      CLOSE(LUGSS,STATUS='DELETE')
      ENDIF
C
      RETURN
      END
