      SUBROUTINE MKBASMAP(IMAP,FACTOR,IZATOM,
     &                    NSSHELL,NPSHELL,NDSHELL,
     &                    NSTOT,NPTOT,NDTOT,TYPE)
      IMPLICIT NONE
C-----------------------------------------------------------------------
      CHARACTER*6 TYPE
      DOUBLE PRECISION FACTOR
      INTEGER IMAP,IZATOM,NSTOT,NPTOT,NDTOT
C-----------------------------------------------------------------------
      INTEGER I,IOFF,ISHELL,NS,NP,NSSHELL,NPSHELL,NDSHELL
      LOGICAL FIRSTROW,SECNDROW,THIRDROW
C-----------------------------------------------------------------------
      DIMENSION IMAP(1),FACTOR(1)
      DIMENSION NS(7),NP(6)
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C     Set NSSHELL, NPSHELL, NDSHELL based on atomic number.
C
C     NSSHELL --- Number of s shells for this atom. The value is 1 for
C                 H,He; it is 2 for Li-Ne; it is 3 for Na-Ar.
C     NPSHELL --- Number of p shells for this atom. The value is 0 for
C                 H,He; it is 1 for Li-Ne; it is 2 for Na-Ar.
C     NDSHELL --- 0 for the moment.
C
C     Note that these variables can serve as the numbers of functions of
C     different types in a minimal basis set.
C-----------------------------------------------------------------------
C
      FIRSTROW = .FALSE.
      SECNDROW = .FALSE.
      THIRDROW = .FALSE.
C
      IF(IZATOM .EQ.  1 .OR. IZATOM .EQ. 2)THEN
       NSSHELL = 1
       NPSHELL = 0
       NDSHELL = 0
       FIRSTROW = .TRUE.
      ENDIF
C
      IF(IZATOM .GE.  3 .AND. IZATOM .LE. 10)THEN
       NSSHELL = 2
       NPSHELL = 1
       NDSHELL = 0
       SECNDROW = .TRUE.
      ENDIF
C
      IF(IZATOM .GE. 11 .AND. IZATOM .LE. 18)THEN
       NSSHELL = 3
       NPSHELL = 2
       NDSHELL = 0
       THIRDROW = .TRUE.
      ENDIF
C
C-----------------------------------------------------------------------
C     Set NS, NP based on type of basis set.
C
C     NS      --- The number of s functions in the basis set represent-
C                 ing the present atomic shell.
C     NP      --- The number of p functions in the basis set represent-
C                 ing the present atomic shell.
C-----------------------------------------------------------------------
C
      IF(TYPE .EQ. 'SZCDZV')THEN
       IF(FIRSTROW) NS(1) = 2
       IF(SECNDROW) NS(1) = 1
       IF(SECNDROW) NS(2) = 2
       IF(THIRDROW) NS(1) = 1
       IF(THIRDROW) NS(2) = 1
       IF(THIRDROW) NS(3) = 2
C
       IF(SECNDROW) NP(1) = 2
       IF(THIRDROW) NP(1) = 1
       IF(THIRDROW) NP(2) = 2
      ENDIF
C
      IF(TYPE .EQ. 'DZ    ')THEN
       IF(FIRSTROW) NS(1) = 2
       IF(SECNDROW) NS(1) = 2
       IF(SECNDROW) NS(2) = 2
       IF(THIRDROW) NS(1) = 2
       IF(THIRDROW) NS(2) = 2
       IF(THIRDROW) NS(3) = 2
C
       IF(SECNDROW) NP(1) = 2
       IF(THIRDROW) NP(1) = 2
       IF(THIRDROW) NP(2) = 2
      ENDIF
C
      IF(TYPE .EQ. 'SZCTZV')THEN
       IF(FIRSTROW) NS(1) = 3
       IF(SECNDROW) NS(1) = 1
       IF(SECNDROW) NS(2) = 3
       IF(THIRDROW) NS(1) = 1
       IF(THIRDROW) NS(2) = 1
       IF(THIRDROW) NS(3) = 3
C
       IF(SECNDROW) NP(1) = 3
       IF(THIRDROW) NP(1) = 1
       IF(THIRDROW) NP(2) = 3
      ENDIF
C-----------------------------------------------------------------------
C
      IOFF = 0
      DO 20 ISHELL=1,NSSHELL
C
      DO 10 I=1,NS(ISHELL)
      IMAP(IOFF+I) = ISHELL
   10 CONTINUE
C
      IF(NS(ISHELL) .EQ. 1)THEN
       FACTOR(IOFF+1) = 1.0D+00
      ELSEIF(NS(ISHELL) .EQ. 2)THEN
       FACTOR(IOFF+1) =        1.0D+00 / DSQRT(2.0D+00)
       FACTOR(IOFF+2) =        1.0D+00 / DSQRT(2.0D+00)
      ELSEIF(NS(ISHELL) .EQ. 3)THEN
       FACTOR(IOFF+1) =        1.0D+00 / DSQRT(3.0D+00)
       FACTOR(IOFF+2) =        1.0D+00 / DSQRT(3.0D+00)
       FACTOR(IOFF+3) =        1.0D+00 / DSQRT(3.0D+00)
      ELSEIF(NS(ISHELL) .EQ. 4)THEN
       FACTOR(IOFF+1) = 0.25D+00
       FACTOR(IOFF+2) = 0.25D+00
       FACTOR(IOFF+3) = 0.25D+00
       FACTOR(IOFF+4) = 0.25D+00
      ENDIF
C
      IOFF = IOFF + NS(ISHELL)
   20 CONTINUE
C
      IF(IOFF .NE. NSTOT)THEN
       write(6,*) ' @MKBASMAP-F, Offset mismatch. ',IOFF,NSTOT
       call errex
      ENDIF
C
      IF(NPSHELL .GE. 1)THEN
C
       DO 40 ISHELL=1,NPSHELL
C
       DO 30 I=1,NP(ISHELL)
       IMAP(IOFF+(I-1)*3 + 1) = NSSHELL + (ISHELL-1)*3 + 1
       IMAP(IOFF+(I-1)*3 + 2) = NSSHELL + (ISHELL-1)*3 + 2
       IMAP(IOFF+(I-1)*3 + 3) = NSSHELL + (ISHELL-1)*3 + 3
   30  CONTINUE
C
       IF(NP(ISHELL) .EQ. 1)THEN
        FACTOR(IOFF+1) = 1.0D+00
        FACTOR(IOFF+2) = 1.0D+00
        FACTOR(IOFF+3) = 1.0D+00
       ELSEIF(NP(ISHELL) .EQ. 2)THEN
        FACTOR(IOFF+1) =       1.0D+00 /DSQRT(2.0D+00)
        FACTOR(IOFF+2) =       1.0D+00 /DSQRT(2.0D+00)
        FACTOR(IOFF+3) =       1.0D+00 /DSQRT(2.0D+00)
        FACTOR(IOFF+4) =       1.0D+00 /DSQRT(2.0D+00)
        FACTOR(IOFF+5) =       1.0D+00 /DSQRT(2.0D+00)
        FACTOR(IOFF+6) =       1.0D+00 /DSQRT(2.0D+00)
       ELSEIF(NP(ISHELL) .EQ. 3)THEN
        FACTOR(IOFF+1) =       1.0D+00 /DSQRT(3.0D+00)
        FACTOR(IOFF+2) =       1.0D+00 /DSQRT(3.0D+00)
        FACTOR(IOFF+3) =       1.0D+00 /DSQRT(3.0D+00)
        FACTOR(IOFF+4) =       1.0D+00 /DSQRT(3.0D+00)
        FACTOR(IOFF+5) =       1.0D+00 /DSQRT(3.0D+00)
        FACTOR(IOFF+6) =       1.0D+00 /DSQRT(3.0D+00)
        FACTOR(IOFF+7) =       1.0D+00 /DSQRT(3.0D+00)
        FACTOR(IOFF+8) =       1.0D+00 /DSQRT(3.0D+00)
        FACTOR(IOFF+9) =       1.0D+00 /DSQRT(3.0D+00)
       ELSEIF(NP(ISHELL) .EQ. 4)THEN
        FACTOR(IOFF+1) = 0.25D+00
        FACTOR(IOFF+2) = 0.25D+00
        FACTOR(IOFF+3) = 0.25D+00
        FACTOR(IOFF+4) = 0.25D+00
        FACTOR(IOFF+5) = 0.25D+00
        FACTOR(IOFF+6) = 0.25D+00
        FACTOR(IOFF+7) = 0.25D+00
        FACTOR(IOFF+8) = 0.25D+00
        FACTOR(IOFF+9) = 0.25D+00
        FACTOR(IOFF+10) = 0.25D+00
        FACTOR(IOFF+11) = 0.25D+00
        FACTOR(IOFF+12) = 0.25D+00
       ENDIF
C
       IOFF = IOFF + 3*NP(ISHELL)
   40  CONTINUE
C
       IF(IOFF .NE. NSTOT + 3*NPTOT)THEN
        write(6,*) ' @MKBASMAP-F, Offset mismatch. 40 ',IOFF,NSTOT,NPTOT
        call errex
       ENDIF
C
      ENDIF
C
      IF(IZATOM .EQ. 1 .OR. IZATOM .EQ. 2)THEN
       IF(NPTOT .GT. 0)THEN
        DO 50 I=1,NPTOT
        IMAP(IOFF+I  ) = 999
        IMAP(IOFF+I+1) = 999
        IMAP(IOFF+I+2) = 999
C
        FACTOR(IOFF+I  ) = 0.0D+00
        FACTOR(IOFF+I+1) = 0.0D+00
        FACTOR(IOFF+I+2) = 0.0D+00
C
        IOFF = IOFF + 3
C
   50   CONTINUE
       ENDIF
      ENDIF
C
      IF(NDTOT .GT. 0)THEN
       DO 60 I     =1,NDTOT
       IMAP(IOFF+I  ) = 999
       IMAP(IOFF+I+1) = 999
       IMAP(IOFF+I+2) = 999
       IMAP(IOFF+I+3) = 999
       IMAP(IOFF+I+4) = 999
       IMAP(IOFF+I+5) = 999
C
       FACTOR(IOFF+I  ) = 0.0D+00
       FACTOR(IOFF+I+1) = 0.0D+00
       FACTOR(IOFF+I+2) = 0.0D+00
       FACTOR(IOFF+I+3) = 0.0D+00
       FACTOR(IOFF+I+4) = 0.0D+00
       FACTOR(IOFF+I+5) = 0.0D+00
C
       IOFF = IOFF + 6
C
   60  CONTINUE
      ENDIF
C
      RETURN
      END
