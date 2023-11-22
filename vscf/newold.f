      SUBROUTINE NEWOLD(A,N,EVEC,SCR,NOCC,IRREP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL ROHF
      DIMENSION A(N,N),EVEC(N,N),SCR(N)
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /POPUL/ NOC(8,2)
C
C     Subroutine to examine the matrix
C
C               T
C     A = C(OLD)  * S * C(NEW)
C
C     and reorder the eigenvectors C(NEW) so that they look as much like
C     C(OLD) as possible. Strategy focusses on occupied part of C(OLD).
C     By examining rows of A, we search for the members of C(NEW) which
C     look most like the occupied orbitals of C(OLD).
C
      ROHF = .FALSE.
      ROHF = IFLAGS(11).EQ.2
C
      IF(NOCC.LE.0) RETURN
C
C     For UHF exchanges among the occupied orbitals are of no signif-
C     icance. Only exchanges between the occupied and virtual orbitals
C     are significant.
C
C     For ROHF exchanges among just the doubly occupied orbitals and
C     among just the partially occupied orbitals are of no significance.
C     D-P, D-V, P-V (D=doubly occupied, P=partially occupied, and V=
C     virtual orbitals) exchanges are significant.
C
C     Examine contributions of old occupied orbitals to new ones.
C
      IF (IFLAGS(1) .GE. 10) THEN
         DO   20 J=1,NOCC
            DO   10 I=1,NOCC
               WRITE(6,1010) I,J,A(I,J)
 1010 FORMAT(' NEWOLD-I, Contribution of old I to new J ',2I4,F12.7)
   10       CONTINUE
   20    CONTINUE
      ENDIF
C
C     RHF/UHF
C
      IF(.NOT.ROHF)THEN
C
      DO   50 I=1,NOCC
      AMAXI = 0.0D+00
      DO   30 J=1,NOCC
      IF(DABS(A(I,J)).GT.AMAXI)THEN
      AMAXI = DABS(A(I,J))
      JMAX  = J
      ENDIF
   30 CONTINUE
C
C     If no new occupied orbital looks very much like old orbital I,
C     we must search the virtual orbitals too.
C
ctest      IF(AMAXI.LT.0.2D+00)THEN
C
c      WRITE(6,1020) I
 1020 FORMAT(' NEWOLD-I, Old orbital ',I3,' does not contribute. ')
C
      AMAXI = 0.0D+00
      JMAX  = 0
      DO   40 J=1,N
      IF(DABS(A(I,J)).GT.AMAXI)THEN
      AMAXI = DABS(A(I,J))
      JMAX  = J
      ENDIF
   40 CONTINUE
      WRITE(6,1030) I,JMAX
 1030 FORMAT(' @NEWOLD-I, ',I3,' contributes most to ',I3)
      CALL SWAPORB(EVEC,SCR,N,I,JMAX)
      CALL SWAPORB(   A,SCR,N,I,JMAX)
C
ctest      ENDIF
   50 CONTINUE
C
      ELSE
C
C     ROHF
C
C     Straightforward, unsophisticated algorithm. No attempt to prevent
C     redundant D-D, P-P mixing at this stage.
C
Ctest --- investigate possibility of just locking partially occupied orbitals.
Ctest --- sloppy : we assume number of alpha > number of beta, also assume
C                  all shells have some open (FeCl4 test). clean up later.
C
      DO  100 I=1,NOCC
c      DO  100 I=NOC(IRREP,2)+1,NOC(IRREP,1)
      AMAXI = 0.0D+00
      DO   80 J=1,NOCC
      IF(DABS(A(I,J)).GT.AMAXI)THEN
      AMAXI = DABS(A(I,J))
      JMAX  = J
      ENDIF
   80 CONTINUE
C
      IF(AMAXI.LT.0.2D+00)THEN
      WRITE(LUOUT,1020) I
      ENDIF
C
      AMAXI = 0.0D+00
      JMAX  = 0
      DO   90 J=1,N
      IF(DABS(A(I,J)).GT.AMAXI)THEN
      AMAXI = DABS(A(I,J))
      JMAX  = J
      ENDIF
   90 CONTINUE
      WRITE(LUOUT,1030) I,JMAX
      CALL SWAPORB(EVEC,SCR,N,I,JMAX)
      CALL SWAPORB(   A,SCR,N,I,JMAX)
  100 CONTINUE
C
      ENDIF
C
      RETURN
      END
