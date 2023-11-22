
C THIS SUBROUTINE READS A SET OF EIGENVECTORS FROM JOBARC,
C EXPANDS THEM TO THE FULL PRIMITIVE AO BASIS AND THEN
C PROCEEDS TO PROJECT OUT THE TOTALLY SYMMETRIC PART
C OF THE DENSITY MATRIX.
C
C FOR SPHERICAL HARMONIC CASES:
C
C NBASX = THE NUMBER OF COMPUTATIONAL BASIS FUNCTIONS (5d, 7f etc)
C NBASIS= THE ACTUAL NUMBER OF PURE CARTESIAN FUNCTIONS, AS USED
C         IN THE ANALYSIS.

      SUBROUTINE PRJDENS(ICORE,MAXCOR,NBASX,IUHF)
      IMPLICIT INTEGER (A-Z)
C
      LOGICAL YESNO,HANDUP
      DIMENSION ICORE(MAXCOR)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
C

      IONE=1
      CALL GETREC(20,'JOBARC','NAOBASFN',IONE,NBASIS)
      CALL GETREC(20,'JOBARC','NATOMS  ',IONE,NATOMS)
      CALL GETREC(20,'JOBARC','FULLORDR',IONE,IORDGRP)
      I000=1
      I010=I000+NBASIS
      I020=I010+NBASIS
      I030=I020+NATOMS
      I040=I030+NATOMS
      I050=I040+5*NATOMS
      I060=I050+IORDGRP*NATOMS
      IF(MOD(I060-I000+1,2).NE.1)I060=I060+1
      I061=I060+9*IORDGRP*IINTFP
      I062=I061+NBASIS*NBASIS*IINTFP
      I063=I062+NBASIS*NBASIS*IINTFP
      I080=I063+9*IINTFP
C
C  Check first to see if this basis includes f or higher angular momentum
C  functions.  If so, we can't handle it yet, so write out a message
C  and return.
C
      CALL GETREC(20,'JOBARC','ANMOMBF0',NBASIS,ICORE(I010))
      HANDUP=.FALSE.
      DO 100 I=1,NBASIS
        IOFF=I-1
        IF(ICORE(I010+IOFF).GE.5) HANDUP=.TRUE.
  100 CONTINUE
      IF(HANDUP) THEN
        WRITE(LUOUT,5000)
 5000   FORMAT(T3,'@PRJDENS: Basis set contains h and/or higher ',
     &            'angular momentum functions.',/,
     &         T14,'Analysis of density matrix symmetry not ',
     &             'available.',/,
     &         T14,'Check irreps of the molecular orbitals in the ',
     &             'full point group',/,
     &         T14,'for presence of an asymmetric solution.',/)
        RETURN
      ENDIF
C
C  Guess there weren't any f or higher functions; let's do the analysis.
C
      CALL SETUPPRJ(NATOMS,NBASIS,IORDGRP,ICORE(I000),ICORE(I010),
     &              ICORE(I020),ICORE(I030),ICORE(I040),ICORE(I050),
     &              ICORE(I060),ICORE(I061),ICORE(I080))
      ICORR=0
cYAU      WRITE(LUOUT,5100)
cYAU 5100 FORMAT(T3,'@PRJDENS: Analyzing reference function density.')
      CALL PROJECT(NATOMS,NBASIS,IORDGRP,ICORE(I000),ICORE(I010),
     &            ICORE(I020),ICORE(I030),ICORE(I040),ICORE(I050),
     &            ICORE(I060),ICORE(I063),ICORE(I061),ICORE(I062),
     &            IUHF,ICORR,NBASX)
C

      RETURN
      END

