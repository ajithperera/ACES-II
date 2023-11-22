      SUBROUTINE EHTDIM(NBAS,NBASX)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C     Arguments.
C-----------------------------------------------------------------------
      INTEGER NBAS,NBASX
C-----------------------------------------------------------------------
C     Parameter.
C-----------------------------------------------------------------------
      INTEGER MXNATOM
C-----------------------------------------------------------------------
C     Local variables.
C-----------------------------------------------------------------------
      INTEGER NATOMS,IZ,IATOM
C-----------------------------------------------------------------------
C     Common block variables.
C-----------------------------------------------------------------------
      INTEGER IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C-----------------------------------------------------------------------
      PARAMETER(MXNATOM=50)
C-----------------------------------------------------------------------
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C-----------------------------------------------------------------------
      DIMENSION IZ(MXNATOM)
C-----------------------------------------------------------------------
C
      CALL GETREC(20,'JOBARC','NATOMS',1,NATOMS)
      CALL GETREC(20,'JOBARC','ATOMCHRG',NATOMS,IZ)
C
      NBAS = 0
      DO  50 IATOM=1,NATOMS
C
      IF(IZ(IATOM).GE. 1 .AND. IZ(IATOM).LE. 2) NBAS = NBAS + 1
      IF(IZ(IATOM).GE. 3 .AND. IZ(IATOM).LE.10) NBAS = NBAS + 5
      IF(IZ(IATOM).GE.11 .AND. IZ(IATOM).LE.18) NBAS = NBAS + 9
C
      IF(IZ(IATOM).GT.18)THEN
       write(6,*) ' @EHTDIM-F, Charge cannot be greater than 18 '
       call errex
      ENDIF
C
   50 CONTINUE
C
      NBASX = NBAS
C
      RETURN
      END
