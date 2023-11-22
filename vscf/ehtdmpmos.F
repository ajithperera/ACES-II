      SUBROUTINE EHTDMPMOS(CZAO,EVAL,NOCC,NBFIRR,NBASX,NBAS,NIRREP,
     &                     NATOMS,IUHF)
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
C     Arguments.
C-----------------------------------------------------------------------
      DOUBLE PRECISION CZAO,EVAL
      INTEGER NOCC,NBFIRR,NBASX,NBAS,NIRREP,NATOMS,IUHF
C-----------------------------------------------------------------------
C     Common block variables.
C-----------------------------------------------------------------------
      INTEGER IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C-----------------------------------------------------------------------
C      Local variables.
C-----------------------------------------------------------------------
      DOUBLE PRECISION COORD
      INTEGER ISPIN,IOCC,IBAS,IRREP,I,NCOLS,NFULL,NLEFT,IBLOCK,IROW,J,
     &        IJUNK,ANGTYP,
     &        IATOM
      LOGICAL YESNO
C-----------------------------------------------------------------------
C
      DIMENSION CZAO(NBASX*NBAS),EVAL((IUHF+1)*NBAS)
      DIMENSION NOCC(8,2),NBFIRR(8)
C
#include "maxbasfn.par"
      dimension iocc(maxbasfn,2),coord(3,50),angtyp(maxbasfn),ijunk(50)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
C-----------------------------------------------------------------------
C     Evaluate IOCC (array of occupation numbers).
C-----------------------------------------------------------------------
      DO 30 ISPIN=1,IUHF+1
      IBAS = 0
      DO 20 IRREP=1,NIRREP
      IF(NBFIRR(IRREP).GT.0)THEN
       DO 10 I    =1,NBFIRR(IRREP)
       IBAS = IBAS + 1
       IOCC(IBAS,ISPIN) = 0
       IF(I.LE.NOCC(IRREP,ISPIN)) IOCC(IBAS,ISPIN) = 2 - IUHF
   10  CONTINUE
      ENDIF
   20 CONTINUE
   30 CONTINUE
C
C-----------------------------------------------------------------------
C     Open file EHTMOFILE for writing EHT MOs in AO basis.
C-----------------------------------------------------------------------
      INQUIRE(FILE='EHTMOFILE',EXIST=YESNO)
      IF(YESNO)THEN
       OPEN(71,FILE='EHTMOFILE',STATUS='OLD',ACCESS='SEQUENTIAL',
     &      FORM='FORMATTED')
       CLOSE(UNIT=71,STATUS='DELETE')
      ENDIF
      OPEN(71,FILE='EHTMOFILE',STATUS='NEW',ACCESS='SEQUENTIAL',
     &      FORM='FORMATTED')
C
C-----------------------------------------------------------------------
C     Opening line.
C-----------------------------------------------------------------------
       WRITE(71,       *) NBAS,NBASX,IUHF,(NBFIRR(IRREP),IRREP=1,NIRREP)
C
C-----------------------------------------------------------------------
C     For each spin write MO coefficients, occupation numbers, eigen-
C     values, and angular type array.
C-----------------------------------------------------------------------
      DO 200 ISPIN=1,IUHF+1
C
      NCOLS = 4
      NFULL = NBAS/NCOLS
      NLEFT = NBAS - NCOLS*NFULL
      DO 150 IBLOCK=1,NFULL
      DO 140 IROW  =1,NBASX
      WRITE(71,1020)
     & (CZAO((IBLOCK-1)*NBASX*NCOLS + (J-1)*NBASX + IROW),J=1,NCOLS)
  140 CONTINUE
  150 CONTINUE
C
      IF(NLEFT.NE.0)THEN
       DO 160 IROW  =1,NBASX
       WRITE(71,1020)
     &  (CZAO(     NFULL*NBASX*NCOLS + (J-1)*NBASX + IROW),J=1,NLEFT)
  160  CONTINUE
      ENDIF
C
C     Write occupation numbers and eigenvalues.
C
      DO 170 IBAS=1,NBAS
      WRITE(71,1030) IOCC(IBAS,ISPIN),EVAL((ISPIN-1)*NBAS+IBAS)
  170 CONTINUE
C
C     Write angular type array (length is NBASX). [Fudge here: No
C     use of ANGTYP is anticipated as yet. ANGTYP is used when
C     AOBASMOS (OLDAOMOS) is used. It is only in EHTMOFILE becuase the
C     format is the same as AOBASMOS. Obviously, ANGTYP will have to
C     be properly calculated if a use is anticipated.
C
      DO 180 IBAS=1,NBASX
      ANGTYP(IBAS) = 0
      WRITE(71,*) ANGTYP(IBAS)
  180 CONTINUE
C
  200 CONTINUE
C
C     Write coordinates.
C
      CALL GETREC(20,'JOBARC','COORD',3*NATOMS*IINTFP,COORD)
      DO 250 IATOM=1,NATOMS
      WRITE(71,*) (COORD(I,IATOM),I=1,3)
  250 CONTINUE
C
C     Write computational to ZMAT mapping.
C
      CALL GETREC(20,'JOBARC','MAP2ZMAT',NATOMS,IJUNK)
      DO 260 IATOM=1,NATOMS
      WRITE(71,*) IJUNK(IATOM)
  260 CONTINUE
C
      CLOSE(UNIT=71,STATUS='KEEP')
 1020 FORMAT(4F20.10)
 1030 FORMAT(I3,F20.10)
      RETURN
      END
