

      SUBROUTINE GSO2AO(ICORE,MAXCOR,NBAS,NAMIN,NAMOUT)
C
C
C  This routine takes any square nbas*nbas quantity and takes it
C  from the SO representation to the pure AO representation.  It is
C  required mainly with generating data for use with the CAChe
C  workstation, at least for now.
C
C      NAMIN = Name of record on JOBARC containing the SO quantity.
C
C     NAMOUT = Name of record on JOBARC to place the AO quantity.
C
CEND
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*8 NAMIN,NAMOUT
      DIMENSION ICORE(MAXCOR)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA ONE /1.0/
      DATA ZILCH /0.0/
C
      I000=1
      I010=I000+NBAS*NBAS*IINTFP
      I020=I010+NBAS*NBAS*IINTFP
      I030=I020+NBAS*NBAS*IINTFP
C
      CALL GETREC(20,'JOBARC','AO2SO   ',NBAS*NBAS*IINTFP,ICORE(I000))
      CALL GETREC(20,'JOBARC',NAMIN,NBAS*NBAS*IINTFP,ICORE(I010))
C
      CALL XGEMM('N','N',NBAS,NBAS,NBAS,ONE,ICORE(I000),NBAS,
     &           ICORE(I010),NBAS,ZILCH,ICORE(I020),NBAS)
      CALL PUTREC(20,'JOBARC',NAMOUT,NBAS*NBAS*IINTFP,ICORE(I020))
C
      RETURN
      END
