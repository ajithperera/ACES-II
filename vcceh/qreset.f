C
      SUBROUTINE QRESET(ICORE, MAXCOR, IUHF)
C
C This routine resets the spin adapted ring lists. On the entry
C list 56 and 58 have the following structure.
C
C  2*Hbar(AB,AB) - Hbar(AB,BA) on list 56
C  Hbar(AB,BA) on list 58
C
C The purpose of this routine is to write the following lists
C
C  Hbar(AB,AB) on list 56
C  Hbar(AB,BA) on list 58
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,HALF,ONEM
      DIMENSION ICORE(MAXCOR)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      DATA HALF /0.5/
      DATA ONEM /-1.0/
      DATA ONE /1.0/
C
      ISIZE = ISYMSZ(ISYTYP(1, 56), ISYTYP(2, 56))
C
      I000 = 1
      I010 = I000 + IINTFP*ISIZE
      I020 = I010 + IINTFP*ISIZE
C      
      CALL GETALL(ICORE(I000), ISIZE, 1, 56)
      CALL GETALL(ICORE(I010), ISIZE, 1, 58)
C
      CALL SAXPY (ISIZE, ONE, ICORE(I010), 1, ICORE(I000), 1)
      CALL SSCAL (ISIZE, HALF, ICORE(I000), 1)
      CALL PUTALL(ICORE(I000), ISIZE, 1, 56)
C
      CALL SAXPY (ISIZE, ONEM, ICORE(I000), 1, ICORE(I010), 1)
      CALL PUTALL(ICORE(I010), ISIZE, 1, 54)
C
      RETURN
      END
