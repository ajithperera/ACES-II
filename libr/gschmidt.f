      SUBROUTINE GSCHMIDT(VEC,VORTH,NSIZE,NDIM,TMP,RESID,TOL)
C
C THIS PROJECTS OUT ALL PARTS OF AN INPUT VECTOR (VEC)
C WHICH LIE IN THE SPACE SPANNED BY THE ORTHOGONAL BASIS
C VORTH.
C
C   |v'> = |v> - SUM <i|v> |i>
C                 i 
C
C WHERE THE |i> ARE NORMALIZED BASIS VECTORS FOR THE SPACE VORTH
C
C INPUT:
C       VEC : THE VECTOR WHICH IS TO BE ORTHOGONALIZED TO
C             THE EXISTING BASIS.  *IT IS ASSUMED THAT VEC
C             IS NORMALIZED ON INPUT)
C     VORTH : THE BASIS VECTORS FOR THE EXISTING ORTHOGONAL
C             BASIS
C     NSIZE : THE LENGTH OF THE BASIS VECTORS
C      NDIM : THE DIMENSION OF THE ORTHOGONAL SPACE
C       TMP : A SCRATCH VECTOR OF LENGTH NSIZE
C     RESID : THE NORM OF VEC, AFTER ORTHONALIZATION AND
C             BEFORE NORMALIZATION
C       TOL : TOLERANCE FOR RENORMALIZATION
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VEC(NSIZE),VORTH(NSIZE,NDIM),TMP(NSIZE)
C
      DATA ONE /1.0/
C
      IF(NDIM.EQ.0)THEN
       RESID=SNRM2(NSIZE,VEC,1)
       RETURN
      ENDIF
C
      CALL SCOPY(NSIZE,VEC,1,TMP,1)
      DO 10 I=1,NDIM
       FACT=SDOT(NSIZE,VORTH(1,I),1,VEC,1)
       CALL SAXPY(NSIZE,-FACT,VORTH(1,I),1,TMP,1)
10    CONTINUE
      CALL SCOPY(NSIZE,TMP,1,VEC,1)
C
C RENORMALIZE THE RESIDUAL
C
      RESID=SNRM2(NSIZE,VEC,1)
      IF(RESID.GT.TOL)THEN
       X=ONE/RESID
       CALL SSCAL(NSIZE,X,VEC,1)
      ENDIF 
      RETURN
      END
