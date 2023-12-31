      SUBROUTINE CHKSYM(A,N,TOL)
C
C ROUTINE INSPECTS AN NxN MATRIX A AND CHECKS TO SEE
C  IF IT IS SYMMETRIC TO WITHIN A TOLERANCE TOL
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N,N)
      DO 10 IROW=2,N
       DO 20 ICOL=1,IROW-1
        XIJ=A(IROW,ICOL)
        XJI=A(ICOL,IROW)
        Z=ABS(XIJ-XJI)
        IF(ABS(Z).GT.TOL)THEN
         WRITE(6,1000)IROW,ICOL,XIJ,XJI
1000     FORMAT(2I5,2F20.10)
        ENDIF
20     CONTINUE
10    CONTINUE
      RETURN
      END
