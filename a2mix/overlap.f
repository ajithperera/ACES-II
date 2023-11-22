      SUBROUTINE OVERLAP(L, M, N, GAMA, V)
C
C Calculate the Overlap integral
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DIMENSION DFTR(11)
C
      DFTR(1) = 1.0D0 
      DO I = 2, 11
         DFTR(I) = DFTR(I - 1)*(1.0D0 + 2.0D0*(DFLOAT(I) - 2.0D0))
      ENDDO
      PI = DATAN(1.0D0+00)*4.0D+00
C
      LH = L/2
      MH = M/2
      NH = N/2
C
      IF (2*LH - L) 50, 1, 50
 1    IF (2*MH - M) 50, 2, 50
 2    IF (2*NH - N) 50, 3, 50
C
3     V = (SQRT(PI/GAMA))**3*(0.5/GAMA)**(LH+MH+NH)*DFTR(LH+1)*
     &     DFTR(MH+1)*DFTR(NH+1)
C
      RETURN
C
50    V = 0.D0
C
      RETURN
      END
