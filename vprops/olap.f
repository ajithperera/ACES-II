C
      FUNCTION OLAP(L, M, N, GAMA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /VARS/ FNU(20),DP(20),FN(15),FD(15),PI,DFTR(11)
      LH=L/2
      MH=M/2
      NH=N/2
      IF(2*LH-L) 50,1,50
1     IF(2*MH-M) 50,2,50
2     IF(2*NH-N) 50,3,50
3     OLAP=( SQRT(PI/GAMA))**3*(.5/GAMA)**(LH+MH+NH)*DFTR(LH+1)*DFTR(MH+
     11)*DFTR(NH+1)
      RETURN
50    OLAP=0.
      RETURN
      END
