      SUBROUTINE SEC_MOMENT(L,M,N,GA,V,NT,D)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION FNU(9),V(6),D(3)
      NT=6
      FNU(1)=OLAP(L,M,N,GA)
      FNU(2)=OLAP(L+1,M,N,GA)
      FNU(3)=OLAP(L,M+1,N,GA)
      FNU(4)=OLAP(L,M,N+1,GA)
      FNU(5)=OLAP(L+2,M,N,GA)
      FNU(6)=OLAP(L,M+2,N,GA)
      FNU(7)=OLAP(L,M,N+2,GA)
      FNU(8)=OLAP(L+1,M+1,N,GA)
      FNU(9)=OLAP(L+1,M,N+1,GA)
      FNU(10)=OLAP(L,M+1,N+1,GA)
      V(1)=FNU(5)+2.*D(1)*FNU(2)+D(1)**2*FNU(1)
      V(2)=FNU(6)+2.*D(2)*FNU(3)+D(2)**2*FNU(1)
      V(3)=FNU(7)+2.*D(3)*FNU(4)+D(3)**2*FNU(1)
      V(4)=FNU(8)+D(1)*FNU(3)+D(2)*FNU(2)+D(1)*D(2)*FNU(1)
      V(5)=FNU(9)+D(1)*FNU(4)+D(3)*FNU(2)+D(1)*D(3)*FNU(1)
      V(6)=FNU(10)+D(2)*FNU(4)+D(3)*FNU(3)+D(2)*D(3)*FNU(1)

      RETURN
      END

