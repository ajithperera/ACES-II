      SUBROUTINE DIP_MOMENT(L,M,N,GA,V,NT,D)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DIMENSION FNU(4), V(3), D(3)
C
      FNU(1)=OLAP(L,M,N,GA)
      FNU(2)=OLAP(L+1,M,N,GA)
      FNU(3)=OLAP(L,M+1,N,GA)
      FNU(4)=OLAP(L,M,N+1,GA)
      V(1)=-(FNU(2)+D(1)*FNU(1))
      V(2)=-(FNU(3)+D(2)*FNU(1))
      V(3)=-(FNU(4)+D(3)*FNU(1))
C
      RETURN
      END

