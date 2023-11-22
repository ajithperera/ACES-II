      SUBROUTINE EVALREPLP(L,M,N,GA,V,D)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PRTBUG/ IPRTFS
C
      DIMENSION V(10),D(10), VS(10)
      DIMENSION DD(10),DS(10)
C
      COMMON /VARS/ FNU(20),DP(20),FN(15),FD(15),PI,DFTR(11)
C
      DATA FQR/1.3333333333333D0/
      DATA SQRPI/1.772453851D0/
C
C Set the factorial arrays
C      
      PI = DATAN(1.0D0+00)*4.0D+00
      FN(1) = 1.0D0
      FD(1) = 1.0D0
      DO I = 2, 15
         FN(I) = (I - 1)*FN(I - 1)
         FD(I) = 1.0D0/FN(I)
      ENDDO
C
      NT=1
      ITOT=L+M+N
      ITOTH=ITOT/2
      PRE=(2.*PI/GA)*FN(L+1)*FN(M+1)*FN(N+1)
C
      IF (2*ITOTH - ITOT) 1, 2, 1
C
 1    PRE=-PRE
 2    DEL=.25/GA
      X=GA*(D(1)**2+D(2)**2+D(3)**2)
      XX = X + X
C
      CALL FMC(12,X,EXPMX,FMCH)
C
      FNU(13)=FMCH
      FNU(12)=(EXPMX+XX*FNU(13))/23.
      FNU(11)=(EXPMX+XX*FNU(12))/21.
      FNU(10)=(EXPMX+XX*FNU(11))/19.
      FNU(9)=(EXPMX+XX*FNU(10))/17.
      FNU(8)=(EXPMX+XX*FNU(9))/15.
      FNU(7)=(EXPMX+XX*FNU(8))/13.
      FNU(6)=(EXPMX+XX*FNU(7))/11.
      FNU(5)=(EXPMX+XX*FNU(6))/9.
      FNU(4)=(EXPMX+XX*FNU(5))/7.
      FNU(3)=(EXPMX+XX*FNU(4))/5.
      FNU(2)=(EXPMX+XX*FNU(3))/3.
      FNU(1)=EXPMX+XX*FNU(2)
C
      DP(1)=1.0D0
      DP(2)=DEL
      DP(3)=DEL**2
      DP(4)=DP(3)*DEL
      DP(5)=DP(4)*DEL
      DP(6)=DP(5)*DEL
      DP(7)=DP(6)*DEL
      V(1) =PRE*AAINER(0,0,0,L,M,N,D)
C$$      WRITE(6,*) V(1)
C

      RETURN
      END
