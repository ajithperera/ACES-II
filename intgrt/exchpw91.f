c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      SUBROUTINE EXCHPW91(D,S,U,V,EX,VX,LPOT)
C  GGA91 EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
C  INPUT D : DENSITY
C  INPUT S:  ABS(GRAD D)/(2*KF*D)
C  INPUT U:  (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KF)**3)
C  INPUT V: (LAPLACIAN D)/(D*(2*KF)**2)
C  INPUT LPOT: Added by S. Ivanov if LPOT=0, potential is not needed.
C  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter(thrd=0.333333333333D0,thrd4=1.33333333333D0)
      parameter(a1=0.19645D0,a2=0.27430D0,a3=0.15084D0,a4=100.d0)
      parameter(ax=-0.7385588D0,a=7.7956D0,b1=0.004d0)

c     for Becke exchange, set a3=b1=0
c      parameter(a1=0.19645D0,a2=0.27430D0,a3=0.0D0,a4=100.d0)
c      parameter(ax=-0.7385588D0,a=7.7956D0,b1=0.0d0)

      FAC = AX*D**THRD
      S2 = S*S
      S3 = S2*S
      S4 = S2*S2
      P0 = 1.D0/DSQRT(1.D0+A*A*S2)
      P1 = DLOG(A*S+1.D0/P0)
      P2 = DEXP(-A4*S2)
      P3 = 1.D0/(1.D0+A1*S*P1+B1*S4)
      P4 = 1.D0+A1*S*P1+(A2-A3*P2)*S2
      F = P3*P4
      EX = FAC*F
      IF (LPOT.EQ.0) RETURN
C  LOCAL EXCHANGE OPTION
C     EX = FAC
C  ENERGY DONE. NOW THE POTENTIAL:
      P5 = B1*S2-(A2-A3*P2)
      P6 = A1*S*(P1+A*S*P0)
      P7 = 2.D0*(A2-A3*P2)+2.D0*A3*A4*S2*P2-4.D0*B1*S2*F
      FS = P3*(P3*P5*P6+P7)
      P8 = 2.D0*S*(B1-A3*A4*P2)
      P9 = A1*P1+A*A1*S*P0*(3.D0-A*A*S2*P0*P0)
      P10 = 4.D0*A3*A4*S*P2*(2.D0-A4*S2)-8.D0*B1*S*F-4.D0*B1*S3*FS
      P11 = -P3*P3*(A1*P1+A*A1*S*P0+4.D0*B1*S3)
      FSS = P3*P3*(P5*P9+P6*P8)+2.D0*P3*P5*P6*P11+P3*P10+P7*P11
      VX = FAC*(THRD4*F-(U-THRD4*S3)*FSS-V*FS)
C  LOCAL EXCHANGE OPTION:
C     VX = FAC*THRD4
      RETURN
      END
