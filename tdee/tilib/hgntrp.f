      SUBROUTINE HGNTRP (TOUT, Y, N0, Y0)
C
C
C     Description of routine.
C
C     Subroutine HGNTRP computes interpolated values of the dependent
C     variable Y and stores them in Y0.  the interpolation is to the
C     point T = TOUT, and uses the Nordsieck history array Y, as follows
C     
C                NQ
C     Y0(I)  =  sum  Y(I,J+1)*S**J ,
C               J=0
C     
C     where S = -(T-TOUT)/H.
C
C     Original code: copyright A. C. Hindmarsh and LLL, 1975.
C     Modifications: copyright Erik Deumens and QTP, 1986, 1992.
C
      IMPLICIT NONE

      INTEGER I, L, J, N0
      DOUBLE PRECISION Y0(N0),Y(N0,1),T,H,HMIN,HMAX,EPSC,UROUND
      INTEGER N, MF,KFLAG, JSTART
      COMMON /GEAR1/ T,H,HMIN,HMAX,EPSC,UROUND,N,MF,KFLAG,JSTART
      DOUBLE PRECISION S, S1, TOUT
C-----------------------------------------------------------------------
      DO 10 I = 1,N
 10     Y0(I) = Y(I,1)
      L = JSTART + 1
      S = (TOUT - T)/H
      S1 = 1.D0
      DO 30 J = 2,L
        S1 = S1*S
        DO 20 I = 1,N
 20       Y0(I) = Y0(I) + S1*Y(I,J)
 30     CONTINUE
      RETURN
      END
