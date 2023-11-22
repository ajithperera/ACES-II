      SUBROUTINE ROTMAT (J,I,JORBS,IORBS,NUMAT,COORD,R,YY)
C     *
C     ROTATION MATRIX FOR A GIVEN ATOM PAIR I-J (I.GT.J).
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     J,I       NUMBERS OF ATOMS IN PAIR I-J (I).
C     JORBS     NUMBER OF BASIS FUNCTIONS AT ATOM J (I).
C     IORBS     NUMBER OF BASIS FUNCTIONS AT ATOM I (I).
C     NUMAT     NUMBER OF ATOMS IN ARRAY COORD (I).
C     COORD()   CARTESIAN COORDINATES, IN ANGSTROM (I).
C     R         INTERATOMIC DISTANCE, IN ATOMIC UNITS (O).
C     YY()      PRECOMBINED ELEMENTS OF THE ROTATION MATRIX (O).
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (SMALL=1.0D-07)
      PARAMETER (PT5SQ3=0.8660254037841D0)
c      COMMON
c     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
c     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
      DATA A0    /  0.529167D0/
      DATA AFACT / 57.29577951308232D0/
      DATA EV    / 27.21D0/
      DATA EVCAL / 23.061D0/
      DATA PI    /  3.141592653589793D0/
      DATA W1    / 14.399D0/
      DATA W2    /  7.1995D0/
      DATA BIGEXP/ 50.0D0/
C *** NUMERICAL CONSTANTS.
      DATA ZERO  /  0.0D0/
      DATA ONE   /  1.0D0/
      DATA TWO   /  2.0D0/
      DATA THREE /  3.0D0/
      DATA FOUR  /  4.0D0/
      DATA PT5   /  0.5D0/
      DATA PT25  /  0.25D0/
      DIMENSION COORD(3,NUMAT),YY(15,45)
      DIMENSION P(3,3),D(5,5)
      DIMENSION INDX(9)
      DATA INDX /0,1,3,6,10,15,21,28,36/
C *** CALCULATE GEOMETRIC DATA AND INTERATOMIC DISTANCE.
C     CA  = COS(PHI)    , SA  = SIN(PHI)
C     CB  = COS(THETA)  , SB  = SIN(THETA)
C     C2A = COS(2*PHI)  , S2A = SIN(2*PHI)
C     C2B = COS(2*THETA), S2B = SIN(2*PHI)
      X11    = COORD(1,J)-COORD(1,I)
      X22    = COORD(2,J)-COORD(2,I)
      X33    = COORD(3,J)-COORD(3,I)
      B      = X11*X11+X22*X22
      R      = SQRT(B+X33*X33)
      SQB    = SQRT(B)
      SB     = SQB/R
C     CHECK FOR SPECIAL CASE (BOTH ATOMS ON Z AXIS).
      IF(SB.GT.SMALL) THEN
         CA  = X11/SQB
         SA  = X22/SQB
         CB  = X33/R
      ELSE
         SA  = ZERO
         SB  = ZERO
         IF(X33.LT.ZERO) THEN
            CA  =-ONE
            CB  =-ONE
         ELSE IF(X33.GT.ZERO) THEN
            CA  = ONE
            CB  = ONE
         ELSE
            CA  = ZERO
            CB  = ZERO
         ENDIF
      ENDIF
C     CONVERT DISTANCE TO ATOMIC UNITS.
      R      = R/A0
C *** CALCULATE ROTATION MATRIX ELEMENTS.
      P(1,1) = CA*SB
      P(2,1) = CA*CB
      P(3,1) =-SA
      P(1,2) = SA*SB
      P(2,2) = SA*CB
      P(3,2) = CA
      P(1,3) = CB
      P(2,3) =-SB
      P(3,3) = ZERO
      IF(IORBS.GE.9 .OR. JORBS.GE.9) THEN
         C2A    = TWO*CA*CA-ONE  
         C2B    = TWO*CB*CB-ONE  
         S2A    = TWO*SA*CA
         S2B    = TWO*SB*CB
         D(1,1) = PT5SQ3*C2A*SB*SB
         D(2,1) = PT5*C2A*S2B
         D(3,1) =-S2A*SB
         D(4,1) = C2A*(CB*CB+PT5*SB*SB)
         D(5,1) =-S2A*CB
         D(1,2) = PT5SQ3*CA*S2B
         D(2,2) = CA*C2B
         D(3,2) =-SA*CB
         D(4,2) =-PT5*CA*S2B
         D(5,2) = SA*SB
         D(1,3) = CB*CB-PT5*SB*SB
         D(2,3) =-PT5SQ3*S2B
         D(3,3) = ZERO
         D(4,3) = PT5SQ3*SB*SB
         D(5,3) = ZERO
         D(1,4) = PT5SQ3*SA*S2B
         D(2,4) = SA*C2B
         D(3,4) = CA*CB
         D(4,4) =-PT5*SA*S2B
         D(5,4) =-CA*SB
         D(1,5) = PT5SQ3*S2A*SB*SB
         D(2,5) = PT5*S2A*S2B
         D(3,5) = C2A*SB
         D(4,5) = S2A*(CB*CB+PT5*SB*SB)
         D(5,5) = C2A*CB
      ENDIF
C     *
C     PRECOMBINE ROTATION MATRIX ELEMENTS.
C     *
C     THE FIRST INDEX OF YY(IJ,KL) IS CONSECUTIVE. AS MANY ELEMENTS
C     AS NEEDED ARE DEFINED (1 FOR KL=SS, 3 FOR KL=PS, 6 FOR KL=PP,
C     5 FOR KL=DS, 15 FOR KL=DP, AND 15 FOR KL=DD).
C     THE SECOND INDEX OF YY(IJ,KL) IS A STANDARD PAIR INDEX.
C     KL=(K*(K-1)/2+L , ORDER OF K AND L AS IN INTEGRAL EVALUATION.
C     *

c      print*,'indx',indx

C     S-S
      YY(1,1)   = ONE
C     P-S
      DO 10 K=1,3
      KL        = INDX(K+1)+1
c      print*,'first set',kl
      YY(1,KL)  = P(K,1)
      YY(2,KL)  = P(K,2)
   10 YY(3,KL)  = P(K,3)
C     P-P
      DO 20 K=1,3
      KL        = INDX(K+1)+K+1
c      print*,'second set',kl
      YY(1,KL)  = P(K,1)*P(K,1)
      YY(2,KL)  = P(K,1)*P(K,2)
      YY(3,KL)  = P(K,2)*P(K,2)
      YY(4,KL)  = P(K,1)*P(K,3)
      YY(5,KL)  = P(K,2)*P(K,3)
   20 YY(6,KL)  = P(K,3)*P(K,3)
      DO 30 K=2,3
      DO 30 L=1,K-1
      KL        = INDX(K+1)+L+1
c      print*,'third set',kl
      YY(1,KL)  = P(K,1)*P(L,1)*TWO
      YY(2,KL)  = P(K,1)*P(L,2)+P(K,2)*P(L,1)
      YY(3,KL)  = P(K,2)*P(L,2)*TWO
      YY(4,KL)  = P(K,1)*P(L,3)+P(K,3)*P(L,1)
      YY(5,KL)  = P(K,2)*P(L,3)+P(K,3)*P(L,2)
   30 YY(6,KL)  = P(K,3)*P(L,3)*TWO
      IF(IORBS.LE.4 .AND. JORBS.LE.4) RETURN
C     D-S
      DO 40 K=1,5
      KL        = INDX(K+4)+1
      YY(1,KL)  = D(K,1)
      YY(2,KL)  = D(K,2)
      YY(3,KL)  = D(K,3)
      YY(4,KL)  = D(K,4)
   40 YY(5,KL)  = D(K,5)
C     D-P
      DO 50 K=1,5
      DO 50 L=1,3
      KL        = INDX(K+4)+L+1
      YY(1,KL)  = D(K,1)*P(L,1)
      YY(2,KL)  = D(K,1)*P(L,2)
      YY(3,KL)  = D(K,1)*P(L,3)
      YY(4,KL)  = D(K,2)*P(L,1)
      YY(5,KL)  = D(K,2)*P(L,2)
      YY(6,KL)  = D(K,2)*P(L,3)
      YY(7,KL)  = D(K,3)*P(L,1)
      YY(8,KL)  = D(K,3)*P(L,2)
      YY(9,KL)  = D(K,3)*P(L,3)
      YY(10,KL) = D(K,4)*P(L,1)
      YY(11,KL) = D(K,4)*P(L,2)
      YY(12,KL) = D(K,4)*P(L,3)
      YY(13,KL) = D(K,5)*P(L,1)
      YY(14,KL) = D(K,5)*P(L,2)
   50 YY(15,KL) = D(K,5)*P(L,3)
C     D-D
      DO 60 K=1,5
      KL        = INDX(K+4)+K+4
      YY(1,KL)  = D(K,1)*D(K,1)
      YY(2,KL)  = D(K,1)*D(K,2)
      YY(3,KL)  = D(K,2)*D(K,2)
      YY(4,KL)  = D(K,1)*D(K,3)
      YY(5,KL)  = D(K,2)*D(K,3)
      YY(6,KL)  = D(K,3)*D(K,3)
      YY(7,KL)  = D(K,1)*D(K,4)
      YY(8,KL)  = D(K,2)*D(K,4)
      YY(9,KL)  = D(K,3)*D(K,4)
      YY(10,KL) = D(K,4)*D(K,4)
      YY(11,KL) = D(K,1)*D(K,5)
      YY(12,KL) = D(K,2)*D(K,5)
      YY(13,KL) = D(K,3)*D(K,5)
      YY(14,KL) = D(K,4)*D(K,5)
   60 YY(15,KL) = D(K,5)*D(K,5)
      DO 70 K=2,5
      DO 70 L=1,K-1
      KL        = INDX(K+4)+L+4
      YY(1,KL)  = D(K,1)*D(L,1)*TWO
      YY(2,KL)  = D(K,1)*D(L,2)+D(K,2)*D(L,1)
      YY(3,KL)  = D(K,2)*D(L,2)*TWO
      YY(4,KL)  = D(K,1)*D(L,3)+D(K,3)*D(L,1)
      YY(5,KL)  = D(K,2)*D(L,3)+D(K,3)*D(L,2)
      YY(6,KL)  = D(K,3)*D(L,3)*TWO
      YY(7,KL)  = D(K,1)*D(L,4)+D(K,4)*D(L,1)
      YY(8,KL)  = D(K,2)*D(L,4)+D(K,4)*D(L,2)
      YY(9,KL)  = D(K,3)*D(L,4)+D(K,4)*D(L,3)
      YY(10,KL) = D(K,4)*D(L,4)*TWO
      YY(11,KL) = D(K,1)*D(L,5)+D(K,5)*D(L,1)
      YY(12,KL) = D(K,2)*D(L,5)+D(K,5)*D(L,2)
      YY(13,KL) = D(K,3)*D(L,5)+D(K,5)*D(L,3)
      YY(14,KL) = D(K,4)*D(L,5)+D(K,5)*D(L,4)
   70 YY(15,KL) = D(K,5)*D(L,5)*TWO
      RETURN
      END
