      SUBROUTINE DETSPINF
C
C IN RGF EOMCC CALCULATIONS WE HAVE TO EVALUATE 'DOT-PRODUCTS' OF THE FORM
C    A = (2 L(IJ,AB) - L(IJ,BA)) R(IJ,AB)
C IF WE DENOTE THE DIRECT TERM LD (OR RD) AND THE EXCHANGE TERM LX (RX)
C WE CAN WRITE THE CONTRIBUTION TO THE DOT PRODUCT FROM THESE TERMS IN
C MATRIX FORMAT
C
C   (LD LX)  [ 2  - 1 ] ( RD)  = (L . SA R)
C            [ -1  2  ] (RX)
C
C THIS PRODUCT CAN BE EVALUATED IN A SYMMETRIC FORM WHICH IS MORE EFFICIENT.
C  DIAGONALIZE SA - > EVEC, EVAL
C  THEN DEFINE RT = SQRT(EVAL) EVEC^-1 R
C              LT = SQRT(EVAL) EVEC^-1 L
C 
C AND SPINADAPTED DOT PRODUCTS CAN BE EVALUATED AS (LT . RT)
C
C  IN THIS ROUTINE THE FACTORS ARE CALCULATED TO PERFORM THIS TRANSFORMATION
C  AND ITS INVERSE
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION L, LT
      DIMENSION SA(2,2), EVEC(2,2), SQRTEVAL(2), DUM(2,2), L(2), R(2),
     &   LT(2),RT(2)
C
      COMMON /SPINF/ T1F, T2FD, T2FX, T1FI, T2FDI, T2FXI
C
      T1F = SQRT(2.0D0)
      T1FI = 1.0D0 / T1F
C
      SA(1,1) = 2.0D0
      SA(1,2) = - 1.0D0
      SA(2,1) = -1.0D0
      SA(2,2) = 2.0D0
C
      CALL SCOPY(4, SA, 1, DUM, 1)
      CALL EIG(DUM, EVEC, 0, 2, 2)
C
C  MAKE SURE THE EIGENVECTORS HAVE THE PROPER SIGN
C
      
      DO I = 1, 2
        SQRTEVAL(I) = SQRT(DUM(I,I))
      ENDDO
C
C MAKE DIAGONAL TRANSFORMATION T = SQRTEVAL * EVEC^-1 (= EVEC^T)
C
      DO I = 1, 2
        DO J = 1,2
          DUM(I,J) = SQRTEVAL(I) * EVEC(J,I)
        ENDDO
      ENDDO
C
      T2FD = DUM(1,1)
      T2FX = DUM(1,2)
C
C  CHECK TRANSFORMATION
C
C  CREATE RANDOM VECTORS L AND R
C
      L(1) = 2.5
      L(2) = 3.1
      R(1) = 4.2
      R(2) = 5.7
C
      A1 = 0.0
      DO I = 1, 2
        DO J=1,2
          A1 = A1 + L(I) * SA(I,J) * R(J)
        ENDDO
      ENDDO
C
C  TRANSFORM L AND R
C
      CALL ZERO(RT, 2)
      CALL ZERO(LT, 2)
C
      DO I = 1, 2
        DO J=1,2
          RT(I) = RT(I) + DUM(I,J) * R(J)
          LT(I) = LT(I) + DUM(I,J) * L(J)
        ENDDO
      ENDDO
C
      A2 = SDOT(2,RT, 1, LT, 1)
      IF (ABS(A1-A2) . GT. 1.0D-8) THEN
        WRITE(6,*) ' SOMETHING WRONG IN DETSPINF'
        CALL ERREX
      ENDIF
C
C DETERMINE INVERSE COEFFICIENTS
C
      DO I = 1, 2
        DO J = 1,2
          DUM(I,J) = EVEC(I,J) / SQRTEVAL(J)
        ENDDO
      ENDDO
C
      T2FDI = DUM(1,1)
      T2FXI = DUM(1,2)      
C
C CHECK INVERTED COEFFICIENTS
C
      CALL ZERO(LT,2)
      DO I = 1, 2
        DO J= 1, 2
          LT(I) = LT(I) + DUM(I,J) * RT(J)
        ENDDO
      ENDDO
C
      DIFF = 0.0
      DO I=1, 2
        DIFF = DIFF + R(I) - LT(I)
      ENDDO
      IF (ABS(DIFF) . GT. 1.0D-8) THEN
        WRITE(6,*)'  SOMETHING WRONG WITH INVERTED DETSPINF'
        CALL ERREX
      ENDIF
C
      WRITE(6,*) ' SPIN-FACTORS IN SPIN-ADAPTATION ALGORITHM'
      WRITE(6,*) ' T1F, T2FD, T2FX ', T1F, T2FD, T2FX
      WRITE(6,*) ' T1FI, T2FDI, T2FXI ', T1FI, T2FDI, T2FXI
      WRITE(6,*)
C
      RETURN
      END
