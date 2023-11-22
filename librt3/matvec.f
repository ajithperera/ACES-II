      SUBROUTINE MATVEC(A,B,C,NI,NJ,IACC,IFUN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NI,NJ),B(NJ),C(NI)
C
      IF(IFUN.EQ.0)THEN
      SIGN =  1.0D+00
      ELSE
      SIGN = -1.0D+00
      ENDIF
C
      IF(IACC.EQ.0)THEN
      FACT = 1.0D+00
      ELSE
      FACT = 0.0D+00
      ENDIF
C
      CALL XGEMM('N','N',NI,1,NJ,
     1           SIGN,
     1           A,NI,B,NJ,FACT,C,NI)
CJ      DO   20 J=1,NJ
CJ      DO   10 I=1,NI
CJ      C(I) = FACT*C(I) + SIGN * A(I,J) * B(J)
CJ   10 CONTINUE
CJ   20 CONTINUE
      RETURN
      END
