      SUBROUTINE SYTRAB1(GAMMA,CL,CR,NAOL,NAOR,NMOL,NMOR,
     &                   BUF,RESULT,IACC)
C
C THIS ROUTINE PERFORMS THE TRANSFORMATION
C
C                                       + 
C                Q(XX) = C(X,P) Z(P,Q) C (X,Q) 
C
C
C USED IN TRANSFORMING THE T2 VECTOR TO THE ATOMIC ORBITAL BASIS
C
C THE RESULT IS RETURNED IN VECTOR RESULT, AND BUF IS A SCRATCH
C  VECTOR OF LENGTH NAOA*NAOB.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CL(NAOL*NMOL),CR(NAOR*NMOR),GAMMA(NMOL*NMOR),BUF(1)
      DIMENSION RESULT(NAOL*NAOR)
C
      DATA ZILCH /0.0/
      DATA ONE   /1.0/
C
      IF(IACC.EQ.0)CALL ZERO(RESULT,NAOL*NAOR)
C
      IF(MIN(NAOL,NAOR).NE.0)THEN
       CALL XGEMM('N','N',NAOL,NMOR,NMOL,ONE,CL,NAOL,GAMMA,NMOL,
     &            ZILCH,BUF,NAOL)
       CALL XGEMM('N','T',NAOL,NAOR,NMOR,ONE,BUF,NAOL,CR,NAOR,
     &            ONE,RESULT,NAOL)
      ENDIF
C
      RETURN
      END
