
      SUBROUTINE READD(DOOA,DOOB,DVVA,DVVB,DVOA,DVOB,IUHF)
C
C  THIS ROUTINE READS THE RELAXED DENSITY MATRIX FROM
C  THE GAMLAM FILE
C
CEND
C
C CODED SEP/91 JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER POP,VRT
C
      DIMENSION DOOA(1),DOOB(1),DVVA(1),DVVB(1),
     &          DVOA(1),DVOB(1)
C
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
C
C READ IN RELAXED DENSITY MATRICES  AND INTERMEDIATES
C
      CALL GETLST(DOOA,1,1,1,1,160)
      CALL GETLST(DVVA,1,1,1,3,160)
      CALL GETLST(DVOA,1,1,1,5,160)
      IF(IUHF.NE.0) THEN
       CALL GETLST(DOOB,1,1,1,2,160)
       CALL GETLST(DVVB,1,1,1,4,160)
       CALL GETLST(DVOB,1,1,1,6,160)
      ENDIF
C
      RETURN
C
      END 
