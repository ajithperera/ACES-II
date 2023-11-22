
      SUBROUTINE UPDAT1(T1INC,T1OLD,NSIZ,LIST1,LIST2)
C
C THIS ROUTINE UPDATES THE T1 VECTOR (OR OTHER TWO-INDEX QUANTITY)
C  WITH AN INCREMENT.
C
C       T1INC - INCREMENT TO BE ADDED TO VALUES ON LIST (LIST1,LIST2).
C       T1OLD - SCRATCH ARRAY TO HOLD ORIGINAL VALUES.
C       NSIZ  - THE LENGTH OF THE LIST.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T1OLD(NSIZ),T1INC(NSIZ)
      DATA ONE /1.0/
      CALL GETLST(T1OLD,1,1,1,LIST1,LIST2)
      CALL SAXPY(NSIZ,ONE,T1INC,1,T1OLD,1)
      CALL PUTLST(T1OLD,1,1,1,LIST1,LIST2)
      RETURN
      END
