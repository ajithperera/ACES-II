      SUBROUTINE DUMPIT(MATRIX, LENGTH)
C
C Prints out the first elements up to and including LENGTH of MATRIX 
C
      DIMENSION MATRIX(*)
C
      CALL PRVECR(MATRIX(1),LENGTH)
C
      RETURN
      END

