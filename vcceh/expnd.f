C
      SUBROUTINE EXPND(WPACK, WFULL, NDIM, IANTI)
C
C This routine expands a taringular packed vector
C of numbers into a square matrix
C 
C  WPACK(NDIM*(NDIM+1))/2) ----> WFULL(NDIM, NDIM)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION WPACK((NDIM*(NDIM + 1))/2), WFULL(NDIM, NDIM)
C
      ITHRU = 0
C
      DO 10 I = 1, NDIM
CDIR$ IVDEP
*VOCL LOOP,NOVREC
C
         DO 20 J = 1, I
            ITHRU = ITHRU + 1
            WFULL(I, J) = WPACK(ITHRU)
C
            IF (IANTI .EQ. 1) THEN
               WFULL(J, I) = -WPACK(ITHRU)
            ELSE
               WFULL(J, I) = WPACK(ITHRU)
            ENDIF
C     
 20      CONTINUE
 10   CONTINUE
C
      RETURN
      END
