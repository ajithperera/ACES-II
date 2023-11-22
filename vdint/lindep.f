      SUBROUTINE LINDEP(AMAT,THRESH,OK,NROW,NCOL,NOLD,NNEW,IPRINT)
C
C     Checks input vector for linear dependency against old vectors
C     and adds input vector to old ones if independent.
C
C     tuh Nov 1988
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      LOGICAL OK
      DIMENSION AMAT(NROW,NCOL)
C
      IF (NOLD + NNEW .GT. NCOL) THEN
         OK = .FALSE.
      ELSE
         OK = .TRUE.
         DO 100 INEW = NOLD + 1, NOLD + NNEW
            DO 200 IOLD = 1, INEW - 1
               COEF = - SDOT(NROW,AMAT(1,IOLD),1,AMAT(1,INEW),1)
               CALL SAXPY(NROW,COEF,AMAT(1,IOLD),1,AMAT(1,INEW),1)
  200       CONTINUE
            PROD = SDOT(NROW,AMAT(1,INEW),1,AMAT(1,INEW),1)
            PROD = SQRT(PROD)
            IF (PROD .GT. THRESH) THEN
               FACTOR = D1/PROD
            ELSE
               FACTOR = D0
            END IF
            CALL SSCAL(NROW,FACTOR,AMAT(1,INEW),1)
            OK = OK .AND. PROD .GT. THRESH
  100    CONTINUE
      END IF
      IF (IPRINT .GE. 15) THEN
         CALL TITLER('Output from LINDEP','*',103)
         WRITE (LUPRI,'(A,2I5)') ' NROW, NCOL  ', NROW, NCOL
         WRITE (LUPRI,'(A,2I5)') ' NOLD, NNEW  ', NOLD, NNEW
         WRITE (LUPRI,'(A,L5)') ' OK ', OK
         CALL AROUND('AMAT - LINDEP')
         CALL OUTPUT(AMAT,1,NROW,1,NOLD+NNEW,NROW,NCOL,1,LUPRI)
      END IF
      RETURN
      END
