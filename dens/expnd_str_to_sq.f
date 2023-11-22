      SUBROUTINE EXPND_STR_TO_SQ(TRANGL, SQUARE, NTOTAL, NBFIRR, NIRREP)

C
C This is a very simple utility routine to extract a square matrix from 
C trangular symmetry pack vector. Note that this work only for 
C symmetric matrices at the moment, and extensions to other cases
C is a trivial matter. Ajith Perera 07/2000.
C
C  TRANGL : input symmetry packed array. The overall dimension is 
C           N*(N +1)/2 (where N is the total number of functions).
C           The (n*(n+1)/2 elements (n is the number of functions 
C           for each irrep) are stored in upper triangular pack form.
C
C   SQUARE : The standard N X N matrix. 
C
C   NIRREP, NBFIRR : are the number of irreps and the number of 
C                    basis functions for each irrep.        
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
      DIMENSION SQUARE(NTOTAL, NTOTAL), TRANGL(NTOTAL*(NTOTAL+1)/2),
     &          NBFIRR(NIRREP)
C
      DATA IONE /1/, ITWO /2/, IZERO /0/
C
      IBEGNR = IONE
      IBEGNL = IONE
      IBEGN  = IONE
      IOFFST = IZERO
C
      CALL ZERO(SQUARE, NTOTAL*NTOTAL)

      DO 10 IRREP = 1, NIRREP
C
         NELMNTS = NBFIRR(IRREP)
         NELMNTP  = IZERO
         IF (IRREP .NE. 1) NELMNTP = NBFIRR(IRREP - 1)
C
         IBEGNR = IBEGNR + NELMNTP
         IBEGNL = IBEGNL + NELMNTP
C     
         DO 20 IELMNTS = 1, NELMNTS
C    
            IOOFR  = IBEGNL 
            IOOFL  = IBEGNR +  (IELMNTS - IONE)
            LENGTH = IELMNTS
C
            CALL SCOPY(LENGTH, TRANGL(IBEGN), IONE, SQUARE(IOOFR, 
     &                 IOOFL), IONE)
C
            IOFFST = IOFFST + LENGTH
            IBEGN  = IOFFST + IONE
C
 20      CONTINUE
C
 10   CONTINUE
C
C Symmetrize the square matrix. 
C
      DO 30 IROW = 1, NTOTAL
         DO 40 ICOL = 1, NTOTAL
C
            IF (ICOL .NE. IROW) SQUARE(ICOL, IROW) = SQUARE(IROW, ICOL)
C
 40      CONTINUE
 30   CONTINUE

      RETURN
      END
