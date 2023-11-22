      SUBROUTINE WRITECORLNPOT(NRADPT, NANGPTS, NCNTR, IUATMS, NATOMS,
     &                         XGRDPT, YGRDPT, ZGRDPT, COORD, 
     &                         CORLN_POT)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*10 FILENAME
      CHARACTER*2 LABEL(9)
C
      DIMENSION XGRDPT(NATOMS,NRADPT,194), YGRDPT(NATOMS,NRADPT,194),
     &          ZGRDPT(NATOMS,NRADPT,194)
C  
      DATA LABEL /"1", "2", "3", "4", "5", "6", "7", "8", "9"/
C
C This simply write out a formatted file that can be viewed by
C using a standard plotting program. The variables are radial distance
C and correlation potential corresponding to that point.
C
      IF (NCNTR .GT. 10) THEN
         WRITE(6, *) "@Writecorpot - Two Many Symmetry Unique Atoms"
         CALL ERREX
      ENDIF
C
      FILENAME = "CORLNPOT"//LABEL(NCNTR)
      IUNIT = 9 + NCNTR
C
      OPEN (UNIT=IUNIT, FILE=FILENAME, FORM="FORMATTED")
C
      DO IRADPT = 1, NRADPT
         
         DO IANGPT = 1, NANGPTS
C
            XPOINT = XGRDPT(NCNTR, IRADPT, IANGPT)  + 
     &                      COORD((NCNTR-1)*3 + 1)
            YPOINT = YGRDPT(NCNTR, IRADPT, IANGPT)  + 
     &                      COORD((NCNTR-1)*3 + 2)
            ZPOINT = ZGRDPT(NCNTR, IRADPT, IANGPT)  +
     &                      COORD((NCNTR-1)*3 + 3)

            R = DSQRT(XPOINT**2 + YPOINT**2 + ZPOINT**2)

            WRITE(IUNIT, 100) R, CORLN_POT
C
         ENDDO
      ENDDO
C
 100  FORMAT(1X, F12.6, 3X, F12.6)
C
      CLOSE(IUNIT)
C
      RETURN
      END
