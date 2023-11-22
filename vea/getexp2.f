      SUBROUTINE GETEXP2(Z, N, NUM, DIS, LIST, POP1, DISB)
C     
C  READ IN THE COMPLETE LIST, ANTISYMMETRIC IN THE ROW-INDICES 
C  AND SIMULTANEOUSLY EXPAND THE LIST SUCH THAT ALL ENTRIES ARE GIVEN
C
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION Z
      logical print
      DIMENSION Z(N), NUM(8), DIS(8), DISB(8), POP1(8)
C
C  NUM AND DIS CONTAIN THE LENGTH OF THE ROWS AND COLUMNS, PER IRREP OF THE
C  MATRIX TO BE READ IN, AFTER EXPANSION. POP1 IS THE POPULATION OF THE
C  INDICES TO BE ANTISYMMETRIZED, DISB IS THE LENGTH OF THE DISTRIBUTIONS
C  IN ANTISYMMETRIZED FORM.
C
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
C
      print = .false.
           ISTART = 1
              if (print) then
                 write(6,*)' getexp2 info is printed'
              endif
           DO 10 IRREP = 1, NIRREP
              NUMDSW = NUM(IRREP)
              DISSYW = DIS(IRREP)
              IF ((NUMDSW*DISSYW).GT.0) THEN
              DISA = DISB(IRREP)
              IF (DISA .GT. 0) THEN
                 CALL GETLST(Z(ISTART), 1, NUMDSW, 1, IRREP, LIST)
              ENDIF
              CALL SYMEXP2(IRREP, POP1, DISSYW, DISA, NUMDSW,
     $           Z(ISTART), Z(ISTART))
              if (print) then
                 write(6,*)' final expanded array, irrep = ', irrep
                 call output(Z(istart), 1, dissyw, 1, numdsw,
     $              dissyw, numdsw, 1)
              endif
              ISTART = ISTART + NUMDSW * DISSYW
              IF (ISTART.GT.N+1) THEN
                 WRITE(6,*)'SOMETHING WRONG IN GETEXP'
                 CALL ERREX
              ENDIF
              ENDIF
 10        CONTINUE
C
           RETURN
           END
