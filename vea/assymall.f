C
C **********************************************************
C  GENERAL PROCEDURES TO MANIPULATE WITH LISTS THAT ARE 
C  NOT INCLUDED IN THE ACES II UTILITIES LIBRARY
C **********************************************************
C
      SUBROUTINE ASSYMALL(Z, N, NUM, DIS, POP1, ICORE, MAXCOR)
C     
C  GIVEN IS A COMPLETE LIST (ALL IRREPS) IN ARRAY Z. THIS ARRAY IS
C  ANTISYMMETRIZED (PER IRREP) BUT NOT SQUEEZED!
C
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION Z
      logical print
      DIMENSION Z(N), NUM(8), DIS(8), POP1(8), ICORE(MAXCOR)
C
C  NUM AND DIS CONTAIN THE LENGTH OF THE ROWS AND COLUMNS, PER IRREP OF THE
C  MATRIX READ IN, N IS THE TOTAL LENGHT OF Z
C  POP1 CONTAINS THE POPULATION OF THE INDICES TO BE ANTISYMMETRIZED
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
C
      print = .false.
           ISTART = 1
           I000 = 1
           DO 10 IRREP = 1, NIRREP
              NUMDSW = NUM(IRREP)
              DISSYW = DIS(IRREP)
              I010 = I000 + NUMDSW * IINTFP
              I020 = I010 + NUMDSW * IINTFP
              IF (I020 .GT. MAXCOR) THEN
                 WRITE(6,*) 'INSUFFICIENT MEMORY ASSYMALL'
                 CALL ERREX
              ENDIF
              IF ((NUMDSW*DISSYW).GT.0) THEN
                 if (print) then
                    write(6,*) 'before antisymmetrization,
     $                 irrep = ', irrep
                    call output(Z(istart), 1, dissyw, 1, numdsw,
     $                 dissyw, numdsw, 1)
                 endif
              CALL ASSYM2A(IRREP, POP1, DISSYW, NUMDSW, Z(ISTART),
     $              ICORE(I000), ICORE(I010))
              if (print) then
                 write(6,*)' after antisymmetrization'
                 call output(Z(ISTART),1, dissyw, 1, numdsw,
     $              dissyw, numdsw, 1)
              endif
              ISTART = ISTART + NUMDSW * DISSYW
              IF (ISTART.GT.N+1) THEN
                 WRITE(6,*)'SOMETHING WRONG IN PUTEXP'
                 CALL ERREX
              ENDIF
              ENDIF
 10        CONTINUE
C
           RETURN
           END
