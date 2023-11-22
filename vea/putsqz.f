      SUBROUTINE PUTSQZ(Z, N, NUM, DIS, LIST, POP1, NUMB, ICORE, MAXCOR)
C     
C  GIVEN IS A COMPLETE LIST (ALL IRREPS) IN ARRAY Z. THIS ARRAY IS
C  ANTISYMMETRY-COMPRESSED (PER IRREP) AND WRITTEN TO DISK.
C
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION Z
      logical print
      DIMENSION Z(N), NUM(8), DIS(8), NUMB(8), POP1(8), ICORE(MAXCOR)
C
C  NUM AND DIS CONTAIN THE LENGTH OF THE ROWS AND COLUMNS, PER IRREP OF THE
C  MATRIX READ IN, BEFORE COMPRESSION. N IS THE TOTAL LENGHT OF Z
C  NUMB CONTAINS THE LENGHTHS OF THE ROW AFTER COMPRESSION AND POP1
C  CONTAINS THE POPULATION OF THE INDICES TO BE COMPRESSED.
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
C
      print = .false.
           ISTART = 1
           I000 = 1
           DO 10 IRREP = 1, NIRREP
              DISSYB = NUMB(IRREP)
              NUMDSW = NUM(IRREP)
              DISSYW = DIS(IRREP)
              I010 = I000 + NUMDSW * IINTFP
              IF ((NUMDSW*DISSYB).GT.0) THEN
                 if (print) then
                    write(6,*) 'before contraction, irrep = ', irrep
                    call output(Z(istart), 1, dissyw, 1, numdsw,
     $                 dissyw, numdsw, 1)
                 endif
              CALL SQSYM(IRREP, POP1,DISSYB, DISSYW, NUMDSW,
     $              ICORE(I000), Z(ISTART))
              if (print) then
                 write(6,*)' after antisymmetrization'
                 call output(ICORE(I000),1, dissyb, 1, numdsw,
     $              dissyb, numdsw, 1)
              endif
              CALL PUTLST(ICORE(I000), 1, NUMDSW, 1, IRREP, LIST)
              ENDIF
              ISTART = ISTART + NUMDSW * DISSYW
              IF (ISTART.GT.N+1) THEN
                 WRITE(6,*)'SOMETHING WRONG IN PUTSQZ'
                 CALL ERREX
              ENDIF
 10        CONTINUE
C
           RETURN
           END
