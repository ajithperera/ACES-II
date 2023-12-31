      SUBROUTINE PRINTALL(A,N, DIS, NUM)
C
C  PRINTS OUT AL IRREPS OF DOUBLE PRECISION ARRAY A
C  DIS AND NUM GIVE THE DISTRIBUTION SIZE AND THE NUMBER OF 
C  DISTRIBUTIONS IN A RESPECTIVELY, PER IRREP.
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION A
      DIMENSION A(N), DIS(8), NUM(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
C
      ISTART = 1
      DO 10 XIRREP = 1, NIRREP
         NR = DIS(XIRREP)
         NC = NUM(XIRREP)
         IF ((NR * NC) .GT. 0) THEN
            WRITE(6,*) '  XIRREP :  ', XIRREP
            CALL OUTPUT(A(ISTART), 1, NR, 1, NC, NR, NC, 1)
            ISTART = ISTART + NR * NC
         ENDIF
 10   CONTINUE
C
      RETURN
      END
