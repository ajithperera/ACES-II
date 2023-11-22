      SUBROUTINE DEFOPN(LU)
C
C 29-Sep-1986 Hans Joergen Aa. Jensen
C
C DEFAULT OPEN :
C
C     OPEN(LU,STATUS='UNKNOWN',FORM='UNFORMATTED') with
C     default file name, following CWBOPN by C.W.Bauschlicher.
C
C     This is especially useful for CRAY machines which do not
C     allow for default names in the OPEN statement (i.e.
C     FILE='name' must be specified if you use OPEN).
C
      REWIND (LU)
      READ   (LU,ERR=1,END=1) A
    1 REWIND (LU)
      RETURN
      END
