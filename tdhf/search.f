      SUBROUTINE SEARCH (A,NTAPE)
CNCHAR   1
C     REAL*8 A,B,C,D,X
      CHARACTER*8 A,B,C,D,X
C     DIMENSION B(4)
      DATA C/'********'/
      MTAPE=IABS(NTAPE)
      D=A
    1 READ (MTAPE,END=3) X
    2 IF (X.NE.C) GO TO 1
C     BACKSPACE MTAPE
      READ (MTAPE) B
      IF (NTAPE.GE.0) GO TO 5
C   5 IF (B(4).NE.D) GO TO 1
    5 IF (B.NE.D) GO TO 1
      RETURN
    3 WRITE (6,4) A,NTAPE
    4 FORMAT (5X,A8,18H NOT FOUND ON UNIT ,I3,1H.)
      call c_exit(1)
      END
