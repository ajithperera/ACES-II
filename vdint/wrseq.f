      SUBROUTINE WRSEQ(LUFIL,BUF,IBUF,NUT,LEN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BUF(LEN),IBUF(LEN)
      WRITE(LUFIL)BUF,IBUF,NUT
      NUT=0
      RETURN
      END
