

      SUBROUTINE PLUNK(LUFIL,BUF,IBUF,ICHAIN,NUT,LEN,NREC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BUF(LEN),IBUF(LEN)
      WRITE(LUFIL,REC=NREC)BUF,IBUF,NUT,ICHAIN
      ICHAIN=NREC
      NREC=NREC+1
      NUT=0
      RETURN
      END
