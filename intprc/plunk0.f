      SUBROUTINE PLUNK0(LUFIL,BUF,IBUF,ICHAIN,NUT,LEN,NREC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BUF(LEN),IBUF(LEN)
      IF(MOD(NREC,2).EQ.0) THEN
       MREC=NREC/2
       WRITE(LUFIL,REC=MREC)BUF,IBUF,NUT,ICHAIN
      ELSE
       MREC=NREC/2+1
       WRITE(LUFIL+1,REC=MREC)BUF,IBUF,NUT,ICHAIN
      ENDIF
      ICHAIN=NREC
      NREC=NREC+1
      NUT=0
      RETURN
      END
