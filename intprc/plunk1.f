      SUBROUTINE PLUNK1(LUFIL,BUF,IBUF,ICHAIN,NUT,LEN,lent,NREC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BUF(LENt),IBUF(LENt)
c      IF(MOD(NREC,2).EQ.0) THEN
c       MREC=NREC/2
C
C DETERMINE ACTUAL NUMBER OF RECORDS REQUIRED
C
       MAXREC=(LENT-1)/LEN+1
       IOFF=1
       NREC0=NREC
       DO 100 IREC=1,MAXREC
        JREC=NREC+1
        MUT=MIN(NUT,LEN)  
        NUT=NUT-MUT
        IF(NUT.EQ.0) JREC=ICHAIN
c        WRITE(LUFIL,REC=NREC)BUF(IOFF),IBUF(IOFF),MUT,JREC
        call plunk2(lufil,buf(ioff),ibuf(ioff),len,mut,
     &              nrec,jrec)
        NREC=NREC+1
        IOFF=IOFF+LEN
        IF(NUT.EQ.0) THEN
         ICHAIN=NREC0
         RETURN
        ENDIF
100    CONTINUE
c      ELSE
c       MREC=NREC/2+1
c       WRITE(LUFIL+1,REC=MREC)BUF,IBUF,NUT,ICHAIN
c      ENDIF
      RETURN
      END
