
C READS A RECORD FROM DIRECT ACCESS UNIT IUNIT.
C
C  INPUT:
C        IUNIT - UNIT NUMBER FOR DIRECT ACCESS FILE.
C        IREC  - NUMBER OF RECORD TO BE READ.
C        LENGTH- LENGTH OF VECTOR IN *INTEGER* WORDS.
C        IMOD  - RETURNED AS ZERO.
C
C        IVEC  - CONTENTS OF RECORD IREC.

      SUBROUTINE RDDIR(IUNIT,IREC,IVEC,LENGTH,IMOD)
      IMPLICIT INTEGER (A-Z)
      DIMENSION IVEC(LENGTH)
      READ(IUNIT,REC=IREC,ERR=555,IOSTAT=IER)IVEC
      IMOD=0
      RETURN
555   CALL IOERR('RDDIR',IUNIT,IER)
      END