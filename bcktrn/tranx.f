      SUBROUTINE TRANX(ICORE,MAXCOR,XROW,XCOL,UNITX,UNITXT)
C
C THIS ROUTINE TRANSPOSES A GENERAL MATRIX X(A,B) TO X(B,A).
C ON ENTRY, EACH COLUMN OF X IS STORED ON A DIRECT ACCESS RECORD
C LABELLED BY THE ROW INDEX (B).  THIS ROUTINE PROCEEDS TO WRITE
C THE TRANSPOSED MATRIX TO FILE NUMBER UNITRN.  THIS ROUTINE IS
C NOT A GENERAL YOSHIMINE SORT, BUT RATHER ASSUMES THAT ONLY A
C FEW PASSES ARE REQUIRED.
C
CEND
      IMPLICIT INTEGER (A-Z)
      LOGICAL PASS1
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      PASS1=.TRUE.
      NSTART=1
      NPASS=0
C
C CALCULATE HOW MUCH OF X CAN BE HELD IN CORE SIMULTANEOUSLY
C
      RECLEFT=XCOL
      RECSIZX=XROW*IINTFP
      RECSIZXT=XCOL*IINTFP
      NUMRECX=XCOL
      NUMRECXT=XROW
      WRKSPACE=MAXCOR-RECSIZXT
      IWORK0=WRKSPACE*IINTFP+1
      IWORK =IWORK0
      MAXREAD=WRKSPACE/RECSIZX
1     NREAD=MIN(MAXREAD,RECLEFT)
      NPASS=NPASS+1
C
C READ IN NUMREC RECORDS OF INPUT MATRIX
C
      IOFF=1
      DO 10 I=NSTART,NSTART+NREAD-1
       CALL ACES_IO_READ(UNITX,I,ICORE(IOFF),RECSIZX)
       IOFF=IOFF+RECSIZX
10    CONTINUE
      NSTART=NSTART+NREAD
      IOFF=1
      DO 20 J=1,XROW
       IF(.NOT.PASS1)THEN
        CALL ACES_IO_READ(UNITXT,J,ICORE(IWORK0),RECSIZXT)
       ELSE
        CALL ZERO(ICORE(IWORK0),XCOL)
       ENDIF
       CALL SCOPY(NREAD,ICORE(IOFF),XROW,ICORE(IWORK),1)
       CALL ACES_IO_WRITE(UNITXT,J,ICORE(IWORK0),RECSIZXT)
       IOFF=IOFF+IINTFP
20    CONTINUE
      PASS1=.FALSE.
      IWORK=IWORK+NREAD*IINTFP
      RECLEFT=RECLEFT-NREAD
      IF(RECLEFT.NE.0)GOTO 1
C
      WRITE(6,100)NPASS
100   FORMAT(T3,'@TRANX-I, Transposition of half-transformed ',
     &          'gammas required ',I5,' passes.')
C
      RETURN
      END
