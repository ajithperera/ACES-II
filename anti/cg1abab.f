      SUBROUTINE CG1ABAB(W,BUF,IBUF,ICHAIN,DISSIZ,NDBCK,
     &                   NBKINT,NBUCK,IRECL,LIST1,LIST2,
     &                   IRREPDO,ALASKA)
C
C PICKS UP ABAB MO GAMMAS FROM SORT FILE AND WRITES THEM TO MOINTS.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL MOEQAO,ALASKA
      INTEGER DISSIZ,DIRPRD,VRT,POP,DISFIRST,DISLEFT,DISWRITE
      INTEGER AND
      DIMENSION W(DISSIZ,NDBCK),BUF(NBKINT),IBUF(NBKINT)
      DIMENSION ICHAIN(NBUCK)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON/SYM/POP(8,2),VRT(8,2),NJUNK(6)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON /SIZES/ MOEQAO
C
      IUPKI(INT)=AND(INT,IALONE)
      IUPKJ(INT)=AND(ISHFT(INT,-IBITWD),IALONE)
C
C LOOP OVER THE NUMBER OF BUCKETS
C
      LUSRT=20
      NUMDIS=DISSIZ
      DISFIRST=1
      DO 10 IDIS=1,NBUCK
       NOFF=(IDIS-1)*NDBCK
       CALL ZERO(W,DISSIZ*NDBCK)
C
C CHAIN IN A BUCKETFULL OF GAMMAS HERE
C
1      READ(LUSRT,REC=ICHAIN(IDIS))BUF,IBUF,NUT,ICHAN
       ICHAIN(IDIS)=ICHAN
       DO 20 INT=1,NUT
        I=IUPKI(IBUF(INT))
        J=IUPKJ(IBUF(INT))-NOFF
        W(I,J)=W(I,J)+BUF(INT)
20     CONTINUE
       IF(ICHAIN(IDIS).NE.0)GOTO 1
C
C  IF WE REACH THIS POINT, THE ENTIRE BUCKET IS IN CORE.  DUMP THEM
C  TO THE MOINTS FILE.
C
       DISLEFT=NUMDIS+1-DISFIRST
       DISWRITE=MIN(NDBCK,DISLEFT)
CJDW 9/13/96. See JFS notes.
C      IF(LIST2.NE.1.OR.MOEQAO)THEN
       IF(MOEQAO .OR. ALASKA)THEN
        CALL PUTLST(W,DISFIRST,DISWRITE,1,LIST1,LIST2)
       ELSE
CJDW 9/13/96. See JFS notes.
C       DO 21 I=DISFIRST,DISWRITE
        DO 21 I=1,DISWRITE
         IDISREAL=DISFIRST+I-1
         CALL PUTLST(W(1,I),IDISREAL,1,1,LIST1,LIST2)
21      CONTINUE
       ENDIF       
       DISFIRST=DISFIRST+NDBCK
10    CONTINUE
C
      CLOSE(UNIT=LUSRT,STATUS='DELETE')
      RETURN
      END
