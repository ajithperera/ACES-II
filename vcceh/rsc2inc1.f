      SUBROUTINE RSC2INC1(ICORE,MAXCOR,IUHF,IRREPX,ISIDE)
C
C This routine calculates the parts of
C   C2 x Hbar -> C1 that get calculated as part of Hbar x C2 -> C2.
C   This routine is for those cases when DT2INT2 is not run.  Much of
C   this code is lifted from DT2INT2.
C SG 6/95 MN 6/95
C
      IMPLICIT NONE
      INTEGER ICORE, MAXCOR, IUHF, IRREPX, ISIDE, LISTT2IN, LISTT2RS
      DIMENSION ICORE(MAXCOR)
C
      DOUBLE PRECISION ONE, ONEM
      PARAMETER (ONE=1.0D0, ONEM=-1.0D0)
      LOGICAL THREEBOD
C
C THREE BODY TERM
C
      THREEBOD = .TRUE.
      IF (THREEBOD) THEN
C
        IF (ISIDE .EQ. 2) THEN
C
C left hand C2 into C1 contribution
C
          CALL GFORMG (IRREPX,1,444,44,400,ICORE,MAXCOR,0,ONE,IUHF)
          CALL GINC1L (IRREPX,ICORE,MAXCOR,IUHF)
        ENDIF
      ENDIF
C
      RETURN
      END
