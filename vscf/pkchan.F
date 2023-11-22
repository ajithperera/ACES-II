      SUBROUTINE PKCHAN(ICHAIN,PK,IDISSIZ,IDISNUM,BUF,IBUF,IBKSTRT,
     &                  IBKOF,MAXBUK,IBKDIS,IBKSIZ,NBUCK,NBKINT,IUHF)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER(LUSRT1=21)
      PARAMETER(LUSRT2=22)
      PARAMETER(LUSRT3=23)
C
      DIMENSION ICHAIN(3*NBUCK),IDISNUM(IBKDIS,2)
      DIMENSION PK(IBKSIZ),IDISSIZ(IBKDIS)
      DIMENSION IBKSTRT(MAXBUK),IBKOF(MAXBUK)
      DIMENSION BUF(NBKINT),IBUF(NBKINT)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
#include "symm2.com"
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
C
      INDX(I,J)=J+(I*(I-1))/2
C
      CALL TIMER(1)
C
C  Initialize all of the PK lists before we get the actual numbers.
C
c      CALL UPDMOI(0,0,0,0,1,0)
      call aces_io_remove(50,'MOINTS')
      DO 10 I=1,NIRREP
        DO 11 J=1,I
          IF(I.EQ.J) THEN
            IJUNK=INDX(I,I)
            CALL UPDMOI(IRPDS1(IJUNK),IRPDS1(IJUNK),1,IJUNK,0,0)
            IF(IUHF.NE.0) THEN
              CALL UPDMOI(IRPDS1(IJUNK),IRPDS1(IJUNK),2,IJUNK,0,0)
            ENDIF
          ELSE
            IJUNK=INDX(I-1,J)
            CALL UPDMOI(IRPDS2(2*IJUNK),IRPDS2(2*IJUNK-1),1,
     &                  36+IJUNK,0,0)
            IF(IUHF.NE.0) THEN
              CALL UPDMOI(IRPDS2(2*IJUNK),IRPDS2(2*IJUNK-1),2,
     &                    36+IJUNK,0,0)
            ENDIF
          ENDIF
   11   CONTINUE
   10 CONTINUE
C
      DO 50 IBK=1,NBUCK
        CALL ZERO(PK,IBKSIZ)
        IF(ICHAIN(IBK).EQ.0) GOTO 198
   99   READ(LUSRT1,REC=ICHAIN(IBK))BUF,IBUF,NUT,ICHAN
        ICHAIN(IBK)=ICHAN
        DO 100 INT=1,NUT
          IADR=IBUF(INT)-IBKOF(IBK)
          PK(IADR)=PK(IADR)+BUF(INT)
  100   CONTINUE
        IF(ICHAIN(IBK).NE.0) GOTO 99
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  If this is not an RHF calculation, then dump out the J integrals as
C  a separate list in addition to their contribution to the P list.
C
        IF(IUHF.NE.0) THEN
          ISTART=IBKSTRT(IBK)
          IEND=IBKSTRT(IBK+1)-1
          IOFF=1
          DO 500 I=ISTART,IEND
            CALL PUTLST(PK(IOFF),IDISNUM(I,1),1,2,2,IDISNUM(I,2))
            IOFF=IOFF+IDISSIZ(I)
  500     CONTINUE
        ENDIF
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
  198   IF(ICHAIN(NBUCK+IBK).EQ.0) GOTO 298
  199   READ(LUSRT2,REC=ICHAIN(NBUCK+IBK))BUF,IBUF,NUT,ICHAN
        ICHAIN(NBUCK+IBK)=ICHAN
        DO 200 INT=1,NUT
          IADR=IBUF(INT)-IBKOF(IBK)
          PK(IADR)=PK(IADR)+BUF(INT)
  200   CONTINUE
        IF(ICHAIN(NBUCK+IBK).NE.0) GOTO 199
C
  298   IF(ICHAIN(2*NBUCK+IBK).EQ.0) GOTO 310
  299   READ(LUSRT3,REC=ICHAIN(2*NBUCK+IBK))BUF,IBUF,NUT,ICHAN
        ICHAIN(2*NBUCK+IBK)=ICHAN
        DO 300 INT=1,NUT
          IADR=IBUF(INT)-IBKOF(IBK)
          PK(IADR)=PK(IADR)+BUF(INT)
  300   CONTINUE
        IF(ICHAIN(2*NBUCK+IBK).NE.0) GOTO 299
C
  310   CONTINUE
C
        ISTART=IBKSTRT(IBK)
        IEND=IBKSTRT(IBK+1)-1
        IOFF=1
        DO 400 I=ISTART,IEND
          CALL PUTLST(PK(IOFF),IDISNUM(I,1),1,2,1,IDISNUM(I,2))
          IOFF=IOFF+IDISSIZ(I)
  400   CONTINUE
C
   50 CONTINUE
C
      CLOSE(UNIT=LUSRT1,STATUS='DELETE')
      CLOSE(UNIT=LUSRT2,STATUS='DELETE')
      CLOSE(UNIT=LUSRT3,STATUS='DELETE')
C
      CALL TIMER(1)
      WRITE(LUOUT,5000)TIMENEW
 5000 FORMAT(T3,'@PKCHAN-I, Chaining of PK lists required ',F10.3,
     &          ' seconds.',/)
C
      RETURN
      END
