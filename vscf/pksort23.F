      SUBROUTINE PKSORT23(IBKNUM,ICHAIN,BUCK,IBUCK,NINBCK,BUF,IBUF,
     &                    IBKNMOF,LDIM1,IRECL,NBKINT,ILNBUF,IBKDIS,
     &                    NBUCK,IUHF)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER(LUINT=10)
      PARAMETER(LUSRT3=23)
C
      INTEGER AND,OR
      CHARACTER*80 FNAME
C
      DIMENSION BUCK(NBKINT,3*NBUCK),IBUCK(NBKINT,3*NBUCK)
      DIMENSION NINBCK(3*NBUCK),ICHAIN(3*NBUCK)
      DIMENSION BUF(ILNBUF),IBUF(ILNBUF)
      DIMENSION IBKNMOF(LDIM1),IBKNUM(IBKDIS)
C
      DIMENSION IWHEREJ(1200),IWHEREK1(1200),IWHEREK2(1200)
      DIMENSION BUFJ(1200),BUFK1(1200),BUFK2(1200)
      DIMENSION IBUFJ(1200),IBUFK1(1200),IBUFK2(1200)
      DIMENSION ILOOKUP(600)
C
      DIMENSION IPKOFF(73)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS/ IFLAGS(100)
#include "symm2.com"
      COMMON /PKOFF/ IPKOFF
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
C
      DATA TWO /2.0/
C
      IUPKI(INT)=AND(INT,IALONE)
      IUPKJ(INT)=AND(ISHFT(INT,-IBITWD),IALONE)
      IUPKK(INT)=AND(ISHFT(INT,-2*IBITWD),IALONE)
      IUPKL(INT)=AND(ISHFT(INT,-3*IBITWD),IALONE)
      IPACK(I,J,K,L)=OR(OR(OR(I,ISHFT(J,IBITWD)),ISHFT(K,2*IBITWD)),
     &                  ISHFT(L,3*IBITWD))
      INDX(I,J)=J+(I*(I-1))/2
      INDX2(I,J,N)=I+(J-1)*N
C
      CALL TIMER(1)
C
      IF(IUHF.EQ.0) THEN
        SCALE=0.25
      ELSE
        SCALE=0.5
      ENDIF
C
      NREC1=1
      NREC2=1
      NREC3=1
      CALL GFNAME('AOK2    ',FNAME,ILENGTH)
      OPEN(UNIT=LUSRT3,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &     ACCESS='DIRECT',RECL=IRECL)
      CALL TIMER(1)
      CALL ZERO(BUCK,NBKINT*3*NBUCK)
      CALL IZERO(IBUCK,NBKINT*3*NBUCK)
      CALL IZERO(NINBCK,3*NBUCK)
      NAOBUF=0
      NAOINT=0
C
C  Read in the AO integrals one buffer at a time.  Make sure that we
C  have the indices in canonical order, i.e. i>=j, k>=l, and ij>=kl.
C
      CALL GFNAME('IIII    ',FNAME,ILENGTH)
      OPEN(LUINT,FILE=FNAME(1:ILENGTH),STATUS='OLD',FORM='UNFORMATTED',
     &     ACCESS='SEQUENTIAL')
      CALL LOCATE(LUINT,'TWOELSUP')
    1 READ(LUINT)BUF,IBUF,NUT
      NAOBUF=NAOBUF+1
      ICNT=0
      CALL ZERO(BUFK2,1200)
      CALL IZERO(IBUFK2,1200)
      CALL IZERO(IWHEREK2,1200)
      DO 10 INT=1,NUT
        IX=IUPKI(IBUF(INT))
        JX=IUPKJ(IBUF(INT))
        KX=IUPKK(IBUF(INT))
        LX=IUPKL(IBUF(INT))
        I=MAX(IX,JX)
        J=MIN(IX,JX)
        K=MAX(KX,LX)
        L=MIN(KX,LX)
        IND1=INDX(I,J)
        IND2=INDX(K,L)
        IF(IND1.LT.IND2) THEN
          IX=I
          JX=J
          I=K
          J=L
          K=IX
          L=JX
        ENDIF
        X=BUF(INT)
C
C  Now deal with the integrals based on what symmetry type it is, i.e.,
C  (aa|aa), (aa|bb), or (ab|ab).
C
        IND3=INDX(INEWVC(I),INEWVC(K))
        IND4=INDX(MAX(INEWVC(J),INEWVC(L)),MIN(INEWVC(J),
     &            INEWVC(L)))
        IND5=INDX(INEWVC(I),INEWVC(L))
        IND6=INDX(MAX(INEWVC(J),INEWVC(K)),MIN(INEWVC(J),
     &            INEWVC(K)))
        IF(IND4.GT.IND3) THEN
          ICRAP=IND3
          IND3=IND4
          IND4=ICRAP
        ENDIF
        IF(IND6.GT.IND5) THEN
          ICRAP=IND5
          IND5=IND6
          IND6=ICRAP
        ENDIF
        IREPOFF=INDX(IDXVEC(I),IDXVEC(I))
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
        IF((I.NE.J).AND.(I.NE.K).AND.(I.NE.L).AND.(J.NE.K).AND.
     &     (J.NE.L).AND.(K.NE.L).AND.(J.LT.K).AND.(J.GT.L)) THEN
C
          IBK=IBKNMOF(IREPOFF)+IND4
          IBKET=2*NBUCK+IBKNUM(IBK)
C
          ICNT=ICNT+1
          IWHEREK2(ICNT)=IBKET
          IBUFK2(ICNT)=IPKOFF(IREPOFF)+
     &                 INDX2(IND3,IND4,IRPDS1(IREPOFF))
          BUFK2(ICNT)=-SCALE*BUF(INT)
C
        ENDIF
C
   10 CONTINUE
C
      DO 440 ITYPE=2*NBUCK+1,3*NBUCK
        CALL IZERO(ILOOKUP,600)
        CALL WHENEQ(ILNBUF,IWHEREK2,1,ITYPE,ILOOKUP,NVAL)
        IOFF=1
        INBUCK=NINBCK(ITYPE)
  441   IDUMP=MIN(NVAL,NBKINT-INBUCK)
        CALL GATHER(IDUMP,BUCK(INBUCK+1,ITYPE),BUFK2,ILOOKUP(IOFF))
        CALL IGATHER(IDUMP,IBUCK(INBUCK+1,ITYPE),IBUFK2,ILOOKUP(IOFF))
        INBUCK=INBUCK+IDUMP
        IF(INBUCK.EQ.NBKINT) THEN
          CALL PLUNK(LUSRT3,BUCK(1,ITYPE),IBUCK(1,ITYPE),ICHAIN(ITYPE),
     &               INBUCK,NBKINT,NREC3)
          NVAL=NVAL-IDUMP
          IOFF=IOFF+IDUMP
          IF(NVAL.NE.0) GOTO 441
        ENDIF
        NINBCK(ITYPE)=INBUCK
  440 CONTINUE
C
      ICNT=0
      CALL ZERO(BUFK2,1200)
      CALL IZERO(IBUFK2,1200)
      CALL IZERO(IWHEREK2,1200)
      DO 11 INT=1,NUT
        IX=IUPKI(IBUF(INT))
        JX=IUPKJ(IBUF(INT))
        KX=IUPKK(IBUF(INT))
        LX=IUPKL(IBUF(INT))
        I=MAX(IX,JX)
        J=MIN(IX,JX)
        K=MAX(KX,LX)
        L=MIN(KX,LX)
        IND1=INDX(I,J)
        IND2=INDX(K,L)
        IF(IND1.LT.IND2) THEN
          IX=I
          JX=J
          I=K
          J=L
          K=IX
          L=JX
        ENDIF
        X=BUF(INT)
C
C  Now deal with the integrals based on what symmetry type it is, i.e.,
C  (aa|aa), (aa|bb), or (ab|ab).
C
        IND3=INDX(INEWVC(I),INEWVC(K))
        IND4=INDX(MAX(INEWVC(J),INEWVC(L)),MIN(INEWVC(J),
     &            INEWVC(L)))
        IND5=INDX(INEWVC(I),INEWVC(L))
        IND6=INDX(MAX(INEWVC(J),INEWVC(K)),MIN(INEWVC(J),
     &            INEWVC(K)))
        IF(IND4.GT.IND3) THEN
          ICRAP=IND3
          IND3=IND4
          IND4=ICRAP
        ENDIF
        IF(IND6.GT.IND5) THEN
          ICRAP=IND5
          IND5=IND6
          IND6=ICRAP
        ENDIF
        IREPOFF=INDX(IDXVEC(I),IDXVEC(I))
C
        IF((J.EQ.K).AND.(I.NE.J).AND.(I.NE.L).AND.(J.NE.L)) THEN
          SCALEN=TWO*SCALE
        ELSE
          SCALEN=SCALE
        ENDIF
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C
        IF((((I.NE.J).AND.(I.NE.K).AND.(I.NE.L).AND.(J.NE.K).AND.
     &      (J.NE.L).AND.(K.NE.L).AND..NOT.((J.LT.K).AND.(J.GT.L))).OR.
     &      ((I.EQ.J).AND.(I.EQ.K).AND.(I.EQ.L)).OR.
     &      ((I.EQ.K).AND.(J.EQ.L).AND.(I.NE.J)).OR.
     &      ((J.EQ.K).AND.(J.EQ.L).AND.(I.NE.J)).OR.
     &      ((I.EQ.J).AND.(I.EQ.K).AND.(I.NE.L)).OR.
     &      ((J.EQ.K).AND.(I.NE.J).AND.(I.NE.L).AND.(J.NE.L)).OR.
     &      ((I.NE.K).AND.(J.EQ.L).AND.(I.NE.J)).OR.
     &      ((I.EQ.K).AND.(I.NE.J).AND.(J.NE.L)))) THEN
C
          IBK=IBKNMOF(IREPOFF)+IND6
          IBKET=2*NBUCK+IBKNUM(IBK)
C
          ICNT=ICNT+1
          IWHEREK2(ICNT)=IBKET
          IBUFK2(ICNT)=IPKOFF(IREPOFF)+
     &                 INDX2(IND5,IND6,IRPDS1(IREPOFF))
          BUFK2(ICNT)=-SCALEN*BUF(INT)
C
        ENDIF
C
   11 CONTINUE
C
      DO 442 ITYPE=2*NBUCK+1,3*NBUCK
        CALL IZERO(ILOOKUP,600)
        CALL WHENEQ(ILNBUF,IWHEREK2,1,ITYPE,ILOOKUP,NVAL)
        IOFF=1
        INBUCK=NINBCK(ITYPE)
  443   IDUMP=MIN(NVAL,NBKINT-INBUCK)
        CALL GATHER(IDUMP,BUCK(INBUCK+1,ITYPE),BUFK2,ILOOKUP(IOFF))
        CALL IGATHER(IDUMP,IBUCK(INBUCK+1,ITYPE),IBUFK2,ILOOKUP(IOFF))
        INBUCK=INBUCK+IDUMP
        IF(INBUCK.EQ.NBKINT) THEN
          CALL PLUNK(LUSRT3,BUCK(1,ITYPE),IBUCK(1,ITYPE),ICHAIN(ITYPE),
     &               INBUCK,NBKINT,NREC3)
          NVAL=NVAL-IDUMP
          IOFF=IOFF+IDUMP
          IF(NVAL.NE.0) GOTO 443
        ENDIF
        NINBCK(ITYPE)=INBUCK
  442 CONTINUE
C
      IF(NUT.EQ.ILNBUF) GOTO 1
C
      NAOINT=(NAOBUF-1)*ILNBUF+NUT
      CLOSE(LUINT,STATUS='KEEP')
C
  200 WRITE(LUOUT,5000)NAOINT
 5000 FORMAT(T3,'@PKSORT23-I, ',I8,' AO integrals sorted as PK ',
     &          'list elements.')
C
C  Flush remaining buffers
C
      DO 50 I=1,NBUCK
        IBK3=2*NBUCK+I
        CALL PLUNK(LUSRT3,BUCK(1,IBK3),IBUCK(1,IBK3),ICHAIN(IBK3),
     &             NINBCK(IBK3),NBKINT,NREC3)
   50 CONTINUE
C
      CALL TIMER(1)
C
      NREC3=NREC3-1
      WRITE(LUOUT,5010)NREC3, TIMENEW
 5010 FORMAT(T3,'@PKSORT23-I, Wrote out ',I8,' records of K ',
     &          'integrals.',/,
     &       T14,'Integral sort required ',F10.3,' seconds.',/)
C
C  Now we need to reopen the IIII file because other routines need info
C  from this file.
C
      CALL GFNAME('IIII    ',FNAME,ILENGTH)
      OPEN(LUINT,FILE=FNAME(1:ILENGTH),STATUS='OLD',FORM='UNFORMATTED',
     &     ACCESS='SEQUENTIAL')
C
      RETURN
      END
