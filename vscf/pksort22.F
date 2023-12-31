      SUBROUTINE PKSORT22(IBKNUM,ICHAIN,BUCK,IBUCK,NINBCK,BUF,IBUF,
     &                    IBKNMOF,LDIM1,IRECL,NBKINT,ILNBUF,IBKDIS,
     &                    NBUCK,IUHF)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER(LUINT=10)
      PARAMETER(LUSRT2=22)
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
      CALL GFNAME('AOK1    ',FNAME,ILENGTH)
      OPEN(UNIT=LUSRT2,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &     ACCESS='DIRECT',RECL=IRECL)
      CALL TIMER(1)
      CALL ZERO(BUCK,NBKINT*3*NBUCK)
      CALL IZERO(IBUCK,NBKINT*3*NBUCK)
      CALL IZERO(NINBCK,3*NBUCK)
      ICNT=0
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
      ICNT2=ILNBUF
      CALL ZERO(BUFK1,1200)
      CALL IZERO(IBUFK1,1200)
      CALL IZERO(IWHEREK1,1200)
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
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  This handles the most general integral, i.e., those where
C  all four indices are different.
C
        IF((I.NE.J).AND.(I.NE.K).AND.(I.NE.L).AND.(J.NE.K).AND.
     &     (J.NE.L).AND.(K.NE.L).AND.(J.LT.K).AND.(J.GT.L)) THEN
C
          INDA=IND5
          INDB=IND6
C
        ELSE
          INDA=IND3
          INDB=IND4
C
        ENDIF
C
        IF(
     &     (.NOT.((I.NE.J).AND.(I.NE.K).AND.(I.NE.L).AND.(J.NE.K).AND.
     &            (J.NE.L).AND.(K.NE.L)).AND.
     &      .NOT.(((I.EQ.J).AND.(I.EQ.K).AND.(I.EQ.L)).OR.
     &            ((I.EQ.J).AND.(K.EQ.L).AND.(I.NE.K))) ).AND.
     &     ((((I.EQ.K).AND.(J.EQ.L).AND.(I.NE.J)) .OR.
     &       ((I.NE.K).AND.(J.EQ.L).AND.(I.NE.J)) .OR.
     &       ((I.EQ.K).AND.(J.NE.L).AND.(I.NE.J))) .AND..NOT.
     &       ((J.EQ.K).AND.(J.EQ.L).AND.(I.NE.J)))  ) THEN
C
          SCALEN=TWO*SCALE
        ELSE
          SCALEN=SCALE
        ENDIF
C
        IBK=IBKNMOF(IREPOFF)+INDB
        IBKET=NBUCK+IBKNUM(IBK)
C
        ICNT=ICNT+1
        IWHEREK1(ICNT)=IBKET
        IBUFK1(ICNT)=IPKOFF(IREPOFF)+
     &               INDX2(INDA,INDB,IRPDS1(IREPOFF))
        BUFK1(ICNT)=-SCALEN*BUF(INT)
C
   10 CONTINUE
C
C  Now we put the integrals into the appropriate buckets.
C
      CALL IZERO(ILOOKUP,600)
      DO 420 ITYPE=NBUCK+1,2*NBUCK
        CALL WHENEQ(ILNBUF,IWHEREK1,1,ITYPE,ILOOKUP,NVAL)
        IOFF=1
        INBUCK=NINBCK(ITYPE)
  421   IDUMP=MIN(NVAL,NBKINT-INBUCK)
        CALL GATHER(IDUMP,BUCK(INBUCK+1,ITYPE),BUFK1,ILOOKUP(IOFF))
        CALL IGATHER(IDUMP,IBUCK(INBUCK+1,ITYPE),IBUFK1,ILOOKUP(IOFF))
        INBUCK=INBUCK+IDUMP
        IF(INBUCK.EQ.NBKINT) THEN
          CALL PLUNK(LUSRT2,BUCK(1,ITYPE),IBUCK(1,ITYPE),ICHAIN(ITYPE),
     &               INBUCK,NBKINT,NREC2)
          NVAL=NVAL-IDUMP
          IOFF=IOFF+IDUMP
          IF(NVAL.NE.0) GOTO 421
        ENDIF
        NINBCK(ITYPE)=INBUCK
  420 CONTINUE
C
      IF(NUT.EQ.ILNBUF) GOTO 1
C
      NAOINT=(NAOBUF-1)*ILNBUF+NUT
      CLOSE(LUINT,STATUS='KEEP')
C
      IF(NIRREP.EQ.1) GOTO 200
C
      CALL GFNAME('IJIJ    ',FNAME,ILENGTH)
      OPEN(LUINT,FILE=FNAME(1:ILENGTH),STATUS='OLD',FORM='UNFORMATTED',
     &     ACCESS='SEQUENTIAL')
      CALL LOCATE(LUINT,'TWOELSUP')
      NAOBUF=0
    3 READ(LUINT)BUF,IBUF,NUT
      NAOBUF=NAOBUF+1
      CALL ZERO(BUFK1,1200)
      CALL IZERO(IBUFK1,1200)
      CALL IZERO(IWHEREK1,1200)
      DO 30 INT=1,NUT
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
        IF(I.EQ.K.OR.J.EQ.L) THEN
          X=TWO*X
        ENDIF
C
        IND1=INDX(INEWVC(I),INEWVC(K))
        IND2=INDX(MAX(INEWVC(J),INEWVC(L)),MIN(INEWVC(J),
     &            INEWVC(L)))
        IREPOFF=INDX(IDXVEC(I),IDXVEC(J))
        IJUNK=INDX(IDXVEC(I)-1,IDXVEC(J))
        IBK=IBKNMOF(IREPOFF)+IND2
        IBKET=NBUCK+IBKNUM(IBK)
C
        IWHEREK1(INT)=IBKET
        IBUFK1(INT)=IPKOFF(IREPOFF)+INDX2(IND1,IND2,IRPDS2(2*IJUNK-1))
        BUFK1(INT)=-SCALE*X
C
   30 CONTINUE
C
      CALL IZERO(ILOOKUP,600)
      DO 520 ITYPE=NBUCK+1,2*NBUCK
        CALL WHENEQ(ILNBUF,IWHEREK1,1,ITYPE,ILOOKUP,NVAL)
        IOFF=1
        INBUCK=NINBCK(ITYPE)
  521   IDUMP=MIN(NVAL,NBKINT-INBUCK)
        CALL GATHER(IDUMP,BUCK(INBUCK+1,ITYPE),BUFK1,ILOOKUP(IOFF))
        CALL IGATHER(IDUMP,IBUCK(INBUCK+1,ITYPE),IBUFK1,ILOOKUP(IOFF))
        INBUCK=INBUCK+IDUMP
        IF(INBUCK.EQ.NBKINT) THEN
          CALL PLUNK(LUSRT2,BUCK(1,ITYPE),IBUCK(1,ITYPE),ICHAIN(ITYPE),
     &               INBUCK,NBKINT,NREC2)
          NVAL=NVAL-IDUMP
          IOFF=IOFF+IDUMP
          IF(NVAL.NE.0) GOTO 521
        ENDIF
        NINBCK(ITYPE)=INBUCK
  520 CONTINUE
C
      IF(NUT.EQ.ILNBUF) GOTO 3
C 
      NAOINT=NAOINT+(NAOBUF-1)*ILNBUF+NUT
      CLOSE(LUINT,STATUS='KEEP')
C
  200 WRITE(LUOUT,5000)NAOINT
 5000 FORMAT(T3,'@PKSORT22-I, ',I8,' AO integrals sorted as PK ',
     &          'list elements.')
C
C  Flush remaining buffers
C
      DO 50 I=1,NBUCK
        IBK2=NBUCK+I
        CALL PLUNK(LUSRT2,BUCK(1,IBK2),IBUCK(1,IBK2),ICHAIN(IBK2),
     &             NINBCK(IBK2),NBKINT,NREC2)
   50 CONTINUE
C
      CALL TIMER(1)
C
      NREC2=NREC2-1
      WRITE(LUOUT,5010)NREC2,TIMENEW
 5010 FORMAT(T3,'@PKSORT22-I, Wrote out ',I8,' records of K ',
     &          'integrals.',/,
     &       T14,'Integral sort required ',F10.3,' seconds.',/)
C
      RETURN
      END
