
        

      
      SUBROUTINE ONETRIP(T2,Z,BUF,IBUF,IK0,IL0,JK0,JL0,KI0,KJ0,LI0,LJ0,
     &                   ITYPE,IWHERE,ILNBUF,NAO,NOCC1,NOCC2,LENIJ,
     &                   LUINT)
C
C THIS ROUTINE LOADS THE AO INTEGRALS FROM THE IIII FILE. 
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR
      DIMENSION BUF(ILNBUF),IBUF(ILNBUF)
      DIMENSION Z(NAO*NAO,*),T2(NAO*NAO,LENIJ)
      DIMENSION IK0(ILNBUF),IL0(ILNBUF),JK0(ILNBUF),JL0(ILNBUF)
      DIMENSION KI0(ILNBUF),KJ0(ILNBUF),LI0(ILNBUF),LJ0(ILNBUF)
      DIMENSION ITYPE(ILNBUF),IWHERE(ILNBUF)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      IUPKI(INT)=AND(INT,IALONE)
      IUPKJ(INT)=AND(ISHFT(INT,-IBITWD),IALONE)
      IUPKK(INT)=AND(ISHFT(INT,-2*IBITWD),IALONE)
      IUPKL(INT)=AND(ISHFT(INT,-3*IBITWD),IALONE)
      IPACK(I,J,K,L)=OR(OR(OR(I,ISHFT(J,IBITWD)),ISHFT(K,2*IBITWD)),
     &                  ISHFT(L,3*IBITWD))
      INDX(I,J)=J+(I*(I-1))/2
      INDX3(I,J)=I+(J*(J-1))/2
      INDX2(I,J,N)=I+(J-1)*N
C
      NAOBUF=0
      NUMINT=0
      NAO2=NAO*NAO
      CALL LOCATE(LUINT,'TWOELSUP')
      CALL ZERO(Z,NAO*NAO*LENIJ)
1     READ(LUINT)BUF,IBUF,NUT
      NAOBUF=NAOBUF+1      
      DO 10 INT=1,NUT
       X=BUF(INT)
       INDI=IUPKI(IBUF(INT))
       INDJ=IUPKJ(IBUF(INT))
       INDK=IUPKK(IBUF(INT))
       INDL=IUPKL(IBUF(INT))
       IJ=INDX(INDJ,INDI)
       KL=INDX(INDL,INDK)
       IK0(INT)=INDX2(INDI,INDK,NAO)
       IL0(INT)=INDX2(INDI,INDL,NAO)
       JK0(INT)=INDX2(INDJ,INDK,NAO)
       JL0(INT)=INDX2(INDJ,INDL,NAO)
       KI0(INT)=INDX2(INDK,INDI,NAO)
       KJ0(INT)=INDX2(INDK,INDJ,NAO)
       LI0(INT)=INDX2(INDL,INDI,NAO)
       LJ0(INT)=INDX2(INDL,INDJ,NAO)
       IF(INDI.NE.INDJ.AND.INDK.NE.INDL.AND.IJ.NE.KL)THEN
        ITYPE(INT)=1
       ELSEIF(INDI.NE.INDJ.AND.INDK.NE.INDL.AND.IJ.EQ.KL)THEN
        ITYPE(INT)=2
       ELSEIF(INDI.EQ.INDJ.AND.INDK.NE.INDL)THEN
        ITYPE(INT)=3
       ELSEIF(INDI.NE.INDJ.AND.INDK.EQ.INDL)THEN
        ITYPE(INT)=4
       ELSEIF(INDI.EQ.INDJ.AND.INDK.EQ.INDL.AND.IJ.EQ.KL)THEN
        ITYPE(INT)=5
       ELSEIF(INDI.EQ.INDJ.AND.INDK.EQ.INDL.AND.IJ.NE.KL)THEN
        ITYPE(INT)=6
       ENDIF
10    CONTINUE
C
C PROCESS TYPE 1 INTEGRALS
C
      CALL WHENEQ(NUT,ITYPE,1,1,IWHERE,NMATCH)
      DO 20 INT=1,NMATCH
       IOFF=IWHERE(INT)
       IK=IK0(IOFF)
       IL=IL0(IOFF)
       JK=JK0(IOFF)
       JL=JL0(IOFF)
       KI=KI0(IOFF)
       KJ=KJ0(IOFF)
       LI=LI0(IOFF)
       LJ=LJ0(IOFF)
       X=BUF(IOFF)
       DO 120 I=1,LENIJ
        Z(IK,I)=X*T2(JL,I)+Z(IK,I)
        Z(IL,I)=X*T2(JK,I)+Z(IL,I)
        Z(JK,I)=X*T2(IL,I)+Z(JK,I)
        Z(JL,I)=X*T2(IK,I)+Z(JL,I)
        Z(KI,I)=X*T2(LJ,I)+Z(KI,I)
        Z(KJ,I)=X*T2(LI,I)+Z(KJ,I)
        Z(LI,I)=X*T2(KJ,I)+Z(LI,I)
        Z(LJ,I)=X*T2(KI,I)+Z(LJ,I)
120    CONTINUE
20    CONTINUE
C
C PROCESS TYPE 2 INTEGRALS
C
      CALL WHENEQ(NUT,ITYPE,1,2,IWHERE,NMATCH)
      DO 30 INT=1,NMATCH
       IOFF=IWHERE(INT)
       IK=IK0(IOFF)
       IL=IL0(IOFF)
       JK=JK0(IOFF)
       JL=JL0(IOFF)
       KI=KI0(IOFF)
       KJ=KJ0(IOFF)
       LI=LI0(IOFF)
       LJ=LJ0(IOFF)
       X=BUF(IOFF)
       DO 130 I=1,LENIJ
        Z(IK,I)=X*T2(JL,I)+Z(IK,I)
        Z(IL,I)=X*T2(JK,I)+Z(IL,I)
        Z(JK,I)=X*T2(IL,I)+Z(JK,I)
        Z(JL,I)=X*T2(IK,I)+Z(JL,I)
130    CONTINUE
30    CONTINUE
C
C PROCESS TYPE 3 INTEGRALS
C
      CALL WHENEQ(NUT,ITYPE,1,3,IWHERE,NMATCH)
      DO 40 INT=1,NMATCH
       IOFF=IWHERE(INT)
       IK=IK0(IOFF)
       IL=IL0(IOFF)
       JK=JK0(IOFF)
       JL=JL0(IOFF)
       KI=KI0(IOFF)
       KJ=KJ0(IOFF)
       LI=LI0(IOFF)
       LJ=LJ0(IOFF)
       X=BUF(IOFF)
       DO 140 I=1,LENIJ
c       DO 140 I=1,0
        Z(IK,I)=X*T2(JL,I)+Z(IK,I)
        Z(IL,I)=X*T2(JK,I)+Z(IL,I)
        Z(KI,I)=X*T2(LJ,I)+Z(KI,I)
        Z(LI,I)=X*T2(KJ,I)+Z(LI,I)
140    CONTINUE
40    CONTINUE
C
C PROCESS TYPE 4 INTEGRALS
C
      CALL WHENEQ(NUT,ITYPE,1,4,IWHERE,NMATCH)
      DO 50 INT=1,NMATCH
       IOFF=IWHERE(INT)
       IK=IK0(IOFF)
       IL=IL0(IOFF)
       JK=JK0(IOFF)
       JL=JL0(IOFF)
       KI=KI0(IOFF)
       KJ=KJ0(IOFF)
       LI=LI0(IOFF)
       LJ=LJ0(IOFF)
       X=BUF(IOFF)
       DO 150 I=1,LENIJ
c        Z(IK,I)=X*T2(JL,I)+Z(IK,I)
c        Z(IL,I)=X*T2(JK,I)+Z(IL,I)
c        Z(JK,I)=X*T2(IL,I)+Z(JK,I)
c        Z(JL,I)=X*T2(IK,I)+Z(JL,I)
        Z(IK,I)=X*T2(JL,I)+Z(IK,I)
        Z(JK,I)=X*T2(IL,I)+Z(JK,I)
        Z(KI,I)=X*T2(LJ,I)+Z(KI,I)
        Z(KJ,I)=X*T2(LI,I)+Z(KJ,I)
150    CONTINUE
50    CONTINUE
C
C
C PROCESS TYPE 5 INTEGRALS
C
      CALL WHENEQ(NUT,ITYPE,1,5,IWHERE,NMATCH)
      DO 60 INT=1,NMATCH
       IOFF=IWHERE(INT)
       IK=IK0(IOFF)
       IL=IL0(IOFF)
       JK=JK0(IOFF)
       JL=JL0(IOFF)
       KI=KI0(IOFF)
       KJ=KJ0(IOFF)
       LI=LI0(IOFF)
       LJ=LJ0(IOFF)
       X=BUF(IOFF)
       DO 160 I=1,LENIJ
        Z(IK,I)=X*T2(JL,I)+Z(IK,I)
160    CONTINUE
60    CONTINUE
C
C PROCESS TYPE 6 INTEGRALS
C
      CALL WHENEQ(NUT,ITYPE,1,6,IWHERE,NMATCH)
      DO 70 INT=1,NMATCH
       IOFF=IWHERE(INT)
       IK=IK0(IOFF)
       IL=IL0(IOFF)
       JK=JK0(IOFF)
       JL=JL0(IOFF)
       KI=KI0(IOFF)
       KJ=KJ0(IOFF)
       LI=LI0(IOFF)
       LJ=LJ0(IOFF)
       X=BUF(IOFF)
       DO 170 I=1,LENIJ
        Z(IK,I)=X*T2(JL,I)+Z(IK,I)
        Z(KI,I)=X*T2(LJ,I)+Z(KI,I)
170    CONTINUE
70    CONTINUE
C
      NUMINT=NUMINT+NUT
C
      IF(NUT.NE.-1)GOTO 1
C
c      WRITE(6,*)' processed ',numint,' ao basis integrals ',
c     &          'from ',naobuf,' buffers.'
C
      RETURN
      END
