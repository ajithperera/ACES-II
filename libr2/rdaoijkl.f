      SUBROUTINE RDAOIJKL(T2,Z,BUF,IBUF,IK0,IL0,JK0,JL0,KI0,KJ0,LI0,LJ0,
     &                    ITYPE,IWHERE,ISYM,IAOSYM,IMAP,
     &                    NUMDIS,ILNBUF,LUINT,ISPIN,NAO,IMAX)
C    &
C THIS ROUTINE LOADS THE AO INTEGRALS FROM THE IJKL FILE (ALL FOUR
C INDICES HAVE DIFFERENT SYMMETRIES) AND CONTRACTS THEM WITH THE T2
C AMPLITUDES.  SOME COMPLICATED STUFF!
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,DIRPRD
      DIMENSION BUF(ILNBUF),IBUF(ILNBUF)
      DIMENSION Z(IMAX),T2(IMAX),IMAP(*)
      DIMENSION IK0(ILNBUF),IL0(ILNBUF),JK0(ILNBUF),JL0(ILNBUF)
      DIMENSION KI0(ILNBUF),KJ0(ILNBUF),LI0(ILNBUF),LJ0(ILNBUF)
      DIMENSION ITYPE(ILNBUF),IWHERE(ILNBUF),IAOSYM(*),ISYM(ILNBUF,4)
      DIMENSION IOFFMO(8),NUMDIS(8)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),IOFFV(8,2),IOFFO(8,2),
     &             IRPDPDAO(8),IRPDPDAOS(8),ISTART(8,8),ISTARTMO(8,3)
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
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
      CALL LOCATE(LUINT,'TWOELSUP')
1     READ(LUINT)BUF,IBUF,NUT
      NAOBUF=NAOBUF+1      
      CALL ICOPY(8,ISTARTMO(1,ISPIN),1,IOFFMO,1)
      DO 10 INT=1,NUT
C
       X=BUF(INT)
C
       INDI=IUPKI(IBUF(INT))
       INDJ=IUPKJ(IBUF(INT))
       INDK=IUPKK(IBUF(INT))
       INDL=IUPKL(IBUF(INT))
       IJ=INDX(MAX(INDJ,INDI),MIN(INDJ,INDI))
       KL=INDX(MAX(INDL,INDK),MIN(INDL,INDK))
C
       ISYM(INT,1)=IAOSYM(INDI)
       ISYM(INT,2)=IAOSYM(INDJ)
       ISYM(INT,3)=IAOSYM(INDK)
       ISYM(INT,4)=IAOSYM(INDL)
C
       IK0(INT)=IMAP(INDX2(INDI,INDK,NAO))
       IL0(INT)=IMAP(INDX2(INDI,INDL,NAO))
       JK0(INT)=IMAP(INDX2(INDJ,INDK,NAO))
       JL0(INT)=IMAP(INDX2(INDJ,INDL,NAO))
       KI0(INT)=IMAP(INDX2(INDK,INDI,NAO))
       KJ0(INT)=IMAP(INDX2(INDK,INDJ,NAO))
       LI0(INT)=IMAP(INDX2(INDL,INDI,NAO))
       LJ0(INT)=IMAP(INDX2(INDL,INDJ,NAO))
C
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
       X=BUF(IOFF)
C
       IK=IK0(IOFF)
       IL=IL0(IOFF)
       JK=JK0(IOFF)
       JL=JL0(IOFF)
       KI=KI0(IOFF)
       KJ=KJ0(IOFF)
       LI=LI0(IOFF)
       LJ=LJ0(IOFF)
C
       IRRI=ISYM(IOFF,1)
       IRRJ=ISYM(IOFF,2)
       IRRK=ISYM(IOFF,3)
       IRRL=ISYM(IOFF,4)
       IRRIK=DIRPRD(IRRI,IRRK)
       IRRIL=DIRPRD(IRRI,IRRL)
       IKOFF=IOFFMO(IRRIK)
       ILOFF=IOFFMO(IRRIL)
       IKDSZ=IRPDPDAO(IRRIK)
       ILDSZ=IRPDPDAO(IRRIL)
C
       LENIK=NUMDIS(IRRIK)
       LENIL=NUMDIS(IRRIL)
       DO 120 I=1,LENIK
        Z(IK+IKOFF)=X*T2(JL+IKOFF)+Z(IK+IKOFF)
        Z(KI+IKOFF)=X*T2(LJ+IKOFF)+Z(KI+IKOFF)
        Z(JL+IKOFF)=X*T2(IK+IKOFF)+Z(JL+IKOFF)
        Z(LJ+IKOFF)=X*T2(KI+IKOFF)+Z(LJ+IKOFF)
        IKOFF=IKOFF+IKDSZ
120    CONTINUE
       DO 121 J=1,LENIL
        Z(IL+ILOFF)=X*T2(JK+ILOFF)+Z(IL+ILOFF)
        Z(JK+ILOFF)=X*T2(IL+ILOFF)+Z(JK+ILOFF)
        Z(KJ+ILOFF)=X*T2(LI+ILOFF)+Z(KJ+ILOFF)
        Z(LI+ILOFF)=X*T2(KJ+ILOFF)+Z(LI+ILOFF)
        ILOFF=ILOFF+ILDSZ
121    CONTINUE
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
C
       IRRI=ISYM(IOFF,1)
       IRRJ=ISYM(IOFF,2)
       IRRK=ISYM(IOFF,3)
       IRRL=ISYM(IOFF,4)
       IRRIK=DIRPRD(IRRI,IRRK)
       IRRIL=DIRPRD(IRRI,IRRL)
       IKOFF=IOFFMO(IRRIK)
       ILOFF=IOFFMO(IRRIL)
       IKDSZ=IRPDPDAO(IRRIK)
       ILDSZ=IRPDPDAO(IRRIL)
C
       LENIK=NUMDIS(IRRIK)
       LENIL=NUMDIS(IRRIL)
C
       DO 130 I=1,LENIK
        Z(IK+IKOFF)=X*T2(JL+IKOFF)+Z(IK+IKOFF)
        Z(JL+IKOFF)=X*T2(IK+IKOFF)+Z(JL+IKOFF)
        IKOFF=IKOFF+IKDSZ
130    CONTINUE
       DO 131 I=1,LENIL      
        Z(IL+ILOFF)=X*T2(JK+ILOFF)+Z(IL+ILOFF)
        Z(JK+ILOFF)=X*T2(IL+ILOFF)+Z(JK+ILOFF)
        ILOFF=ILOFF+ILDSZ
131    CONTINUE
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
C
       IRRI=ISYM(IOFF,1)
       IRRJ=ISYM(IOFF,2)
       IRRK=ISYM(IOFF,3)
       IRRL=ISYM(IOFF,4)
       IRRIK=DIRPRD(IRRI,IRRK)
       IRRIL=DIRPRD(IRRI,IRRL)
       IKOFF=IOFFMO(IRRIK)
       ILOFF=IOFFMO(IRRIL)
       IKDSZ=IRPDPDAO(IRRIK)
       ILDSZ=IRPDPDAO(IRRIL)
C
       LENIK=NUMDIS(IRRIK)
       LENIL=NUMDIS(IRRIL)
       DO 140 I=1,LENIK
        Z(IK+IKOFF)=X*T2(JL+IKOFF)+Z(IK+IKOFF)
        Z(KI+IKOFF)=X*T2(LJ+IKOFF)+Z(KI+IKOFF)
        IKOFF=IKOFF+IKDSZ
140    CONTINUE
       DO 141 J=1,LENIL
        Z(IL+ILOFF)=X*T2(JK+ILOFF)+Z(IL+ILOFF)
        Z(LI+ILOFF)=X*T2(KJ+ILOFF)+Z(LI+ILOFF)
        ILOFF=ILOFF+ILDSZ
141    CONTINUE
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
C
       IRRI=ISYM(IOFF,1)
       IRRJ=ISYM(IOFF,2)
       IRRK=ISYM(IOFF,3)
       IRRL=ISYM(IOFF,4)
       IRRIK=DIRPRD(IRRI,IRRK)
       IRRIL=DIRPRD(IRRI,IRRL)
       IKOFF=IOFFMO(IRRIK)
       ILOFF=IOFFMO(IRRIL)
       IKDSZ=IRPDPDAO(IRRIK)
       ILDSZ=IRPDPDAO(IRRIL)
C
       LENIK=NUMDIS(IRRIK)
       LENIL=NUMDIS(IRRIL)
       DO 150 I=1,LENIK
        Z(IK+IKOFF)=X*T2(JL+IKOFF)+Z(IK+IKOFF)
        Z(KI+IKOFF)=X*T2(LJ+IKOFF)+Z(KI+IKOFF)
        Z(JL+IKOFF)=X*T2(IK+IKOFF)+Z(JL+IKOFF)
        Z(LJ+IKOFF)=X*T2(KI+IKOFF)+Z(LJ+IKOFF)
        IKOFF=IKOFF+IKDSZ
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
C
       IRRI=ISYM(IOFF,1)
       IRRJ=ISYM(IOFF,2)
       IRRK=ISYM(IOFF,3)
       IRRL=ISYM(IOFF,4)
       IRRIK=DIRPRD(IRRI,IRRK)
       IRRIL=DIRPRD(IRRI,IRRL)
       IKOFF=IOFFMO(IRRIK)
       ILOFF=IOFFMO(IRRIL)
       IKDSZ=IRPDPDAO(IRRIK)
       ILDSZ=IRPDPDAO(IRRIL)
C
       LENIK=NUMDIS(IRRIK)
       LENIL=NUMDIS(IRRIL)

       DO 160 I=1,LENIK
        Z(IK+IKOFF)=X*T2(JL+IKOFF)+Z(IK+IKOFF)
        IKOFF=IKOFF+IKDSZ
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
C
       IRRI=ISYM(IOFF,1)
       IRRJ=ISYM(IOFF,2)
       IRRK=ISYM(IOFF,3)
       IRRL=ISYM(IOFF,4)
       IRRIK=DIRPRD(IRRI,IRRK)
       IRRIL=DIRPRD(IRRI,IRRL)
       IKOFF=IOFFMO(IRRIK)
       ILOFF=IOFFMO(IRRIL)
       IKDSZ=IRPDPDAO(IRRIK)
       ILDSZ=IRPDPDAO(IRRIL)
C
       LENIK=NUMDIS(IRRIK)
       LENIL=NUMDIS(IRRIL)
C
       DO 170 I=1,LENIK
        Z(IK+IKOFF)=X*T2(JL+IKOFF)+Z(IK+IKOFF)
        Z(KI+IKOFF)=X*T2(LJ+IKOFF)+Z(KI+IKOFF)
        IKOFF=IKOFF+IKDSZ
170    CONTINUE
70    CONTINUE
C
      NUMINT=NUMINT+MAX(NUT,0)
C
      IF(NUT.NE.-1)GOTO 1
C
c      WRITE(6,*)' processed ',numint,' ao basis integrals ',
c     &          'from ',naobuf-1,' buffers.'
C
      RETURN
      END
