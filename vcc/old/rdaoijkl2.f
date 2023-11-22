C
C THIS ROUTINE LOADS THE AO INTEGRALS FROM THE IJKL FILE (ALL FOUR
C INDICES HAVE DIFFERENT SYMMETRIES) AND CONTRACTS THEM WITH THE T2
C AMPLITUDES.  SOME COMPLICATED STUFF!
C
      SUBROUTINE RDAOIJKL2(T2,Z,BUF,IBUF,IK0,IL0,JK0,JL0,KI0,KJ0,
     &                    VALUE,ISYM,ITYPE,IWHERE,IAOSYM,IMAP,
     &                    ILNBUF,LUINT,IUHF,NAO,IMAX,IRREPX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,DIRPRD
      DIMENSION BUF(ILNBUF),IBUF(ILNBUF)
      DIMENSION Z(IMAX),T2(IMAX),IMAP(*)
      DIMENSION IK0(8*ILNBUF),IL0(8*ILNBUF),JK0(8*ILNBUF)
      DIMENSION KI0(8*ILNBUF),KJ0(8*ILNBUF),JL0(8*ILNBUF)
      DIMENSION VALUE(8*ILNBUF),ISYM(8*ILNBUF)
      DIMENSION ITYPE(8*ILNBUF),IWHERE(8*ILNBUF),IAOSYM(*)
      DIMENSION IOFFMO(8)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),IOFFV(8,2),IOFFO(8,2),
     &             IRPDPDAO(8),IRPDPDAOS(8),ISTART(8,8),ISTARTMO(8,3)
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/FLAGS/IFLAGS(100)
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW     
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
      DATA HALF/0.5D0/
C
      TOTAL=0.0D0
      CALL TIMER(1)
      ISTART1=TIMENOW
      NAOBUF=0
      NUMINT=0
      CALL LOCATE(LUINT,'TWOELSUP')
C
      IF(IUHF.EQ.1) THEN
C
C UHF
C
1      READ(LUINT)BUF,IBUF,NUT
       NUT8=8*NUT
       NAOBUF=NAOBUF+1      
       ISTICK=0
       CALL IZERO(ISYM,NUT8)
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
        IK=INDX2(INDI,INDK,NAO)
        KI=INDX2(INDK,INDI,NAO)
        IL=INDX2(INDI,INDL,NAO)
        LI=INDX2(INDL,INDI,NAO)
        JK=INDX2(INDJ,INDK,NAO)
        KJ=INDX2(INDK,INDJ,NAO)
        JL=INDX2(INDJ,INDL,NAO)
        LJ=INDX2(INDL,INDJ,NAO)
C
        ISYM1=DIRPRD(IAOSYM(INDI),IAOSYM(INDK))
        ISYM2=DIRPRD(IAOSYM(INDI),IAOSYM(INDL))
        IF(INDI.NE.INDJ.AND.INDK.NE.INDL.AND.IJ.NE.KL)THEN
         IK0(ISTICK+1)=IMAP(IK)
         JL0(ISTICK+1)=IMAP(JL)
         ISYM(ISTICK+1)=ISYM1
         VALUE(ISTICK+1)=X
         IK0(ISTICK+2)=IMAP(IL)
         JL0(ISTICK+2)=IMAP(JK)
         ISYM(ISTICK+2)=ISYM2
         VALUE(ISTICK+2)=X
         IK0(ISTICK+3)=IMAP(JK)
         JL0(ISTICK+3)=IMAP(IL)
         ISYM(ISTICK+3)=ISYM2
         VALUE(ISTICK+3)=X
         IK0(ISTICK+4)=IMAP(JL)
         JL0(ISTICK+4)=IMAP(IK)
         ISYM(ISTICK+4)=ISYM1
         VALUE(ISTICK+4)=X
         IK0(ISTICK+5)=IMAP(KI)
         JL0(ISTICK+5)=IMAP(LJ)
         ISYM(ISTICK+5)=ISYM1
         VALUE(ISTICK+5)=X
         IK0(ISTICK+6)=IMAP(LI)
         JL0(ISTICK+6)=IMAP(KJ)
         ISYM(ISTICK+6)=ISYM2
         VALUE(ISTICK+6)=X
         IK0(ISTICK+7)=IMAP(KJ)
         JL0(ISTICK+7)=IMAP(LI)
         ISYM(ISTICK+7)=ISYM2
         VALUE(ISTICK+7)=X
         IK0(ISTICK+8)=IMAP(LJ)
         JL0(ISTICK+8)=IMAP(KI)
         ISYM(ISTICK+8)=ISYM1
         VALUE(ISTICK+8)=X
        ELSEIF(INDI.NE.INDJ.AND.INDK.NE.INDL.AND.IJ.EQ.KL)THEN
         IK0(ISTICK+1)=IMAP(IK)
         JL0(ISTICK+1)=IMAP(JL)
         ISYM(ISTICK+1)=ISYM1
         VALUE(ISTICK+1)=X
         IK0(ISTICK+2)=IMAP(IL)
         JL0(ISTICK+2)=IMAP(JK)
         ISYM(ISTICK+2)=ISYM2
         VALUE(ISTICK+2)=X
         IK0(ISTICK+3)=IMAP(JK)
         JL0(ISTICK+3)=IMAP(IL)
         ISYM(ISTICK+3)=ISYM2
         VALUE(ISTICK+3)=X
         IK0(ISTICK+4)=IMAP(JL)
         JL0(ISTICK+4)=IMAP(IK)
         ISYM(ISTICK+4)=ISYM1
         VALUE(ISTICK+4)=X
        ELSEIF(INDI.EQ.INDJ.AND.INDK.NE.INDL)THEN
         IK0(ISTICK+1)=IMAP(IK)
         JL0(ISTICK+1)=IMAP(JL)
         ISYM(ISTICK+1)=ISYM1
         VALUE(ISTICK+1)=X
         IK0(ISTICK+2)=IMAP(IL)
         JL0(ISTICK+2)=IMAP(JK)
         ISYM(ISTICK+2)=ISYM2
         VALUE(ISTICK+2)=X
         IK0(ISTICK+5)=IMAP(KI)
         JL0(ISTICK+5)=IMAP(LJ)
         ISYM(ISTICK+5)=ISYM1
         VALUE(ISTICK+5)=X
         IK0(ISTICK+6)=IMAP(LI)
         JL0(ISTICK+6)=IMAP(KJ)
         ISYM(ISTICK+6)=ISYM2
         VALUE(ISTICK+6)=X
        ELSEIF(INDI.NE.INDJ.AND.INDK.EQ.INDL)THEN
         IK0(ISTICK+1)=IMAP(IK)
         JL0(ISTICK+1)=IMAP(JL)
         ISYM(ISTICK+1)=ISYM1
         VALUE(ISTICK+1)=X
         IK0(ISTICK+4)=IMAP(JL)
         JL0(ISTICK+4)=IMAP(IK)
         ISYM(ISTICK+4)=ISYM1
         VALUE(ISTICK+4)=X
         IK0(ISTICK+5)=IMAP(KI)
         JL0(ISTICK+5)=IMAP(LJ)
         ISYM(ISTICK+5)=ISYM1
         VALUE(ISTICK+5)=X
         IK0(ISTICK+8)=IMAP(LJ)
         JL0(ISTICK+8)=IMAP(KI)
         ISYM(ISTICK+8)=ISYM1
         VALUE(ISTICK+8)=X
        ELSEIF(INDI.EQ.INDJ.AND.INDK.EQ.INDL.AND.IJ.EQ.KL)THEN
         IK0(ISTICK+1)=IMAP(IK)
         JL0(ISTICK+1)=IMAP(JL)
         ISYM(ISTICK+1)=ISYM1
         VALUE(ISTICK+1)=X
        ELSEIF(INDI.EQ.INDJ.AND.INDK.EQ.INDL.AND.IJ.NE.KL)THEN
         IK0(ISTICK+1)=IMAP(IK)
         JL0(ISTICK+1)=IMAP(JL)
         ISYM(ISTICK+1)=ISYM1
         VALUE(ISTICK+1)=X
         IK0(ISTICK+5)=IMAP(KI)
         JL0(ISTICK+5)=IMAP(LJ)
         ISYM(ISTICK+5)=ISYM1
         VALUE(ISTICK+5)=X
        ENDIF
        ISTICK=ISTICK+8
10     CONTINUE
C
       CALL TIMER(1)
       DO 20 IRREPIJ = 1, NIRREP
          IRREPAB = DIRPRD(IRREPIJ,IRREPX)
          CALL WHENEQ(NUT8,ISYM,1,IRREPAB,IWHERE,NMATCH)
       DO 20 ISPIN = 1, 3
          IKOFF  = ISTARTMO(IRREPIJ,ISPIN)
          LENIK  = IRPDPD(IRREPIJ,ISYTYP(2,43+ISPIN))
          LENIK5 = LENIK/5
          DO 21 INT = 1, NMATCH
             IOFF = IWHERE(INT)
             X  = VALUE(IOFF)
             IK = IK0(IOFF)
             JL = JL0(IOFF)
             IOFF1 = 1 + IKOFF + (JL-1)*LENIK
             IOFF2 = 1 + IKOFF + (IK-1)*LENIK
             IREST = MOD(LENIK,5)
             DO 22 I = 1, LENIK5
                Z(IOFF2)   = Z(IOFF2)   + X*T2(IOFF1)
                Z(IOFF2+1) = Z(IOFF2+1) + X*T2(IOFF1+1)
                Z(IOFF2+2) = Z(IOFF2+2) + X*T2(IOFF1+2)
                Z(IOFF2+3) = Z(IOFF2+3) + X*T2(IOFF1+3)
                Z(IOFF2+4) = Z(IOFF2+4) + X*T2(IOFF1+4)
                IOFF1 = IOFF1 + 5
                IOFF2 = IOFF2 + 5
22           CONTINUE
             DO 23 I = 1, IREST
                Z(IOFF2) = Z(IOFF2) + X*T2(IOFF1)
                IOFF1 = IOFF1 + 1
                IOFF2 = IOFF2 + 1
23           CONTINUE
21        CONTINUE
20     CONTINUE
       CALL TIMER(1)
       TCPU  = TIMENEW
       TOTAL = TOTAL + TIMENEW
C
       NUMINT = NUMINT + MAX(NUT,0)
C
       IF (NUT.NE.-1) GOTO 1
C
      ELSE
C
C RHF
C
101    READ(LUINT)BUF,IBUF,NUT
       NUT4=4*NUT
       NAOBUF=NAOBUF+1      
       ISTICK=0
       CALL IZERO(ISYM,NUT4)
       DO 110 INT=1,NUT
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
        IK=INDX2(INDI,INDK,NAO)
        KI=INDX2(INDK,INDI,NAO)
        IL=INDX2(INDI,INDL,NAO)
        LI=INDX2(INDL,INDI,NAO)
        JK=INDX2(INDJ,INDK,NAO)
        KJ=INDX2(INDK,INDJ,NAO)
        JL=INDX2(INDJ,INDL,NAO)
        LJ=INDX2(INDL,INDJ,NAO)
C
        ISYM1=DIRPRD(IAOSYM(INDI),IAOSYM(INDK))
        ISYM2=DIRPRD(IAOSYM(INDI),IAOSYM(INDL))
        IF(INDI.NE.INDJ.AND.INDK.NE.INDL.AND.IJ.NE.KL)THEN
         IK0(ISTICK+1)=IMAP(IK)
         JL0(ISTICK+1)=IMAP(JL)
         ISYM(ISTICK+1)=ISYM1
         VALUE(ISTICK+1)=X
         IK0(ISTICK+2)=IMAP(IL)
         JL0(ISTICK+2)=IMAP(JK)
         ISYM(ISTICK+2)=ISYM2
         VALUE(ISTICK+2)=X
         IK0(ISTICK+3)=IMAP(JK)
         JL0(ISTICK+3)=IMAP(IL)
         ISYM(ISTICK+3)=ISYM2
         VALUE(ISTICK+3)=X
         IK0(ISTICK+4)=IMAP(JL)
         JL0(ISTICK+4)=IMAP(IK)
         ISYM(ISTICK+4)=ISYM1
         VALUE(ISTICK+4)=X
        ELSEIF(INDI.NE.INDJ.AND.INDK.NE.INDL.AND.IJ.EQ.KL)THEN
         IK0(ISTICK+1)=IMAP(IK)
         JL0(ISTICK+1)=IMAP(JL)
         ISYM(ISTICK+1)=ISYM1
         VALUE(ISTICK+1)=X*HALF
         IK0(ISTICK+2)=IMAP(IL)
         JL0(ISTICK+2)=IMAP(JK)
         ISYM(ISTICK+2)=ISYM2
         VALUE(ISTICK+2)=X
         IK0(ISTICK+4)=IMAP(JL)
         JL0(ISTICK+4)=IMAP(IK)
         ISYM(ISTICK+4)=ISYM1
         VALUE(ISTICK+4)=X*HALF
        ELSEIF(INDI.EQ.INDJ.AND.INDK.NE.INDL)THEN
         IK0(ISTICK+1)=IMAP(IK)
         JL0(ISTICK+1)=IMAP(JL)
         ISYM(ISTICK+1)=ISYM1
         VALUE(ISTICK+1)=X
         IK0(ISTICK+2)=IMAP(IL)
         JL0(ISTICK+2)=IMAP(JK)
         ISYM(ISTICK+2)=ISYM2
         VALUE(ISTICK+2)=X
        ELSEIF(INDI.NE.INDJ.AND.INDK.EQ.INDL)THEN
         IK0(ISTICK+1)=IMAP(IK)
         JL0(ISTICK+1)=IMAP(JL)
         ISYM(ISTICK+1)=ISYM1
         VALUE(ISTICK+1)=X
         IK0(ISTICK+4)=IMAP(JL)
         JL0(ISTICK+4)=IMAP(IK)
         ISYM(ISTICK+4)=ISYM1
         VALUE(ISTICK+4)=X
        ELSEIF(INDI.EQ.INDJ.AND.INDK.EQ.INDL.AND.IJ.EQ.KL)THEN
         IK0(ISTICK+1)=IMAP(IK)
         JL0(ISTICK+1)=IMAP(JL)
         ISYM(ISTICK+1)=ISYM1
         VALUE(ISTICK+1)=X*HALF
        ELSEIF(INDI.EQ.INDJ.AND.INDK.EQ.INDL.AND.IJ.NE.KL)THEN
         IK0(ISTICK+1)=IMAP(IK)
         JL0(ISTICK+1)=IMAP(JL)
         ISYM(ISTICK+1)=ISYM1
         VALUE(ISTICK+1)=X
        ENDIF
        ISTICK=ISTICK+4
110    CONTINUE
C
       CALL TIMER(1)
       DO 120 IRREPIJ=1,NIRREP
        IRREPAB=DIRPRD(IRREPIJ,IRREPX)
        CALL WHENEQ(NUT4,ISYM,1,IRREPAB,IWHERE,NMATCH)
        IKOFF=ISTARTMO(IRREPIJ,3)
        LENIK=IRPDPD(IRREPIJ,ISYTYP(2,46))
        LENIK5=LENIK/5
        DO 121 INT=1,NMATCH
         IOFF=IWHERE(INT)
         X=VALUE(IOFF)
C
         IK=IK0(IOFF)
         JL=JL0(IOFF)
C 
         IOFF1=IKOFF+(JL-1)*LENIK+1
         IOFF2=IKOFF+(IK-1)*LENIK+1
         IREST=MOD(LENIK,5)
         DO 122 I=1,LENIK5
          Z(IOFF2)=Z(IOFF2)+X*T2(IOFF1)
          Z(IOFF2+1)=Z(IOFF2+1)+X*T2(IOFF1+1)
          Z(IOFF2+2)=Z(IOFF2+2)+X*T2(IOFF1+2)
          Z(IOFF2+3)=Z(IOFF2+3)+X*T2(IOFF1+3)
          Z(IOFF2+4)=Z(IOFF2+4)+X*T2(IOFF1+4)
          IOFF1=IOFF1+5
          IOFF2=IOFF2+5
122      CONTINUE
         DO 123 I=1,IREST
          Z(IOFF2)=Z(IOFF2)+X*T2(IOFF1)
          IOFF1=IOFF1+1
          IOFF2=IOFF2+1
123      CONTINUE
C
121     CONTINUE
C
120    CONTINUE
       CALL TIMER(1)
       TCPU=TIMENEW
       TOTAL=TOTAL+TIMENEW
C
       NUMINT=NUMINT+MAX(NUT,0)
C
       IF(NUT.NE.-1)GOTO 101
C
      ENDIF
C
      CALL TIMER(1)
      TEND1=TIMENOW
      TOTAL1=TEND1-TSTART1
      IF(IFLAGS(1).GE.10)THEN
       write(*,*) 'CPU time for contraction : ',total
       write(*,*) 'Total time for ladders : ',total1
       write(*,*)' processed ',numint,' ao basis integrals ',
     &           'from ',naobuf-1,' buffers.'
      ENDIF
C
      RETURN
      END
