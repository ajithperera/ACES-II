      SUBROUTINE DHF23(AOINT,BUF,IBUF,ILNBUF,IRREP1,IRREP2,
     &                 IRREP3,IRREP4,NMO2,NMO3,NMO4,
     &                 NOCC,NSTART1,NEND1,IOFF1,IOFF2,
     &                 IREORD,RHF,ISPIN)
C
CEND
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,AND,OR
      LOGICAL LKVIRT,RHF
C
      DIMENSION AOINT(1),BUF(ILNBUF),IBUF(ILNBUF),NMO(8)
      DIMENSION IOFF1(8),IOFf2(8),NOCC(8),IREORD(1)
C
      COMMON/INTTOL/TOL
      COMMON/HF2FIL/LUHF2,LUHFA,LUHFB
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NDUMMY,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/VTINFO/NPASS1,NPASS2,NPASS3,NPASS4,
     &              NLOAD1,NLOAD2,NLOAD3,NLOAD4,
     &              NWRIT1,NWRIT2,NWRIT3,NWRIT4,
     &              NWRIT1A,NWRIT2A,NWRIT3A,NWRIT4A,
     &              NWRIT1B,NWRIT2B,NWRIT3B,NWRIT4B
C
      IPACK(I,J,K,L)=OR(OR(OR(I,ISHFT(J,IBITWD)),ISHFT(K,2*IBITWD)),
     &               ISHFT(L,3*IBITWD))
C
C LOOP OVER ALL INTEGRALS
C
      IF(RHF.OR.ISPIN.LE.2) THEN
C
       IF(RHF) THEN
        LHF2=LUHF2
       ELSE
        IF(ISPIN.EQ.1) THEN
         LHF2=LUHFA
        ELSE IF(ISPIN.EQ.2) THEN
         LHF2=LUHFB
        ENDIF
       ENDIF
C
       IND=0
       ICNT=0
       NINT=0
       NOCC2=NOCC(IRREP2)
       NOCC3=NOCC(IRREP3)
       NOCC4=NOCC(IRREP4)
C 
       DO 100 I=NSTART1,NEND1
C
        IX=IREORD(I+IOFF1(IRREP1))
C
        DO 200 J=1,NMO2
C
         JX=IREORD(J+IOFF1(IRREP2))
C
         DO 300 K=1,NMO3
C
          KX=IREORD(K+IOFF1(IRREP3))
C
          DO 400 L=1,NMO4
C
           LX=IREORD(L+IOFF1(IRREP4))
C
           LKVIRT=(L.GT.NOCC4).AND.(K.GT.NOCC3)
           IF(KX.GT.LX) THEN
            KXX=KX
            LXX=LX
           ELSE
            KXX=LX
            LXX=KX
           ENDIF
C
C  WRITE OUT ONLY CANONICAL INTEGRALS
C         
           IND=IND+1
           IF(IX.GT.JX) GO TO 400
          
           IF((JX.EQ.KXX).AND.(IX.LT.LXX).AND.(.NOT.LKVIRT)) 
     &        GO TO 400
           IF(JX.GT.KXX.AND.(.NOT.LKVIRT)) GO TO 400
           IF(JX.GT.KXX) THEN
            JXT=KXX
            KXT=JX
            IXT=LXX
            LXT=IX 
           ELSE
            JXT=JX
            KXT=KXX
            IXT=IX
            LXT=LXX 
           ENDIF
           VALUE=AOINT(IND)
           IF(ABS(VALUE).GT.TOL) THEN
C
            ICNT=ICNT+1
            BUF(ICNT)=VALUE
            IBUF(ICNT)=IPACK(IXT,JXT,LXT,KXT)
C
            IF(ICNT.EQ.ILNBUF) THEN
C
C FLUSH BUFFERS
C
             WRITE(LHF2) BUF,IBUF,ILNBUF
C
             NINT=NINT+ILNBUF
             ICNT=0
C
            ENDIF
C
           ENDIF
C
400       CONTINUE
300      CONTINUE
250      CONTINUE
200     CONTINUE
100    CONTINUE
C
C IF THERE ARE STILL INTEGRALS IN THE BUFFER,
C DUMP THEM TO HF2
C
       IF(ICNT.NE.0) WRITE(LHF2) BUF,IBUF,ICNT
       NINT=NINT+ICNT
C
       IF(RHF) THEN
        NWRIT3=NWRIT3+NINT
       ELSE
        IF(ISPIN.EQ.1) THEN
         NWRIT3A=NWRIT3A+NINT
        ELSE
         NWRIT3B=NWRIT3B+NINT
        ENDIF
       ENDIF
C
      ELSE IF(ISPIN.EQ.3) THEN
C
       IND=0
       ICNT=0
       NINT=0
C 
       DO 1100 I=NSTART1,NEND1
C
        IX=IREORD(I+IOFF1(IRREP1))
C
        DO 1200 J=1,NMO2
C
         JX=IREORD(J+IOFF1(IRREP2))
C
         DO 1300 K=1,NMO3
C
          KX=IREORD(K+IOFF2(IRREP3))
C
          DO 1400 L=1,NMO4
C
           LX=IREORD(L+IOFF2(IRREP4))
C
           IF(KX.GT.LX) THEN
            KXX=KX
            LXX=LX
           ELSE
            KXX=LX
            LXX=KX
           ENDIF
C
C  WRITE OUT ONLY CANONICAL INTEGRALS
C         
           IND=IND+1
           IF(IX.GT.JX) GO TO 1400
C          
           VALUE=AOINT(IND)
           IF(ABS(VALUE).GT.TOL) THEN
C
            ICNT=ICNT+1
            BUF(ICNT)=VALUE
            IBUF(ICNT)=IPACK(LXX,KXX,IX,JX)
C
            IF(ICNT.EQ.ILNBUF) THEN
C
C FLUSH BUFFERS
C
             WRITE(LUHF2) BUF,IBUF,ILNBUF
C
             NINT=NINT+ILNBUF
             ICNT=0
C
            ENDIF
C
           ENDIF
C
1400      CONTINUE
1300     CONTINUE
1250     CONTINUE
1200    CONTINUE
1100   CONTINUE
C
C IF THERE ARE STILL INTEGRALS IN THE BUFFER,
C DUMP THEM TO HF2
C
       IF(ICNT.NE.0) WRITE(LUHF2) BUF,IBUF,ICNT
       NINT=NINT+ICNT
C
       NWRIT3=NWRIT3+NINT
C
      ELSE IF(ISPIN.EQ.4) THEN
C
       IND=0
       ICNT=0
       NINT=0
       NOCC3=NOCC(IRREP3)
       NOCC4=NOCC(IRREP4)
C 
       DO 2100 I=NSTART1,NEND1
C
        IX=IREORD(I+IOFF1(IRREP1))
C
        DO 2200 J=1,NMO2
C
         JX=IREORD(J+IOFF1(IRREP2))
C
         DO 2300 K=1,NMO3
C
          KX=IREORD(K+IOFF2(IRREP3))
C
          DO 2400 L=1,NMO4
C
           LX=IREORD(L+IOFF2(IRREP4))
C
           IF(KX.GT.LX) THEN
            KXX=KX
            LXX=LX
           ELSE
            KXX=LX
            LXX=KX
           ENDIF
C
C  WRITE OUT ONLY CANONICAL INTEGRALS
C         
           IND=IND+1
           IF(L.LE.NOCC4.OR.K.LE.NOCC3) GO TO 2400
           IF(IX.GT.JX) GO TO 2400
C          
           VALUE=AOINT(IND)
           IF(ABS(VALUE).GT.TOL) THEN
C
            ICNT=ICNT+1
            BUF(ICNT)=VALUE
            IBUF(ICNT)=IPACK(IX,JX,LXX,KXX)
C
            IF(ICNT.EQ.ILNBUF) THEN
C
C FLUSH BUFFERS
C
             WRITE(LUHF2) BUF,IBUF,ILNBUF
C
             NINT=NINT+ILNBUF
             ICNT=0
C
            ENDIF
C
           ENDIF
C
2400      CONTINUE
2300     CONTINUE
2250     CONTINUE
2200    CONTINUE
2100   CONTINUE
C
C IF THERE ARE STILL INTEGRALS IN THE BUFFER,
C DUMP THEM TO HF2
C
       IF(ICNT.NE.0) WRITE(LUHF2) BUF,IBUF,ICNT
       NINT=NINT+ICNT
C
       NWRIT3=NWRIT3+NINT
C
      ENDIF
      RETURN
      END 
