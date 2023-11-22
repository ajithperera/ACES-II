      SUBROUTINE DHF11(AOINT,BUF,IBUF,ILNBUF,IRREP,NSTART,NEND,
     &                 NMO,IOFF1,IOFF2,IREORD,RHF,ISPIN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,AND,OR
      LOGICAL RHF
C
      DIMENSION AOINT(1),BUF(ILNBUF),IBUF(ILNBUF)
      DIMENSION IREORD(1)
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
      IF(RHF.OR.ISPIN.NE.3) THEN
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
#ifdef _DEBUG_LVLM1
       Write(6,*) 
       Int_length=(NEND-NSTART+1)*NMO*NMO
       Write(6,*) "@-DHF11 transformed ints before write",Int_length
       Write(*,'(4(1X,F12.7)))'), (AOint(I), I=1, Int_length)
#endif
        
       IND=0
       ICNT=0
       NINT=0
       DO 100 I=NSTART,NEND
C
        IX=IREORD(I+IOFF1)
C
        DO 200 J=1,NMO
C
         JX=IREORD(J+IOFF1)
C      
         DO 300 K=1,NMO
C
          KX=IREORD(K+IOFF1)
C
          DO 400 L=1,K
C
           LX=IREORD(L+IOFF1)
C
C  WRITE OUT ONLY CANONICAL INTEGRALS
C         
           IND=IND+1
           IF(IX.LT.JX) GO TO 400
           IF(IX.LT.KX) GO TO 400
           IF(IX.EQ.KX.AND.JX.LT.LX) GO TO 400
           VALUE=AOINT(IND)
#ifdef _DEBUG_LVM1
C         write(*,"(t3,a,1x,I3,1x,i3,1x,i3,1x,i3)"),
C     &         "IX,JX,KX,LX ", IX,JX, KX, LX
         Write(6,"(a, I5, F10.5)") "The transformed Int",IND, VALUE
#endif
           IF(ABS(VALUE).GT.TOL) THEN
C
            ICNT=ICNT+1
            BUF(ICNT)=VALUE
            IBUF(ICNT)=IPACK(LX,KX,JX,IX)
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
200     CONTINUE
100    CONTINUE
C
C IF THERE ARE STILL INTEGRALS IN THE BUFFER,
C DUMP THEM TO HF2
C
#ifdef _DEBUG_LVM1
       Write(6,*)
       Write(6,"(a,I5)") "The number of ints", icnt
#endif

       IF(ICNT.NE.0) WRITE(LHF2) BUF,IBUF,ICNT
       NINT=NINT+ICNT
C
       IF(RHF) THEN
        NWRIT1=NWRIT1+NINT
       ELSE
        IF(ISPIN.EQ.1) THEN
         NWRIT1A=NWRIT1A+NINT
        ELSE
         NWRIT1B=NWRIT1B+NINT
        ENDIF
       ENDIF
C
      ELSE
C
       IND=0
       ICNT=0
       NINT=0
       DO 1100 I=NSTART,NEND
C
        IX=IREORD(I+IOFF1)
C
        DO 1200 J=1,NMO
C
         JX=IREORD(J+IOFF1)
C      
         DO 1300 K=1,NMO
C
          KX=IREORD(K+IOFF2)
C
          DO 1400 L=1,K
C
           LX=IREORD(L+IOFF2)
C
C  WRITE OUT ONLY CANONICAL INTEGRALS
C         
           IND=IND+1
           IF(IX.LT.JX) GO TO 1400
           IF(KX.LT.LX) GO TO 1400
#ifdef _DEBUG_LVLM1
         write(*,"(t3,a,1x,I3,1x,i3,1x,i3,1x,i3)"),
     &         "IX,JX,KX,LX ", IX,JX, KX, LX
         Write(6,"(a, F10.5)") "The transformed Int",VALUE
#endif

           VALUE=AOINT(IND)
           IF(ABS(VALUE).GT.TOL) THEN
C
            ICNT=ICNT+1
            BUF(ICNT)=VALUE
            IBUF(ICNT)=IPACK(LX,KX,JX,IX)
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
1200    CONTINUE
1100   CONTINUE
C
C IF THERE ARE STILL INTEGRALS IN THE BUFFER,
C DUMP THEM TO HF2
C
#ifdef _DEBUG_LVM1
       Write(6,*) 
       Write(6,"(a,I5)") "The number of ints", icnt
#endif
       IF(ICNT.NE.0) WRITE(LUHF2) BUF,IBUF,ICNT
       NINT=NINT+ICNT
C
       NWRIT1=NWRIT1+NINT
C
      ENDIF
      RETURN
      END 
