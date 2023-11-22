C     
      SUBROUTINE SecMOM(VLIST,NOC,CENTR1,EXPS,CONT,LABEL,NPRIM,NCONTR,
     &                  ITRM,ITRN,CTRN,NCONT,BUF,IBUF,AOINT,LBUF,INDX,
     &                  NGROUP,VS,SCR,MXPRM2)
C
C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C &  Evaluates diamagnetic shielding contribution to the NMR         &
C &  chemical shifts. Coded by Ajith 03/94 following vprop.f style   &
C &  (Not the most efficient or most readable program).              &
C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND
      CHARACTER*8 DUML1,LABELD(3),LABELJA(10,25),JUNKLAB
C
CJDW 1/8/98. Modification for f90. TITLE was previously double precision
C            and filled with characters.
C
      CHARACTER*8 TITLE
C
      COMMON /PRTBUG/ IPRTFS
      COMMON /PRTDAT/ IDOSPH, MAXAQN, NSABF
      COMMON /NUCINDX/ ICT
      COMMON /PARM/ISL2,ISL3,ISL4,ISLX,THRESH,MTYP,ILNMCM,INT,INUI,
     &             NCENTS,ID20,NSYM,NAO(12),NMO(12),JUNK,PTOT(3,10), 
     &             INTERF,JOBIPH,RUNMOS
      COMMON /HIGHL/ LMNVAL(3,84), ANORM(84)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /VARS/ FNU(20),DP(20),FN(15),FD(15),PI,DFTR(11)
      COMMON /FILES/LUOUT,MOINTS
C
      DIMENSION VLIST(NOC,10),EXPS(1),CONT(1),LABEL(1),ITRM(1),
     &          ITRN(256,1),CTRN(256,1),NCONTR(1),NPRIM(1),INDX(1)
      DIMENSION IHOLD(10),IADH(10),VS(MXPRM2,10)
      DIMENSION BUF(LBUF,10),IBUF(IINTFP*LBUF,10),AOINT(1)
      DIMENSION VAL(3, 3),TITLE(10,25),TEMP1(3),TEMP2(3,3)
      DIMENSION CENTR1(3),SCR(1),GBUF(600),IGNDX(600)
      DIMENSION LMN(200),JMN(200),LFLAG(3),MFLAG(3),NFLAG(3)
      DIMENSION A(3),B(3),EP(3),DISTQD(3),AAA(200),BBB(200),
     &          V(10), DISTQA(3), C(3), G(3)
C     
      DATA LFLAG / 1, 0, 0/
      DATA MFLAG / 0, 1, 0/
      DATA NFLAG / 0, 0, 1/
C
      DATA (LABELJA(I, 21), I = 1, 10) /'  DSXX  ', '  DSYX  ',
     &      '  DSYY  ','  DSZX  ', '  DSZY  ', '  DSZZ  ',
     &     4*'        '/
C     
      DATA (  TITLE(I, 21), I = 1, 10) /'  DSXX  ', '  DSYX  ',
     &      '  DSYY  ','  DSZX  ', '  DSZY  ', '  DSZZ  ',
     &     4*'        '/
C     DATA (TITLE(I, 21), I = 1, 10)/6H  DSXX, 6H  DSYX, 6H  DSYY,
C    &      6H  DSZX, 6H  DSZY, 6H  DSZZ, 4*6H      /
C                 
      IBTAND(I,J) = IAND(I,J)
      IBTOR(I,J)  = IOR(I,J)
      IBTXOR(I,J) = IEOR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTNOT(I)   = NOT(I)
C
      LXQ = 0
      THR2 = THRESH/100000.D+00
      IWPB = (LBUF-2)/2
      LENBUF = IINTFP*LBUF
      ZERO = 0.0D+00
      ONE  = 1.00D+00
      TWO  = 2.00D+00
      EIGHT = 8.00D+00
C
C Do sum initializations (No nuclear contribution for this term)
C
      CALL SETRHF
      NT = 6
      DO 100 I = 1, 6
         IHOLD(I) = 0
         IBUF(1,I) = 0
         IBUF(2,I) = -1
100   CONTINUE
C
      IEXH  = 0
      ICONT = 0
      ICOFF = 0
      CALL DZERO(G, 3)
C  
C Loop over the distinct group of exponents
C      
      DO 50 I = 1, NGROUP
C     
         ICNT = AND(IBTSHR(LABEL(I), 10), 2**10-1)
         ITYP = IBTAND(LABEL(I), 2**10-1)
C     
         DO 410 II = 1, 3
            A(II) = VLIST(ICNT, II)
 410     CONTINUE

         IE = NPRIM(I)
         ICE = NCONTR(I)
         JEXH = 0
         JCONT = 0
         JCOFF = 0
C     
         DO 49 J = 1, I
C     
            JCNT = IBTAND(IBTSHR(LABEL(J),10), 2**10 - 1)
            JTYP = IBTAND(LABEL(J), 2**10-1)
C
            DO 411 JJ = 1, 3
               B(JJ) = VLIST(JCNT, JJ)
 411        CONTINUE
C 
            JE = NPRIM(J)
            JCE = NCONTR(J)
C
            DO 2 M = 1, NT
               DO 28 MM = 1, MXPRM2
                  VS(MM, M) = ZERO
 28            CONTINUE
 2          CONTINUE
C     
            IEX = IEXH
C
C Loop over the prmitives
C     
            DO 40 II = 1, IE
               IEX = IEX + 1 
               JEX = JEXH
               DO 41 JJ = 1, JE
                  JEX = JEX + 1
                  IJIND = IE*(JJ-1) + II
C     
C Get the common center of gaussian on A and on B.
C
                  CALL CIVPT(A, EXPS(IEX), B, EXPS(JEX), EP, GAMAP
     &                        , EFACTP)
                  CALL RHFTCE(AAA, A, EP, ITYP, ITM, ONE, LMN)
                  CALL RHFTCE(BBB, B, EP, JTYP, JTM, ONE, JMN)
C
                  DO 31 K = 1, 3
                     TEMP1(K) = ZERO
                     DO 32 L = 1, 3
                        VAL(K, L)   = ZERO
CSSS                       TEMP2(K, L) = ZERO
 32                  CONTINUE
 31               CONTINUE
C
C Compute the CA_X*[<|X_C/R_C^3|> contribution and store
C in the first three element of variable VAL.
C This term involves multiflication of usual field integral with
C the distance vector AC_X.
C
              
                  DO 25 K = 1, 3
                     DISTQD(K) = EP(K) - G(K)
                     DISTQA(K) = A(K)  - G(K)
 25               CONTINUE
C
                  DO 210 IIQ = 1, ITM
                     IF(AAA(IIQ) .EQ. 0.0) GO TO 210
                     IL = IBTAND(IBTSHR(LMN(IIQ),20), 2**10 - 1)
                     IM = IBTAND(IBTSHR(LMN(IIQ),10), 2**10 - 1)
                     IN = IBTAND(LMN(IIQ), 2**10 - 1)
C
                     DO 211 JJQ = 1, JTM
                        DD = AAA(IIQ)*BBB(JJQ)
                        IF(ABS(DD) .LT. THR2) GO TO 211
                        JL = IBTAND(IBTSHR(JMN(JJQ),20), 2**10 - 1)
     &                       + IL
                        JM = IBTAND(IBTSHR(JMN(JJQ),10), 2**10 - 1)
     &                       + IM
                        JN = IBTAND(JMN(JJQ), 2**10 - 1) + IN
C
                       CALL DIP_MOMENT(JL, JM, JN, GAMAP ,V, NT,
     &                                 DISTQD)
                        DO 51 L = 1, 3
                           TEMP1(L) = EFACTP*DD*V(L) + TEMP1(L)
 51                     CONTINUE
C     
 211                 CONTINUE
 210              CONTINUE
C
CSSS                  Write(6,*) "IJ", I, J 
CSSS                  Write(6,"(3(F15.8))") (A(L), L=1,3)
CSSS                  Write(6,"(3(F15.8))") (B(L), L=1,3)
CSSS                  Write(6,*) 
CSSS                  Write(6,"(3(F12.8))") (TEMP1(L), L=1,3)

                  DO 52 K = 1, 3
                     DO 53 L = 1, 3
                        VAL(K, L) = DISTQA(K)*TEMP1(L) 
 53                  CONTINUE
 52               CONTINUE
CSSS                  Write(6,*) "PRINT 1"
CSSS                  call output(val, 1, 3, 1, 3, 3, 3, 1)
C
C Now compute the other contribution <|X_C*X_C|1/R_C^3|>. Now we have a 
C different gaussian orbital (one having an angular momentum one higher
C than original) and the center remains unchanged.
C 
C Get the angular momentum quantum numbers of the gaussian centered on A
C     
                  LA = LMNVAL(1, ITYP)
                  MA = LMNVAL(2, ITYP)
                  NA = LMNVAL(3, ITYP)
C
                  DO 212 IIQ = 1, 3
                     DO K=1, 3
                        TEMP2(IIQ, K) = ZERO
                     ENDDO
                     CALL GETTYPE (LA + LFLAG(IIQ), MA + MFLAG(IIQ),
     &                             NA + NFLAG(IIQ), IDTYP)
C
C Now get the new IL, IM and IN and initialize LMN array
C
                     CALL RHFTCE (AAA, A, EP, IDTYP, ITM, ONE, LMN)
                     CALL RHFTCE (BBB, B, EP, JTYP, JTM, ONE, JMN)
                    
                     DO 213 JJQ = 1, ITM
                        IF(AAA(JJQ) .EQ. 0.0) GO TO 213
                        IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                        IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                        IN = IBTAND(LMN(JJQ), 2**10 - 1)
C     
                        DO 214 KKQ = 1, JTM
                           DD = AAA(JJQ)*BBB(KKQ)
                           IF(ABS(DD) .LT. THR2) GO TO 214
                           JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1) 
     &                          + IL
                           JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1)
     &                          + IM
                           JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN
C
                       CALL DIP_MOMENT(JL, JM, JN, GAMAP ,V, NT,
     &                                 DISTQD)
C
                           DO 54 L = 1, NT
                              TEMP2(IIQ, L) = EFACTP*DD*V(L) + 
     &                                        TEMP2(IIQ, L)
 54                        CONTINUE
C     
 214                    CONTINUE
 213                 CONTINUE
C
 212              CONTINUE
C
CSSS             do iiq = 1, 3
CSSS             Write(6,"(3(3F12.8))")(TEMP2(IIQ, L), L=1,3)
CSSS             write(6,*)
CSSS             Write(6,"(3(3F12.8))")(VAL(IIQ, L), L=1,3)
CSSS             enddo 
CSSS          Write(6,*)
                  DO 55 MMM = 1, 3
                     DO 56 NNN = 1, 3
                        VAL(MMM, NNN) =   TEMP2(MMM, NNN)
     &                                  + VAL(MMM, NNN)
CSSS                        VAL(MMM, NNN) =   EFACTP*TEMP2(MMM, NNN) 
 56                  CONTINUE
 55               CONTINUE

CSSS                  CALL DSCAL(9, 0.5D0, VAL, 1)

C&&&&&&&DEBUG&&&&
C                  Write(6,*)
C                  Write(6,*) "Uncontracted intgerals"
C                  call output(val, 1, 3, 1, 3, 3, 3, 1)
C&&&&&&&&&&&&&&&&

C&&&&&&&DEBUG&&&&
C                  write(luout,*)
C                  write(luout,*) 'Center1  ', i ,'     ',ITYP
C                  write(luout,*) 'Center2  ', j ,'     ',JTYP
C                  call tab (luout, val, 3, 3, 3, 3)
C&&&&&&&&&&&&&&&&
                  MNX = 0
                  DO 57 N = 1, 3
                     DO 58 K = 1, N
                        MNX = MNX + 1
                        VS(IJIND, MNX) =VAL(N, K) + 
     &                                         VS(IJIND, MNX)
 58                  CONTINUE
 57               CONTINUE
C     
 41            CONTINUE
 40         CONTINUE
C
C Reset the value of NT     
           NT = 6
           DO 101 M = 1, NT
C
C     PERFORM GENERAL CONTRACTION, SPECIAL CODE FOR 1X1 CASE
C     
              IF (IE*JE .EQ. 1) THEN
C                 CONT(ICONT+1) = 1.0D0
C                 CONT(JCONT+1) = 1.0D0
                 VS(1,M) = VS(1,M)*CONT(ICONT+1)*CONT(JCONT+1)
              ELSE
                 CALL MXM (VS(1,M), IE, CONT(JCONT+1), JE, SCR, JCE)
                 CALL MXMT (CONT(ICONT+1), ICE, SCR, IE, VS(1,M), JCE)
              ENDIF
C     
              DO 111 II = 1, ICE
                 IIND = INDX(ICOFF + II)
                 JTOP = JCE
                 IF (I .EQ. J) JTOP = II
                 DO 112 JJ = 1, JTOP
                    JJND = INDX(JCOFF + JJ)
                    VSINT = VS(ICE*(JJ-1)+II, M)
                    IF(ABS(VSINT) .LT. THRESH) GO TO 112
                    IHOLD(M) = IHOLD(M) + 1
                    IPUT = IHOLD(M) + 2
                    IF(JJND.GT.IIND.AND.MTYP.EQ.16) VSINT = -VSINT
                    BUF(IPUT, M) = VSINT
C&&&&&&&DEBUG
C                   Write (luout, 999) BUF(IPUT,M)
C&&&&&&&DEBUG
                    IOFF = IPUT + IWPB
                    IBUF(IOFF*IINTFP,M) = IBTOR(IBTSHL(MAX0(IIND,JJND)
     &                   ,16), MIN0(IIND,JJND))
                    IF(IHOLD(M).LT.IWPB) GO TO 112
                    IBUF(1, M) = IWPB
                    CALL OUTAPD(IBUF(1,M), LENBUF, ID20, LXQ)
                    IBUF(2, M) = LXQ
                    IHOLD(M) = 0
 112             CONTINUE
 111          CONTINUE
C     
 101       CONTINUE
           JEXH = JEXH + JE
           JCONT = JCONT + JE*JCE
           JCOFF = JCOFF + JCE
C     
C&&&&&&DEBUG
C         Write(luout, 999) (BUF(j, m), m = 1, 3)
C&&&&&&DEBUG
 49     CONTINUE
C     
        IEXH = IEXH + IE
        ICONT = ICONT + IE*ICE
        ICOFF = ICOFF + ICE
 50   CONTINUE
C     
C*******DEBUG
C      call tab(luout, buf, lbuf, 3, lbuf, 3)
C&&&&Debug
C     CLOSE BINS
C     
      DO 102 I=1, NT
         IF(IHOLD(I).EQ.0) GO TO 103
         IBUF(1,I) = IHOLD(I)
         CALL OUTAPD(IBUF(1,I), LENBUF, ID20, LXQ)
         IADH(I) = LXQ
         GO TO 102
 103     CONTINUE
C     
         IADH(I) = IBUF(2,I)
 102  CONTINUE
C     
C     NOW TRANSFORM TO SYMMETRY INTS
C     
      NU2 = NSABF*(NSABF+1)/2
      NU = NU2 + 3
      DO 600 NX = 1, NT
C     
         DO 70 I = 1, NU2
            AOINT(I) = ZERO
 70      CONTINUE
C     
         LL=IADH(NX)
 75      CONTINUE
         IF(LL.EQ.-1) GO TO 76
         CALL INTAPD(IBUF, LENBUF, ID20, LL)
         LOOP = IBUF(1, 1)
         LL = IBUF(2, 1)
         DO 74 IQ = 1, LOOP
            X = BUF(IQ + 2, 1)
            IOFF = IQ + 2 + IWPB
            I = IBTAND(IBTSHR(IBUF(IOFF*IINTFP,1), 16), 2**16 - 1)
            J = IBTAND(IBUF(IOFF*IINTFP, 1), 2**16 - 1)
            IE = ITRM(I)
            DO 702 II =1,IE
               CC = CTRN(II,I)
               IS = ITRN(II,I)
               JE = ITRM(J)
               DO 704 JJ = 1, JE
                  JS = ITRN(JJ,J)
                  IF (IS .LT. JS) GO TO 704
                  XP = CC*CTRN(JJ,J)*X
                  IJF = IS*(IS-1)/2+JS
                  AOINT(IJF) = XP + AOINT(IJF)
 704           CONTINUE
 702        CONTINUE
C     
            IF(I .EQ. J) GO TO 74
            DO 722 II = 1, JE
               CC = CTRN(II, J)
               IS = ITRN(II, J)
               DO 724 JJ = 1, IE
                  JS = ITRN(JJ, I)
                  IF(IS .LT. JS) GO TO 724
                  XP = CC*CTRN(JJ, I)*X
                  IF(MTYP .EQ. 16) XP = -XP
                  IJF = IS*(IS-1)/2+JS
                  AOINT(IJF)=XP+AOINT(IJF)
 724           CONTINUE
 722        CONTINUE
 74      CONTINUE
C     
         GO TO 75
 76      CONTINUE
C     
C     CRAPS "INTERFACE"
C     
         NDPROP = 35
         LCBUF = 600
         DUML1 = '********'
C     
C SG 12/6/96
         DUMLAB = 0.0D0
C
         WRITE(NDPROP) DUML1, DUMLAB, DUMLAB, TITLE(NX, MTYP),
     &      PTOT(1,NX)
         ICOUNT = 0
C     
#ifdef _DEBUG_LVLM
         If  (mtyp .eq. 21) Then
         Write(6,*) "Triangularly Packed AO SM integrals"
         Write(6,"(5F12.8)") (aoint(i), i=1, nu2)
         sum = 0.0d0
         do i = 1, nu2
            sum = sum + aoint(i)*aoint(i)
         enddo
CSSS         write(6, "(a,F15.10)") "The check sum = ", sum
         Endif
#endif

         IF (MTYP .EQ. 4) THEN
            CALL PUTREC(20, 'JOBARC', LABELD(NX), NU2*IINTFP, AOINT)
         ELSEIF (MTYP .EQ. 3 .OR. MTYP .EQ. 9) THEN
            IF (ICT .LE. 9) THEN 
               JUNKLAB = LABELJA(NX, MTYP)
               WRITE(JUNKLAB(8:8),'(I1)') ICT
            ELSEIF (ICT .GT.9 .AND. ICT .LT. 100) THEN
               WRITE(JUNKLAB(7:8),'(I2)') ICT
            ELSE
               WRITE(JUNKLAB(6:8),'(I3)') ICT
            ENDIF
C     
            CALL PUTREC(20, 'JOBARC', JUNKLAB, NU2*IINTFP, AOINT)
         ELSE
            CALL PUTREC(20, 'JOBARC', LABELJA(NX,MTYP), NU2*IINTFP,
     &         AOINT)
         ENDIF
C     
         DO 750 IJT = 1, NU2
            If (ICOUNT .LT. LCBUF) GOTO 751
            LIMIT = LCBUF
            WRITE(NDPROP) GBUF, IGNDX, LIMIT
            ICOUNT = 0
 751        ICOUNT = ICOUNT + 1
            IGNDX(ICOUNT) = IJT
            GBUF(ICOUNT) = AOINT(IJT)
 750     CONTINUE
C     
         IF (ICOUNT .gt. 0) THEN
            LIMIT = ICOUNT
            WRITE(NDPROP) GBUF, IGNDX,LIMIT
         ENDIF
C     
         LIMIT =  -1
         WRITE(NDPROP) GBUF, IGNDX, LIMIT
C     
C     END OF "INTERFACE" 
C     
         IF(ISL2.EQ.0) GO TO 600
C     CALL PRINTM(AOINT,NCONT)
C     
 600  CONTINUE
C     
      TI = T2 - T1
      TS = T3 - T2
      TR = T3 - T1
C     
C     WRITE(6,530) TI, TS, TR
 530  FORMAT(' INTS ',F8.3,' SYMMETRY ',F8.3,' TOTAL ',F8.3)
C     CLOSE(26)
C     
      RETURN
      END
