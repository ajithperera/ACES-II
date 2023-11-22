C
      SUBROUTINE OrbDia(VLIST,NOC,CENTR1,CENTR2,EXPS,CONT,
     &   LABEL,NPRIM,NCONTR,ITRM,ITRN,CTRN,NCONT,BUF,
     &   IBUF,AOINT,LBUF,INDX,NGROUP,VS,SCR,MXPRM2,
     &   ROOT,WEIGHT,VALUEM,NGAUSS)
C     
C     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     &  Evaluates Orbital Diamagnetic contribution to the NMR coupling  &
C     &  constant. Coded by Ajith 08/93 following vprop.f style          &
C     &  (Not the most efficient or most readable program).              &
C     &  Implementation is based on JCP 73(11), 5719, 1980               &
C     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND, NGAUSS
      CHARACTER*8 DUML1,LABELD(3),LABELJA(10,25),JUNKLAB
      CHARACTER*8 TITLE
C     
      COMMON /PRTBUG/ IPRTFS
      COMMON /PRTDAT/ IDOSPH, MAXAQN, NSABF
      COMMON /NUCINDX/ ICT
      COMMON /PARM/ISL2,ISL3,ISL4,ISLX,THRESH,MTYP,ILNMCM,INT,INUI,
     &   NCENTS,ID20,NSYM,NAO(12),NMO(12),JUNK,PTOT(3,10), 
     &   INTERF,JOBIPH,RUNMOS
      COMMON /HIGHL/ LMNVAL(3,84), ANORM(84)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /VARS/ FNU(20),DP(20),FN(15),FD(15),PI,DFTR(11)
      COMMON /FILES/LUOUT,MOINTS
C     
      DIMENSION VLIST(NOC,10),EXPS(1),CONT(1),LABEL(1),ITRM(1),
     &   ITRN(256,1),CTRN(256,1),NCONTR(1),NPRIM(1),INDX(1)
      DIMENSION IHOLD(10),IADH(10), VS(MXPRM2,10)
      DIMENSION BUF(LBUF,10),IBUF(IINTFP*LBUF,10),AOINT(1)
      DIMENSION VAL(3,3),TITLE(10,25),TEMP(3)
      DIMENSION CENTR1(3),CENTR2(3),SCR(1),GBUF(600),IGNDX(600)
      DIMENSION LMN(200),JMN(200),KMN(200)
      DIMENSION A(3),B(3),EP(3),EQ(3),DISTQD(3),AAA(200),BBB(200),
     &   CCC(200),V(10)
      DIMENSION ROOT(24),WEIGHT(24)
C     
      DATA (LABELJA(I,17), I = 1, 10) /'  ODXX  ', '  ODYY  ',
     &   '  ODZZ  ','  ODXY  ', '  ODXZ  ', '  ODYZ  ' ,
     &   4*'        '/
C     
      DATA (TITLE(I,17), I = 1, 10)/'  ODXX  ','  ODYY  ','  ODZZ  ',
     &                              '  ODXY  ','  ODXZ  ','  ODYZ  ',
     &                            4*'        '/
C     DATA (TITLE(I,17), I = 1, 10)/6H  ODXX, 6H  ODYY, 6H  ODZZ
C    &   , 6H  ODXY, 6H  ODXZ, 6H  ODYZ, 4*6H      /
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
      IWPB = (LBUF - 2)/2
      LENBUF = IINTFP*LBUF
      ZERO = 0.0D+00
      ONE  = 1.00D+00
      EIGHT = 8.00D+00
C     
C Do sum initializations
      NT = 3
      CALL SETRHF
C      
      DO 100 I = 1, 3
         IHOLD(I) = 0
         IBUF(1,I) = 0
         IBUF(2,I) = -1
 100  CONTINUE
C     
      IEXH  = 0
      ICONT = 0
      ICOFF = 0
C     
C     Loop over the distinct group of exponents
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
            DO 2 M = 1, 3
               DO 28 MM = 1, MXPRM2
                  VS(MM, M) = ZERO
 28            CONTINUE
 2          CONTINUE
C     
            IEX = IEXH
            
C     Loop over the prmitives
C     
            DO 40 II = 1, IE
               IEX = IEX + 1 
               JEX = JEXH
               DO 41 JJ = 1, JE
                  JEX = JEX + 1
                  IJIND = IE*(JJ-1) + II
C     
C     Note the redundancy here for contracted basis functions
C     since the evaluation of orbital diamagnetic integrals
C     are slower this redundancy will be eliminated with 
C     appropriate scheme
C
                  DO 20 K = 1, 3
                     DO 21 L = 1, 3
                        VAL(K, L) =  ZERO
 21                  CONTINUE
 20               CONTINUE
C     
C     Loop over number of integration points (Gauss-Legendre numerical
C     integration. 
                  DO 1999 IROOT = 1, NGAUSS
C     
C     Evaluate the Alpha using m as 0.8 (defined as a parameter)
C     JCP 73(11), 5719, 1980
C     
                     TVALUE = ROOT(IROOT)
                     WVALUE = WEIGHT(IROOT)
C     
                     ALPHA = VALUEM*(ONE + TVALUE)/(ONE - TVALUE)
                     ALPHAC = ALPHA*ALPHA
C     
C     Factor contributing from this root Fac
C     
                     FAC1 = ((ONE + TVALUE)**2)/((ONE - TVALUE)**4)
                     FAC2 = (VALUEM**3)*WVALUE*EIGHT/DSQRT(PI)
C     
C     Calculate the new center and exponent of product of Gaussian on 
C     Center A and Gaussian on center C (Really this is from the transformation
C     of X_C/R(3)_C and always a P function). Resulting Gaussian is centered 
C     on P. Then take the product of Gaussian on A and on P. Resulting Gaussian
C     is centered on Q.
C     
                     CALL CIVPT(CENTR1, ALPHAC, B, EXPS(JEX), EP, GAMAP
     &                  , EFACTP)
                     CALL CIVPT(A, EXPS(IEX), EP, GAMAP, EQ, GAMAQ,
     &                  EFACTQ)
                     
C     Take the product of factors
C     
                     EFAC = EFACTP*EFACTQ
                     CALL RHFTCE (AAA, A, EQ, ITYP, ITM, ONE, LMN)
                     CALL RHFTCE (BBB, B, EQ, JTYP, JTM, ONE, JMN)
C     
                     DO 42 KK = 1, 3
                        KTYP = KK + 1
                        CALL RHFTCE (CCC,CENTR1,EQ,KTYP,KTM,ONE,KMN )
C     
                        DO 25 K = 1, 3
                           DISTQD(K) = EQ(K)- CENTR2(K)
                           TEMP(K) = ZERO
 25                     CONTINUE
C     
                        DO 210 IIQ = 1, ITM
                           IF (AAA(IIQ) .EQ. 0.0) GO TO 210
                           IL = IBTAND(IBTSHR(LMN(IIQ),20), 2**10 - 1)
                           IM = IBTAND(IBTSHR(LMN(IIQ),10), 2**10 - 1)
                           IN = IBTAND(LMN(IIQ), 2**10 - 1)
C     
                           DO 211 JJQ = 1, JTM
                              JL = IBTAND(IBTSHR(JMN(JJQ),20), 2**10 - 1
     &                           ) + IL
                              JM = IBTAND(IBTSHR(JMN(JJQ),10), 2**10 - 1
     &                           ) + IM
                              JN = IBTAND(JMN(JJQ), 2**10 - 1) + IN
C     
                              DO 212 KKQ = 1, KTM
                                 DD = AAA(IIQ)*BBB(JJQ)*CCC(KKQ)
                                 IF(ABS(DD) .LT. THR2) GO TO 212
                                 KL =IBTAND(IBTSHR(KMN(KKQ),20), 2**10 
     &                              - 1) + JL
                                 KM =IBTAND(IBTSHR(KMN(KKQ),10), 2**10 
     &                              - 1)  + JM
                                 KN = IBTAND(KMN(KKQ), 2**10 - 1) + JN
C     
C     Now everything is ready to do the actual evaluation of field integral.
C     
                                 CALL EFINT(KL,KM,KN,GAMAQ,V,NT,DISTQD)
C     
                                 DO 51 L = 1, NT
                                    TEMP(L) = DD*V(L) + TEMP(L)
 51                              CONTINUE
C     
 212                          CONTINUE
 211                       CONTINUE
 210                    CONTINUE
C     
                        DO 52 L = 1, 3
                           VAL(KK, L) = FAC2*FAC1*EFAC*TEMP(L) +
     &                                  VAL(KK, L)
 52                     CONTINUE
C     
 42                  CONTINUE
C 
 1999             CONTINUE                    
C     
C&&&&&&&DEBUG&&&&
C     write(luout,*)
C     write(luout,*) 'Center1  ', i ,'     ',ITYP
C     write(luout,*) 'Center2  ', j ,'     ',JTYP
C     call tab (luout, val, 3, 3, 3, 3)
C&&&&&&&&&&&&&&&&
                  DO 30 MM = 1, 3 
                     VS(IJIND, MM) = VAL(MM, MM) + VS(IJIND, MM)
 30               CONTINUE
C&&&&&&&DEBUG&&&&
C                  write(luout,*)
C                  write(luout,*) 'Center1  ', i ,'     ',ITYP
C                  write(luout,*) 'Center2  ', j ,'     ',JTYP
C                  Write (luout, 999) (VS(IJIND, IIII), IIII = 1 ,3)
C 999              format(5X,F12.5, 5X,F12.5,5X,F12.5) 
C&&&&&&&&&&&&&&&&
C     
 41            CONTINUE
 40         CONTINUE
C     
            DO 101 M = 1, NT
C               
C     PERFORM GENERAL CONTRACTION, SPECIAL CODE FOR 1X1 CASE
C     
               IF (IE*JE .EQ. 1) THEN
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
C                     Write (luout, 999) BUF(IPUT,M)
C&&&&&&&DEBUG
                     IOFF = IPUT + IWPB
                     IBUF(IOFF*IINTFP,M) = IBTOR(IBTSHL(MAX0(IIND,JJND)
     &                  ,16), MIN0(IIND,JJND))
                     IF(IHOLD(M).LT.IWPB) GO TO 112
                     IBUF(1, M) = IWPB
                     CALL OUTAPD(IBUF(1,M), LENBUF, ID20, LXQ)
                     IBUF(2, M) = LXQ
                     IHOLD(M) = 0
 112              CONTINUE
 111           CONTINUE
C     
 101        CONTINUE
            JEXH = JEXH + JE
            JCONT = JCONT + JE*JCE
            JCOFF = JCOFF + JCE
C     
C&&&&&&DEBUG
C            Write(luout, 999) (BUF(j, m), m = 1, 3)
C&&&&&&DEBUG
 49      CONTINUE
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
C     CLOSE(26)
C     
      RETURN
      END
