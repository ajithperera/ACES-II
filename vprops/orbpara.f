C
      SUBROUTINE OrbPara(VLIST,NOC,CENTR1,EXPS,CONT,LABEL,NPRIM,NCONTR,
     &                  ITRM,ITRN,CTRN,NCONT,BUF,IBUF,AOINT,LBUF,INDX,
     &                  NGROUP,VS,SCR,MXPRM2)
C
C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C &  Evaluates Orbital Paramagnetic contribution to the NMR coupling &
C &  constant. Coded by Ajith 08/93 following vprop.f style          &
C &  (Not the most efficient or most readable program).              &
C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND
      CHARACTER*8 DUML1,LABELD(3),LABELJA(10,25),JUNKLAB
      CHARACTER*8 TITLE
C
      COMMON /PRTBUG/ IPRTFS
      COMMON /PRTDAT/ IDOSPH, MAXAQN, NSABF
      COMMON /NUCINDX/ ICT
      COMMON /PARM/ISL2,ISL3,ISL4,ISLX,THRESH,MTYP,ILNMCM,INT,INUI,
     &             NCENTS,ID20,NSYM,NAO(12),NMO(12),JUNK,PTOT(3,10), 
     &             INTERF,JOBIPH,RUNMOS
      COMMON /HIGHL/ LMNVAL(3,84),ANORM(84)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /VARS/ FNU(20),DP(20),FN(15),FD(15),PI,DFTR(11)
      COMMON /FILES/LUOUT,MOINTS
C
      DIMENSION VLIST(NOC,10),EXPS(1),CONT(1),LABEL(1),ITRM(1),
     &          ITRN(256,1),CTRN(256,1),NCONTR(1),NPRIM(1),INDX(1)
      DIMENSION IHOLD(10),IADH(10),VS(MXPRM2,10),TEMP(2,3)
      DIMENSION BUF(LBUF,10),IBUF(IINTFP*LBUF,10),AOINT(1)
      DIMENSION TITLE(10,25)
      DIMENSION CENTR1(3),SCR(1),GBUF(600),IGNDX(600)
      DIMENSION LMN(200),JMN(200)
      DIMENSION A(3),B(3),EP(3),DISTQD(3),AAA(200),BBB(200),
     &          V(10),OP(3)
C     
      DATA (LABELJA(I, 18), I = 1, 10) /'   OPX  ', '   OPY  ',
     &      '   OPZ  ',7*'        '/
C     
      DATA (TITLE(I, 18), I = 1, 10)/'   OPX  ','   OPY  ','   OPZ  ',
     &                             7*'        '/
C     DATA (TITLE(I, 18), I = 1, 10)/6H   OPX, 6H   OPY, 6H   OPZ
C    &      , 7*6H      /
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
      NT = 3
      DO 100 I = 1, NT
         IHOLD(I) = 0
         IBUF(1,I) = 0
         IBUF(2,I) = -1
100   CONTINUE
C
      IEXH  = 0
      ICONT = 0
      ICOFF = 0
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
 410    CONTINUE

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
 411       CONTINUE
C 
           JE = NPRIM(J)
           JCE = NCONTR(J)
C
           DO 2 M = 1, 3
              DO 28 MM = 1, MXPRM2
                 VS(MM, M) = ZERO
 28           CONTINUE
 2         CONTINUE
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
C Get the common center of gaussian on A and on B. The resulting
C center remains the same for new gaussians resulting from the
C differentiation
                 CALL CIVPT(A, EXPS(IEX), B, EXPS(JEX), EP, GAMAP
     &                      , EFACTP)

C     
C Get the angular momentum quantum numbers of the gaussian
C centered on B
                 LB = LMNVAL(1, JTYP)
                 MB = LMNVAL(2, JTYP)
                 NB = LMNVAL(3, JTYP)
                 
C
C Lets form the three intermediates resulting from the derivatives
C of the Gaussian on B. First with respect to X
C
                 CALL RHFTCE (AAA, A, EP, ITYP, ITM, ONE, LMN)
                 DO 25 K = 1, 3
                    DISTQD(K) = EP(K) - CENTR1(K)
 25              CONTINUE
                 
C
C Get the type of the function results from diferentiation
C with respect to X. The cases where LB, MB and/or NB  are 
C zero handeled as they were non-zero to avoid testing for
C those inside the do loops (see subroutine gettype). 
C Contribution from those terms will vanish when we collect
C those for the total derivative with appropriate factos
C (LB, MB and NB). 
C
                 CALL GETTYPE (LB - 1, MB, NB, JDTYP)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
                 DO 1155 IIA = 1, 3
                    DO 1156 JJA = 1, 2
                       TEMP(JJA, IIA) = ZERO
 1156               CONTINUE
 1155            CONTINUE
C     
                 DO 210 IIQ = 1, ITM  
                    IF (AAA(IIQ) .EQ. 0.0) GO TO 210
                    IL = IBTAND(IBTSHR(LMN(IIQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(IIQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(IIQ), 2**10 - 1)
C     
                    DO 211 JJQ = 1, JTM
                       DD = AAA(IIQ)*BBB(JJQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 211
                       JL = IBTAND(IBTSHR(JMN(JJQ),20), 2**10 - 1) + IL
                       JM = IBTAND(IBTSHR(JMN(JJQ),10), 2**10 - 1) + IM
                       JN = IBTAND(JMN(JJQ), 2**10 - 1) + IN
C
C Now everything is ready to do the actual evaluation of field integral
C     
                       CALL EFINT(JL, JM, JN, GAMAP ,V, NT, DISTQD)
C       
                       DO 51 L = 1, 3
                          TEMP(1, L) = EFACTP*DD*V(L) + TEMP(1,L)
 51                    CONTINUE
 211                CONTINUE
 210             CONTINUE
C
                 CALL GETTYPE (LB + 1, MB, NB, JDTYP)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C     
                 DO 212 IIQ = 1, ITM 
                    IF (AAA(IIQ) .EQ. 0.0) GO TO 212
                    IL = IBTAND(IBTSHR(LMN(IIQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(IIQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(IIQ), 2**10 - 1)
C     
                    DO 213 JJQ = 1, JTM
                       DD = AAA(IIQ)*BBB(JJQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 213
                       JL = IBTAND(IBTSHR(JMN(JJQ),20), 2**10 - 1) + IL
                       JM = IBTAND(IBTSHR(JMN(JJQ),10), 2**10 - 1) + IM
                       JN = IBTAND(JMN(JJQ), 2**10 - 1) + IN
C
                       CALL EFINT(JL, JM, JN, GAMAP, V, NT, DISTQD)
C
                       DO 52 L = 1, 3
                          TEMP(2, L) = EFACTP*DD*V(L) + TEMP(2, L)
 52                    CONTINUE
 213                CONTINUE
 212             CONTINUE
C
C Now evaluate the contribution to X-derivative
C
                 DZX = LB*TEMP(1, 3) - TWO*EXPS(JEX)*TEMP(2, 3)
                 DYX = LB*TEMP(1, 2) - TWO*EXPS(JEX)*TEMP(2, 2)
C
C Evalute the derivative intermediate with respect to Y
C
                 CALL GETTYPE (LB , MB - 1, NB, JDTYP)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
                 DO 1159 IIA = 1, 3
                    DO 1160 JJA = 1, 2
                       TEMP(JJA, IIA) = ZERO
 1160               CONTINUE
 1159            CONTINUE
C     
                 DO 214 IIQ = 1, ITM 
                    IF (AAA(IIQ) .EQ. 0.0) GO TO 214
                    IL = IBTAND(IBTSHR(LMN(IIQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(IIQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(IIQ), 2**10 - 1)
C     
                    DO 215 JJQ = 1, JTM
                       DD = AAA(IIQ)*BBB(JJQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 215
                       JL = IBTAND(IBTSHR(JMN(JJQ),20), 2**10 - 1) + IL
                       JM = IBTAND(IBTSHR(JMN(JJQ),10), 2**10 - 1) + IM
                       JN = IBTAND(JMN(JJQ), 2**10 - 1) + IN
C
                       CALL EFINT(JL, JM, JN, GAMAP, V, NT, DISTQD)
C
                       DO 53 L = 1, 3
                          TEMP(1, L) = EFACTP*DD*V(L) + TEMP(1, L)
 53                    CONTINUE
 215                CONTINUE
 214             CONTINUE
C
                 CALL GETTYPE (LB, MB + 1, NB, JDTYP)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C
                 DO 216 IIQ = 1, ITM 
                    IF (AAA(IIQ) .EQ. 0.0) GO TO 216
                    IL = IBTAND(IBTSHR(LMN(IIQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(IIQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(IIQ), 2**10 - 1)
C     
                    DO 217 JJQ = 1, JTM
                       DD = AAA(IIQ)*BBB(JJQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 217
                       JL = IBTAND(IBTSHR(JMN(JJQ),20), 2**10 - 1) + IL
                       JM = IBTAND(IBTSHR(JMN(JJQ),10), 2**10 - 1) + IM
                       JN = IBTAND(JMN(JJQ), 2**10 - 1) + IN
C
                       CALL EFINT(JL, JM, JN, GAMAP, V, NT, DISTQD)
C
                       DO 54 L = 1, 3
                          TEMP(2, L) = EFACTP*DD*V(L) + TEMP(2, L)
 54                    CONTINUE
 217                CONTINUE
 216             CONTINUE
C
C Now evaluate the contribution to Y-derivative
C
                 DZY = MB*TEMP(1, 3) - TWO*EXPS(JEX)*TEMP(2, 3)
                 DXY = MB*TEMP(1, 1) - TWO*EXPS(JEX)*TEMP(2, 1)
C
C Evalute the derivative intermediate with respect to Z
C
                 CALL GETTYPE (LB , MB, NB - 1, JDTYP)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C
                 DO 1165 IIA = 1, 3
                    DO 1166 JJA = 1, 2
                       TEMP(JJA, IIA) = ZERO
 1166                  CONTINUE
 1165               CONTINUE
C     
                 DO 218 IIQ = 1, ITM 
                    IF (AAA(IIQ) .EQ. 0.0) GO TO 218
                    IL = IBTAND(IBTSHR(LMN(IIQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(IIQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(IIQ), 2**10 - 1)
C     
                    DO 219 JJQ = 1, JTM
                       DD = AAA(IIQ)*BBB(JJQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 219
                       JL = IBTAND(IBTSHR(JMN(JJQ),20), 2**10 - 1) + IL
                       JM = IBTAND(IBTSHR(JMN(JJQ),10), 2**10 - 1) + IM
                       JN = IBTAND(JMN(JJQ), 2**10 - 1) + IN
C
                       CALL EFINT(JL, JM, JN, GAMAP, V, NT, DISTQD)
C
                       DO 55 L = 1, 3
                          TEMP(1, L) = EFACTP*DD*V(L) + TEMP(1, L)
 55                    CONTINUE
 219                CONTINUE
 218             CONTINUE
C
                 CALL GETTYPE (LB, MB, NB + 1, JDTYP)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C     
                 DO 220 IIQ = 1, ITM 
                    IF (AAA(IIQ) .EQ. 0.0) GO TO 220
                    IL = IBTAND(IBTSHR(LMN(IIQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(IIQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(IIQ), 2**10 - 1)
C     
                    DO 221 JJQ = 1, JTM
                       DD = AAA(IIQ)*BBB(JJQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 221
                       JL = IBTAND(IBTSHR(JMN(JJQ),20), 2**10 - 1) + IL
                       JM = IBTAND(IBTSHR(JMN(JJQ),10), 2**10 - 1) + IM
                       JN = IBTAND(JMN(JJQ), 2**10 - 1) + IN
C
                       CALL EFINT(JL, JM, JN, GAMAP, V, NT, DISTQD)
C
                       DO 56 L = 1, 3
                          TEMP(2, L) = EFACTP*DD*V(L) + TEMP(2,L)
 56                    CONTINUE
 221                CONTINUE
 220             CONTINUE
C
C Now evaluate the contribution to the Z-derivative
C
                 DXZ = NB*TEMP(1, 1) - TWO*EXPS(JEX)*TEMP(2, 1)
                 DYZ = NB*TEMP(1, 2) - TWO*EXPS(JEX)*TEMP(2, 2)
C
C Finally sum up appropriate terms to get the actual integral
C 
                 OP(1) = DYZ - DZY
                 OP(2) = DXZ - DZX
                 OP(3) = DXY - DYX
C
C&&&&&&&DEBUG&&&&
C                 write(luout,*)
C                 write(luout,*) 'Center1  ', i ,'     ',ITYP
C                 write(luout,*) 'Center2  ', j ,'     ',JTYP
C                 Write (luout, 999) (OP(IIII), IIII = 1 ,3)
C 999             format(5X,F12.5, 5X,F12.5,5X,F12.5) 
C&&&&&&&&&&&&&&&&
                 DO 30 MM = 1, NT
                    VS(IJIND, MM) = OP(MM) + VS(IJIND, MM)
 30              CONTINUE
C&&&&&&&DEBUG&&&&
C                 write(luout,*)
C                 write(luout,*) 'Center1  ', i ,'     ',ITYP
C                 write(luout,*) 'Center2  ', j ,'     ',JTYP
C                 Write (luout, 999) (VS(IJIND, IIII), IIII = 1 ,3)
C&&&&&&&&&&&&&&&&
C
 41           CONTINUE
 40        CONTINUE
C     
           DO 101 M = 1, NT
C
C PERFORM GENERAL CONTRACTION, SPECIAL CODE FOR 1X1 CASE
C
              IF (IE*JE .EQ. 1) THEN
                 VS(1,M) = VS(1,M)*CONT(ICONT+1)*CONT(JCONT+1)
              ELSE
                 CALL MXM (VS(1,M), IE, CONT(JCONT+1), JE, SCR, JCE)
                 CALL MXMT (CONT(ICONT+1), ICE, SCR, IE, VS(1,M), JCE)
              ENDIF
C&&&&&&&&&DEBUG
C              Write (luout, 999) (VS(1,IIII), IIII = 1 ,3)
C&&&&&&&&&&&&          
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
                    IOFF = IPUT + IWPB
                    IBUF(IOFF*IINTFP, M) = IBTOR(IBTSHL(MAX0(IIND,JJND)
     &                                    ,16), MIN0(IIND,JJND))
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
 49     CONTINUE
        IEXH = IEXH + IE
        ICONT = ICONT + IE*ICE
        ICOFF = ICOFF + ICE
 50   CONTINUE
C     
C     CLOSE BINS
C     
C&&&&&&debug
C      nmad=3
C      do 17 I = 1, 3
C         WRITE (LUOUT, *) 'PSO-INTEGRALS'
C         call tab(luout, buf(NMAD, I), NSABF, NSABF, NSABF, NSABF)
C         NMAD = NMAD+NSABF*NSABF
C 17   continue
C&&&&Debug
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
                  IF(MTYP .EQ. 18) XP = -XP
                  IJF = IS*(IS-1)/2+JS
                  AOINT(IJF)= XP + AOINT(IJF)
 724           CONTINUE
 722        CONTINUE
 74      CONTINUE
C     
         GO TO 75
 76      CONTINUE
C     
C CRAPS "INTERFACE"
C     
         NDPROP = 35
         LCBUF = 600
         DUML1 = '********'
C
C SG 12/6/96
         DUMLAB = 0.0D0
C
         WRITE(NDPROP) DUML1, DUMLAB, DUMLAB, TITLE(NX, MTYP),
     &                 PTOT(1,NX)
         ICOUNT = 0
C     
C&&&&&&&&debug
C         call tab(luout, aoint, NSABF, NSABF, NSABF, NSABF)
C&&&&&&&&&debug
         IF (MTYP .EQ. 4) THEN
            CALL  PUTREC(20, 'JOBARC', LABELD(NX), NU2*IINTFP, AOINT)
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
     &                  AOINT)
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
C  END OF "INTERFACE" 
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
      RETURN
      END
