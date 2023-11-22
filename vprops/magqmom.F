C
      SUBROUTINE Magqmom(VLIST,NOC,CENTR1,EXPS,CONT,LABEL,NPRIM,NCONTR,
     &                  ITRM,ITRN,CTRN,NCONT,BUF,IBUF,AOINT,LBUF,INDX,
     &                  NGROUP,VS,SCR,MXPRM2)
C
C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C &  Evaluates magnetic quadrupole momentum AO integrals             &
C &  Ajith Perera, 22/2020                                           &
C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND
      CHARACTER*8 DUML1,LABELD(9),LABELJA(9,26),JUNKLAB
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
      DIMENSION IHOLD(10),IADH(10),VS(MXPRM2,10),TEMP(4,3)
      DIMENSION BUF(LBUF,10),IBUF(IINTFP*LBUF,10),AOINT(1)
      DIMENSION TITLE(9,26)
      DIMENSION CENTR1(3),SCR(1),GBUF(600),IGNDX(600)
      DIMENSION LMN(200),JMN(200)
      DIMENSION A(3),B(3),EP(3),DISTQD(3),AAA(200),BBB(200),
     &          V(10),DV(9)
C     
      DATA LABELD /'   MQXX ', '   MQYX ','   MQZX ','   MQXY ',
     &             '   MQYY ', '   MQZY ','   MQXZ ','   MQYZ ',
     &             '   MQZZ '/
      DATA (TITLE(I, 26), I = 1, 9)/'   MQXX ','   MQYX ','   MQZX ',
     &                              '   MQXY ','   MQYY ','   MQZY ',
     &                              '   MQXZ ','   MQYZ ','   MQZZ '/
      DATA (LABELJA(I, 26), I = 1, 9)
     &                              /'   MQXX ','   MQYX ','   MQZX ',
     &                               '   MQXY ','   MQYY ','   MQZY ',
     &                               '   MQXZ ','   MQYZ ','   MQZZ '/

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
C Do some initializations (No nuclear contribution for this term)
C
      CALL SETRHF
      NT = 9
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
           DO 2 M = 1, 10
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
C
                 CALL CIVPT(A, EXPS(IEX), B, EXPS(JEX), EP, GAMAP,
     &                      EFACTP)
C     
C Get the angular momentum quantum numbers of the gaussian
C centered on B
C
                 LB = LMNVAL(1, JTYP)
                 MB = LMNVAL(2, JTYP)
                 NB = LMNVAL(3, JTYP)
C
C Lets form the three intermediates resulting from the derivatives
C of the Gaussian on B. First with respect to X
C
                 DO 25 K = 1, 3
                    DISTQD(K) = EP(K) - CENTR1(K)
 25              CONTINUE
C
C Built the (x,y,z)l_x = (x,y,z)(yPz - zpy), First (x,y,z)yPz
C
                 CALL GETTYPE (LB, MB, NB - 1, JDTYP)
                 CALL RHFTCE (AAA, A, EP, ITYP, ITM, ONE, LMN)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C
                 CALL DZERO(TEMP,12)
C
                 DO JJQ = 1, ITM
                    IF (AAA(JJQ) .EQ. 0.0) GO TO 211
                    IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(JJQ), 2**10 - 1)
C
                    DO KKQ = 1, JTM
                       DD = AAA(JJQ)*BBB(KKQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 212
                       JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1)
     &                      + IL
                       JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1)
     &                      + IM
                       JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN
C
C Now everything is ready to do the actual evaluation of overlap integral
C
                       CALL SEC_MOMENT(JL, JM, JN, GAMAP ,V, NT1, 
     &                                 DISTQD)
C
                       TEMP(1, 1) = EFACTP*DD*V(4) + TEMP(1, 1)
                       TEMP(1, 2) = EFACTP*DD*V(2) + TEMP(1, 2)
                       TEMP(1, 3) = EFACTP*DD*V(6) + TEMP(1, 3)
C
  212                  CONTINUE
                    ENDDO
  211               CONTINUE
                 ENDDO
C
                 CALL GETTYPE (LB, MB, NB + 1, JDTYP)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C
                 DO JJQ = 1, ITM
                    IF (AAA(JJQ) .EQ. 0.0) GO TO 213
                    IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(JJQ), 2**10 - 1)
C
                    DO KKQ = 1, JTM
                       DD = AAA(JJQ)*BBB(KKQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 214
                       JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1)
     &                      + IL
                       JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1)
     &                      + IM
                       JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN
C
                       CALL SEC_MOMENT(JL, JM, JN, GAMAP, V, NT1, 
     &                                 DISTQD)
C
                       TEMP(2, 1) = EFACTP*DD*V(4) + TEMP(2, 1)
                       TEMP(2, 2) = EFACTP*DD*V(2) + TEMP(2, 2)
                       TEMP(2, 3) = EFACTP*DD*V(6) + TEMP(2, 3)
C
  214                  CONTINUE
                     ENDDO
  213                CONTINUE
                 ENDDO
C Now (x,y,z)zpy
                 CALL GETTYPE (LB, MB - 1, NB, JDTYP)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C
                 DO JJQ = 1, ITM
                    IF (AAA(JJQ) .EQ. 0.0) GO TO 215
                    IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(JJQ), 2**10 - 1)
C
                    DO KKQ = 1, JTM
                       DD = AAA(JJQ)*BBB(KKQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 216
                       JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1)
     &                      + IL
                       JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1)
     &                      + IM
                       JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN
C
C Now everything is ready to do the actual evaluation of overlap integral
C
                       CALL SEC_MOMENT(JL, JM, JN, GAMAP, V, NT1, 
     &                                 DISTQD)
C
                       TEMP(3, 1) = EFACTP*DD*V(5) + TEMP(3, 1)
                       TEMP(3, 2) = EFACTP*DD*V(6) + TEMP(3, 2)
                       TEMP(3, 3) = EFACTP*DD*V(3) + TEMP(3, 3)
C
  216                  CONTINUE
                    ENDDO
  215               CONTINUE
                 ENDDO
C
                 CALL GETTYPE (LB, MB + 1, NB, JDTYP)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C
                 DO JJQ = 1, ITM
                    IF (AAA(JJQ) .EQ. 0.0) GO TO 217
                    IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(JJQ), 2**10 - 1)
C
                    DO KKQ = 1, JTM
                       DD = AAA(JJQ)*BBB(KKQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 218
                       JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1)
     &                      + IL
                       JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1)
     &                      + IM
                       JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN
C
                       CALL SEC_MOMENT(JL, JM, JN, GAMAP, V, NT1,
     &                                 DISTQD)
C
                       TEMP(4, 1) = EFACTP*DD*V(5) + TEMP(4, 1)
                       TEMP(4, 2) = EFACTP*DD*V(6) + TEMP(4, 2)
                       TEMP(4, 3) = EFACTP*DD*V(3) + TEMP(4, 3)
C
  218                  CONTINUE
                     ENDDO
  217                CONTINUE
                 ENDDO

C
C Now MQ(xx,xy,xz)  = -(xypz - xzpy,yypz - yzpy,zypz - zzpy)
C
                 DV(1) = (NB*TEMP(1, 1) - TWO*EXPS(JEX)*TEMP(2, 1)) -
     &                   (MB*TEMP(3, 1) - TWO*EXPS(JEX)*TEMP(4, 1))
                 DV(2) = (NB*TEMP(1, 2) - TWO*EXPS(JEX)*TEMP(2, 2)) -
     &                   (MB*TEMP(3, 2) - TWO*EXPS(JEX)*TEMP(4, 2))
                 DV(3) = (NB*TEMP(1, 3) - TWO*EXPS(JEX)*TEMP(2, 3)) -
     &                   (MB*TEMP(3, 3) - TWO*EXPS(JEX)*TEMP(4, 3))

C Built the (x,y,z)l_y = -(x,y,z)(xPz - zpx), First (x,y,z)xpz
C
                 CALL GETTYPE (LB, MB, NB - 1, JDTYP)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C
                 CALL DZERO(TEMP,12)
C
                 DO JJQ = 1, ITM
                    IF (AAA(JJQ) .EQ. 0.0) GO TO 219
                    IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(JJQ), 2**10 - 1) 
C
                    DO KKQ = 1, JTM
                       DD = AAA(JJQ)*BBB(KKQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 220
                       JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1)
     &                      + IL
                       JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1)
     &                      + IM
                       JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN
C
C Now everything is ready to do the actual evaluation of overlap integral
C                
                       CALL SEC_MOMENT(JL, JM, JN, GAMAP, V, NT1,
     &                                 DISTQD)
C                
                       TEMP(1, 1) = EFACTP*DD*V(1) + TEMP(1, 1)
                       TEMP(1, 2) = EFACTP*DD*V(4) + TEMP(1, 2)
                       TEMP(1, 3) = EFACTP*DD*V(5) + TEMP(1, 3)
C                          
  220                  CONTINUE
                    ENDDO
  219               CONTINUE
                 ENDDO
C                
                 CALL GETTYPE (LB, MB, NB + 1, JDTYP) 
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C
                 DO JJQ = 1, ITM
                    IF (AAA(JJQ) .EQ. 0.0) GO TO 221
                    IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(JJQ), 2**10 - 1)
C
                    DO KKQ = 1, JTM
                       DD = AAA(JJQ)*BBB(KKQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 222
                       JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1)
     &                      + IL
                       JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1)
     &                      + IM
                       JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN
C
                       CALL SEC_MOMENT(JL, JM, JN, GAMAP, V, NT1,
     &                                 DISTQD)
C
                       TEMP(2, 1) = EFACTP*DD*V(1) + TEMP(2, 1)
                       TEMP(2, 2) = EFACTP*DD*V(4) + TEMP(2, 2)
                       TEMP(2, 3) = EFACTP*DD*V(5) + TEMP(2, 3)
C
 222                   CONTINUE
                     ENDDO
 221                 CONTINUE
                 ENDDO
C
C Now (x,y,z)zpx
C
                 CALL GETTYPE (LB - 1, MB, NB, JDTYP)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C     
                 DO JJQ = 1, ITM  
                    IF (AAA(JJQ) .EQ. 0.0) GO TO 223
                    IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(JJQ), 2**10 - 1)
C     
                    DO KKQ = 1, JTM
                       DD = AAA(JJQ)*BBB(KKQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 224
                       JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1)
     &                      + IL
                       JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1)
     &                      + IM
                       JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN
C
C Now everything is ready to do the actual evaluation of overlap integral
C     
                       CALL SEC_MOMENT(JL, JM, JN, GAMAP, V, NT1,
     &                                 DISTQD)
C     
                       TEMP(3, 1) = EFACTP*DD*V(5) + TEMP(3, 1)
                       TEMP(3, 2) = EFACTP*DD*V(6) + TEMP(3, 2)
                       TEMP(3, 2) = EFACTP*DD*V(3) + TEMP(3, 3)
C                          
  224                  CONTINUE
                    ENDDO
  223               CONTINUE
                 ENDDO
C
                 CALL GETTYPE (LB + 1, MB, NB, JDTYP)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C     
                 DO JJQ = 1, ITM 
                    IF (AAA(JJQ) .EQ. 0.0) GO TO 225
                    IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(JJQ), 2**10 - 1)
C     
                    DO KKQ = 1, JTM
                       DD = AAA(JJQ)*BBB(KKQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 226
                       JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1) 
     &                      + IL
                       JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1) 
     &                      + IM
                       JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN
C
C
                       CALL SEC_MOMENT(JL, JM, JN, GAMAP, V, NT1,
     &                                 DISTQD)
C     
                       TEMP(4, 1) = EFACTP*DD*V(5) + TEMP(4, 1)
                       TEMP(4, 2) = EFACTP*DD*V(6) + TEMP(4, 2)
                       TEMP(4, 3) = EFACTP*DD*V(3) + TEMP(4, 3)
C     
  226                  CONTINUE
                     ENDDO
  225                CONTINUE
                 ENDDO
C
C Now M(yx,yy,yz) = (-(xxpz - zxpx),-(yxpz-yzpx),-(zxpz-zzpx))
C
                 DV(4) = -((NB*TEMP(1, 1) - TWO*EXPS(JEX)*TEMP(2, 1))
     &                    -(LB*TEMP(3, 1) - TWO*EXPS(JEX)*TEMP(4, 1)))
                 DV(5) = -((NB*TEMP(1, 2) - TWO*EXPS(JEX)*TEMP(2, 2))
     &                    -(LB*TEMP(3, 2) - TWO*EXPS(JEX)*TEMP(4, 2)))
                 DV(6) = -((NB*TEMP(1, 3) - TWO*EXPS(JEX)*TEMP(2, 3))
     &                    -(LB*TEMP(3, 3) - TWO*EXPS(JEX)*TEMP(4, 3)))

C  l_z =  xPy - ypx, First (x,y,z)xpy
C
                 CALL DZERO(TEMP,12)

                 CALL GETTYPE (LB, MB - 1, NB, JDTYP)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C
                 DO JJQ = 1, ITM
                    IF (AAA(JJQ) .EQ. 0.0) GO TO 227
                    IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(JJQ), 2**10 - 1) 
C
                    DO KKQ = 1, JTM
                       DD = AAA(JJQ)*BBB(KKQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 228
                       JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1)
     &                      + IL
                       JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1)
     &                      + IM
                       JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN
C
C Now everything is ready to do the actual evaluation of overlap integral
C                
                       CALL DIP_MOMENT(JL, JM, JN, GAMAP, V, NT,
     &                                 DISTQD)
C                
                       TEMP(1, 1) = EFACTP*DD*V(1) + TEMP(1, 1)
                       TEMP(1, 2) = EFACTP*DD*V(4) + TEMP(1, 2)
                       TEMP(1, 3) = EFACTP*DD*V(5) + TEMP(1, 3)
C                          
  228                  CONTINUE
                    ENDDO
  227               CONTINUE
                 ENDDO
C                
                 CALL GETTYPE (LB, MB + 1, NB, JDTYP) 
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C
                 DO JJQ = 1, ITM
                    IF (AAA(JJQ) .EQ. 0.0) GO TO 229
                    IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(JJQ), 2**10 - 1)
C
                    DO KKQ = 1, JTM
                       DD = AAA(JJQ)*BBB(KKQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 230
                       JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1)
     &                      + IL
                       JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1)
     &                      + IM
                       JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN
C
                       CALL DIP_MOMENT(JL, JM, JN, GAMAP, V, NT,
     &                                 DISTQD)
C
                       TEMP(2, 1) = EFACTP*DD*V(1) + TEMP(2, 1)
                       TEMP(2, 2) = EFACTP*DD*V(4) + TEMP(2, 2)
                       TEMP(2, 3) = EFACTP*DD*V(5) + TEMP(2, 3)
C
  230                  CONTINUE
                     ENDDO
  229                CONTINUE
                 ENDDO
C
C Now (x,y,z)ypx
                 CALL GETTYPE (LB - 1, MB, NB, JDTYP)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C     
                 DO JJQ = 1, ITM  
                    IF (AAA(JJQ) .EQ. 0.0) GO TO 231
                    IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(JJQ), 2**10 - 1)
C     
                    DO KKQ = 1, JTM
                       DD = AAA(JJQ)*BBB(KKQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 232
                       JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1)
     &                      + IL
                       JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1)
     &                      + IM
                       JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN
C
C Now everything is ready to do the actual evaluation of overlap integral
C     
                       CALL SEC_MOMENT(JL, JM, JN, GAMAP, V, NT1,
     &                                 DISTQD)
C     
                       TEMP(3, 1) = EFACTP*DD*V(4) + TEMP(3, 1)
                       TEMP(3, 2) = EFACTP*DD*V(2) + TEMP(3, 2)
                       TEMP(3, 3) = EFACTP*DD*V(6) + TEMP(3, 3)
C                          
  232                  CONTINUE
                    ENDDO
  231               CONTINUE
                 ENDDO
C
                 CALL GETTYPE (LB + 1, MB, NB, JDTYP)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C     
                 DO JJQ = 1, ITM 
                    IF (AAA(JJQ) .EQ. 0.0) GO TO 233
                    IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(JJQ), 2**10 - 1)
C     
                    DO KKQ = 1, JTM
                       DD = AAA(JJQ)*BBB(KKQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 234
                       JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1) 
     &                      + IL
                       JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1) 
     &                      + IM
                       JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN
C
                       CALL SEC_MOMENT(JL, JM, JN, GAMAP, V, NT1,
     &                                 DISTQD)
C     
                       TEMP(4, 1) = EFACTP*DD*V(4) + TEMP(4, 1)
                       TEMP(4, 2) = EFACTP*DD*V(2) + TEMP(4, 2)
                       TEMP(4, 3) = EFACTP*DD*V(6) + TEMP(4, 3)
C     
  234                  CONTINUE
                     ENDDO
  233                CONTINUE
                 ENDDO
C
C Now M(zx,zy,zz) = (xxpy - xypz,yxpy - yypz,zypy - zypz)
C
                 DV(7) = (MB*TEMP(1, 1) - TWO*EXPS(JEX)*TEMP(2, 1)) - 
     &                   (LB*TEMP(3, 1) - TWO*EXPS(JEX)*TEMP(4, 1))
                 DV(8) = (MB*TEMP(1, 2) - TWO*EXPS(JEX)*TEMP(2, 2)) - 
     &                   (LB*TEMP(3, 2) - TWO*EXPS(JEX)*TEMP(4, 2))
                 DV(9) = (MB*TEMP(1, 3) - TWO*EXPS(JEX)*TEMP(2, 3)) - 
     &                   (LB*TEMP(3, 3) - TWO*EXPS(JEX)*TEMP(4, 3))
C
#ifdef _DEBUG_LVLM1
                 write(luout,*)
                 write(luout,*) 'Center1  ', i ,'     ',ITYP
                 write(luout,*) 'Center2  ', j ,'     ',JTYP
                 Write (luout, 999) (DV(IIII), IIII = 1 ,3)
 999             format(5X,F12.5, 5X,F12.5,5X,F12.5) 
#endif
                 DO 30 MM = 1, 9
                    VS(IJIND, MM) = DV(MM) + VS(IJIND, MM)
 30              CONTINUE
C
#ifdef _DEBUG_LVLM1
                 write(luout,*)
                 write(luout,*) 'Center1  ', i ,'     ',ITYP
                 write(luout,*) 'Center2  ', j ,'     ',JTYP
                 Write (luout, 999) (VS(IJIND, IIII), IIII = 1 ,10)
#endif
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
C
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
                  IF(MTYP .EQ. 24) XP = -XP
                  IJF = IS*(IS-1)/2+JS
                  AOINT(IJF)=XP+AOINT(IJF)
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
c     
#ifdef _DEBUG_LVLM
         Write(6,*) "Triangularly Packed AO Ang. Mom. integrals"
         Write(6,"(5F12.8)") (aoint(i), i=1, nu2)
#endif
         IF (MTYP .EQ. 4 .OR. MTYP .EQ. 24) THEN
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
