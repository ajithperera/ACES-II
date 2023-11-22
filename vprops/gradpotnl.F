C
      SUBROUTINE GradPotnl(VLIST,NOC,CENTR1,EXPS,CONT,LABEL,NPRIM,
     &                     NCONTR,ITRM,ITRN,CTRN,NCONT,BUF,IBUF,
     &                     AOINT,LBUF,INDX,NGROUP,VS,SCR,CHARGE,
     &                     MXPRM2)
C
C Derivative of the nuclear attaraction term.
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
      DIMENSION IHOLD(10),IADH(10),VS(MXPRM2,10),TEMP(4,1)
      DIMENSION BUF(LBUF,10),IBUF(IINTFP*LBUF,10),AOINT(1)
      DIMENSION TITLE(10,25)
      DIMENSION CENTR1(3),SCR(1),GBUF(600),IGNDX(600)
      DIMENSION LMN(200),JMN(200)
      DIMENSION A(3),B(3),EP(3),DISTQD(3),AAA(200),BBB(200),
     &          V(10),DV(3)
C     
      DATA (LABELJA(I, 23), I = 1, 10) /'   PVP  ',9*'        '/
      DATA (TITLE(I, 23), I = 1, 10)/'   PVP  ',   9*'        '/
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
      FOUR = 4.00D+00
C
C Do some initializations (No nuclear contribution for this term)
C
      CALL SETRHF
      NT = 1
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

        DO 410 II = 1, 3
           A(II) = VLIST(ICNT, II)
 410    CONTINUE
C
        IE = NPRIM(I)
        ICE = NCONTR(I)
        JEXH = 0
        JCONT = 0
        JCOFF = 0
     
        DO 49 J = 1, I
C     
           JCNT = IBTAND(IBTSHR(LABEL(J),10), 2**10 - 1)
           JTYP = IBTAND(LABEL(J), 2**10-1)
C
           DO 411 JJ = 1, 3
              B(JJ) = VLIST(JCNT, JJ)
 411       CONTINUE
C
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
C
                 CALL CIVPT(A, EXPS(IEX), B, EXPS(JEX), EP, GAMAP,
     &                      EFACTP)
C     
C Get the angular momentum quantum numbers of the gaussian
C centered on A and B

                 LA = LMNVAL(1, ITYP)
                 MA = LMNVAL(2, ITYP)
                 NA = LMNVAL(3, ITYP)
C
                 LB = LMNVAL(1, JTYP)
                 MB = LMNVAL(2, JTYP)
                 NB = LMNVAL(3, JTYP)
#ifdef _DEBUG_LVM1
                 Write(6,*)
                 Write(6, "(a,6I2)"), "la,ma,na,lb,mb,nb: ",
     &                                 la,ma,na,lb,mb,nb
                 Write(6,*)
#endif
C
C Lets form the three intermediates resulting from the derivatives
C of the Gaussian on B. First with respect to X
C
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
                 CALL GETTYPE (LA - 1, MA, NA, IDTYP)
                 CALL GETTYPE (LB - 1, MB, NB, JDTYP)
                 CALL RHFTCE (AAA, A, EP, IDTYP, ITM, ONE, LMN)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C     
                 DO JJA = 1, 4
                    TEMP(JJA, 1) = ZERO
                 ENDDO

                 DO JJQ = 1, ITM  
                    IF (AAA(JJQ) .EQ. 0.0) GO TO 210
                    IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(JJQ), 2**10 - 1)
C     
                    DO KKQ = 1, JTM
                       DD = AAA(JJQ)*BBB(KKQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 211
                       JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1)
     &                      + IL
                       JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1)
     &                      + IM
                       JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN
C
C Now everything is ready to do the actual evaluation of repulsion integral
C     
                       CALL ATTRACT(JL, JM, JN, GAMAP ,V, DISTQD)
C     
                       TEMP(1, 1) = EFACTP*DD*V(1) + TEMP(1, 1)
  211                  CONTINUE
                    ENDDO
  210               CONTINUE
                 ENDDO
C
#ifdef _DENBUG_LVM1
                 Write(6,"(a,F10.5)")  "dx,la-1,lb-1: ", TEMP(1,1)
#endif
C
                 CALL GETTYPE (LB + 1, MB, NB, JDTYP)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C     
                 DO JJQ = 1, ITM 
                    IF (AAA(JJQ) .EQ. 0.0) GO TO 212
                    IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(JJQ), 2**10 - 1)
C     
                    DO KKQ = 1, JTM
                       DD = AAA(JJQ)*BBB(KKQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 213
                       JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1) 
     &                      + IL
                       JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1) 
     &                      + IM
                       JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN
C
                       CALL ATTRACT(JL, JM, JN, GAMAP, V, DISTQD)
C     
                       TEMP(2, 1) = EFACTP*DD*V(1) + TEMP(2, 1)
  213                  CONTINUE
                    ENDDO
  212               CONTINUE
                 ENDDO
C
#ifdef _DEBUG_LVM1
                 Write(6,"(a,F10.5)")  "dx,la-1,lb+1: ",TEMP(2,1)
#endif
C
                 CALL GETTYPE (LA + 1, MA, NA, IDTYP)
                 CALL GETTYPE (LB - 1, MB, NB, JDTYP)
                 CALL RHFTCE (AAA, A, EP, IDTYP, ITM, ONE, LMN)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C
                 DO JJQ = 1, ITM
                    IF (AAA(JJQ) .EQ. 0.0) GO TO 214
                    IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(JJQ), 2**10 - 1)
C
                    DO KKQ = 1, JTM
                       DD = AAA(JJQ)*BBB(KKQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 215
                       JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1)
     &                      + IL
                       JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1)
     &                      + IM
                       JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN
C
                       CALL ATTRACT(JL, JM, JN, GAMAP, V, DISTQD)
C
                       TEMP(3, 1) = EFACTP*DD*V(1) + TEMP(3, 1)
C
  215                  CONTINUE
                     ENDDO
  214                CONTINUE
                 ENDDO 
C
#ifdef _DEBUG_LVM1
                 Write(6,"(a,F10.5)") "dx,la-1,lb+1: ", TEMP(3,1)
#endif
C
                 CALL GETTYPE (LA + 1, MA , NA, IDTYP)
                 CALL GETTYPE (LB + 1, MB , NB, JDTYP)
                 CALL RHFTCE (AAA, A, EP, IDTYP, ITM, ONE, LMN)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C
                 DO JJQ = 1, ITM
                    IF (AAA(JJQ) .EQ. 0.0) GO TO 216
                    IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(JJQ), 2**10 - 1)
C
                    DO KKQ = 1, JTM
                       DD = AAA(JJQ)*BBB(KKQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 217
                       JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1)
     &                      + IL
                       JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1)
     &                      + IM
                       JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN
C
                       CALL ATTRACT(JL, JM, JN, GAMAP, V, DISTQD)
C
                       TEMP(4, 1) = EFACTP*DD*V(1) + TEMP(4, 1)
C
  217                  CONTINUE
                   ENDDO
  216              CONTINUE
                 ENDDO
C
#ifdef _DEBUG_LVM1
                 Write(6,"(a,F10.5)") "dx,la+1,lb+1: ",TEMP(4,1)
#endif
C
C Now evaluate the contribution to the X-derivative
C
                 DV(1) = LB*LA*TEMP(1, 1)-TWO*LA*EXPS(JEX)*TEMP(2, 1) 
     &                   -TWO*LB*EXPS(IEX)*TEMP(3, 1)+FOUR*EXPS(JEX)*
     &                   EXPS(IEX)*TEMP(4, 1)
C
#ifdef _DEBUG_LVM1
                 Write(6,*)
                 Write(6,"(a,F10.5)")  "DV(1) = ", DV(1)
#endif
C
C Evalute the derivative intermediate with respect to Y
C     
                 CALL GETTYPE (LA , MA - 1, NA, IDTYP)
                 CALL GETTYPE (LB , MB - 1, NB, JDTYP)
                 CALL RHFTCE (AAA, A, EP, IDTYP, ITM, ONE, LMN)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C
                 DO  JJA = 1, 4
                    TEMP(JJA, 1) = ZERO
                 ENDDO 
C
                 DO JJQ = 1, ITM
                    IF (AAA(JJQ) .EQ. 0.0) GO TO 218
                    IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(JJQ), 2**10 - 1)
C
                    DO KKQ = 1, JTM
                       DD = AAA(JJQ)*BBB(KKQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 219
                       JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1)
     &                      + IL
                       JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1)
     &                      + IM
                       JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN
C
C Now everything is ready to do the actual evaluation of overlap integral
C
                       CALL ATTRACT(JL, JM, JN, GAMAP ,V, DISTQD)
C
                       TEMP(1, 1) = EFACTP*DD*V(1) + TEMP(1, 1)
  219                  CONTINUE
                    ENDDO
  218            CONTINUE
                 ENDDO
C                
#ifdef _DEBUG_LVM1
                 Write(6,*)
                 Write(6,"(a,F10.5)")  "dy,ma-1,mb-1: ", TEMP(1,1)
#endif
C
                 CALL GETTYPE (LB , MB + 1, NB, JDTYP)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C
                 DO JJQ = 1, ITM
                    IF (AAA(JJQ) .EQ. 0.0) GO TO 220
                    IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(JJQ), 2**10 - 1)
C
                    DO KKQ = 1, JTM
                       DD = AAA(JJQ)*BBB(KKQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 221
                       JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1)
     &                      + IL
                       JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1)
     &                      + IM
                       JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN
C
                       CALL ATTRACT(JL, JM, JN, GAMAP, V, DISTQD)
C
                       TEMP(2, 1) = EFACTP*DD*V(1) + TEMP(2, 1)
  221                  CONTINUE
                    ENDDO
  220            CONTINUE
                ENDDO
C
#ifdef _DEBUG_LVM1
                 Write(6,"(a,F10.5)")  "dy,ma-1,mb+1: ", TEMP(2,1)
#endif
C
                 CALL GETTYPE (LA , MA + 1, NA, IDTYP)
                 CALL GETTYPE (LB , MB - 1, NB, JDTYP)
                 CALL RHFTCE (AAA, A, EP, IDTYP, ITM, ONE, LMN)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)

                 DO JJQ = 1, ITM
C...watson...added the if AAA go to 222
                    IF (AAA(JJQ) .EQ. 0.0) GO TO 222
                    IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(JJQ), 2**10 - 1)
C
                    DO KKQ = 1, JTM
                       DD = AAA(JJQ)*BBB(KKQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 223
                       JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1)
     &                      + IL
                       JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1)
     &                      + IM
                       JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN

C
                       CALL ATTRACT(JL, JM, JN, GAMAP, V, DISTQD)
C
                       TEMP(3, 1) = EFACTP*DD*V(1) + TEMP(3, 1)
  223                  CONTINUE
                     ENDDO
  222            CONTINUE
                 ENDDO
C
#ifdef _DEBUG_LVM1
                 Write(6,"(a,F10.5)")  "dy,ma+1,mb-1: ", TEMP(3,1)
#endif
C
                 CALL GETTYPE (LA, MA + 1, NA, IDTYP)
                 CALL GETTYPE (LB, MB + 1, NB, JDTYP)
                 CALL RHFTCE (AAA, A, EP, IDTYP, ITM, ONE, LMN)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C
                 DO JJQ = 1, ITM
                    IF (AAA(JJQ) .EQ. 0.0) GO TO 224
                    IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(JJQ), 2**10 - 1)
C
                    DO KKQ = 1, JTM
                       DD = AAA(JJQ)*BBB(KKQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 225
                       JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1)
     &                      + IL
                       JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1)
     &                      + IM
                       JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN
C
                       CALL ATTRACT(JL, JM, JN, GAMAP, V, DISTQD)
C
                       TEMP(4, 1) = EFACTP*DD*V(1) + TEMP(4, 1)
 225                   CONTINUE
                   ENDDO
 224             CONTINUE
                 ENDDO
C
#ifdef _DEBUG_LVM1
                 Write(6,"(a,F10.5)")  "dy,ma+1,mb+1: ", TEMP(4,1)
#endif

C
C eValute the derivative intermediate with respect to Y
C     
                 DV(2) = MA*MB*TEMP(1, 1)-TWO*MA*EXPS(JEX)*TEMP(2, 1) 
     &                   -TWO*MB*EXPS(IEX)*TEMP(3, 1)+FOUR*EXPS(JEX)*
     &                   EXPS(IEX)*TEMP(4, 1)
C
#ifdef _DEBUG_LVM1
                 Write(6,*)
                 Write(6,"(a,F10.5)")  "DV(2) = ", DV(2)
                 Write(6,*)
#endif
C
C Evaluate the derivative intermediate with respect to Z
C
                 CALL GETTYPE (LA , MA , NA - 1, IDTYP)
                 CALL GETTYPE (LB , MB , NB - 1, JDTYP)
                 CALL RHFTCE (AAA, A, EP, IDTYP, ITM, ONE, LMN)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C
                 DO JJA = 1, 4
                    TEMP(JJA, 1) = ZERO
                 ENDDO
C
                 DO JJQ = 1, ITM
                    IF (AAA(JJQ) .EQ. 0.0) GO TO 226
                    IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(JJQ), 2**10 - 1)
C
                    DO KKQ = 1, JTM
                       DD = AAA(JJQ)*BBB(KKQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 227
                       JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1)
     &                      + IL
                       JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1)
     &                      + IM
                       JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN
C
C Now everything is ready to do the actual evaluation of overlap integral
C
                       CALL ATTRACT(JL, JM, JN, GAMAP ,V, DISTQD)
C
                       TEMP(1, 1) = EFACTP*DD*V(1) + TEMP(1, 1)
  227                  CONTINUE
                    ENDDO
  226            CONTINUE
                 ENDDO
C
#ifdef _DEBUG_LVM1
                 Write(6,"(a,F10.5)")  "dz,na-1,nb-1: ", TEMP(1,1)
#endif
C
                 CALL GETTYPE (LB, MB, NB + 1, JDTYP)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C
                 DO JJQ = 1, ITM
                    IF (AAA(JJQ) .EQ. 0.0) GO TO 228
                    IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(JJQ), 2**10 - 1)
C
                    DO KKQ = 1, JTM
                       DD = AAA(JJQ)*BBB(KKQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 229
                       JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1)
     &                      + IL
                       JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1)
     &                      + IM
                       JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN
C
                       CALL ATTRACT(JL, JM, JN, GAMAP, V, DISTQD)
C
                       TEMP(2, 1) = EFACTP*DD*V(1) + TEMP(2, 1)
  229                  CONTINUE
                    ENDDO
  228            CONTINUE
                 ENDDO
C
#ifdef _DEBUG_LVM1
                 Write(6,"(a,F10.5)")  "dz,na-1,nb+1: ", TEMP(2,1)
#endif
C
                 CALL GETTYPE (LA, MA, NA + 1, IDTYP)
                 CALL GETTYPE (LB, MB, NB - 1, JDTYP)
                 CALL RHFTCE (AAA, A, EP, IDTYP, ITM, ONE, LMN)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C                
                 DO JJQ = 1, ITM
                    IF (AAA(JJQ) .EQ. 0.0) GO TO 230
                    IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(JJQ), 2**10 - 1)
C                
                    DO KKQ = 1, JTM
                       DD = AAA(JJQ)*BBB(KKQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 231
                       JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1)
     &                      + IL
                       JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1)
     &                      + IM
c...watson... added the + IM
                       JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN
C                   
                       CALL ATTRACT(JL, JM, JN, GAMAP, V, DISTQD)
C                   
                       TEMP(3, 1) = EFACTP*DD*V(1) + TEMP(3, 1)
  231                  CONTINUE
                    ENDDO
  230            CONTINUE
                 ENDDO
C
#ifdef _DEBUG_LVM1
                 Write(6,"(a,F10.5)")  "dz,na+1,nb-1: ", TEMP(3,1)
#endif
C 
                 CALL GETTYPE (LA, MA, NA + 1, IDTYP)
                 CALL GETTYPE (LB, MB, NB + 1, JDTYP)
                 CALL RHFTCE (AAA, A, EP, IDTYP, ITM, ONE, LMN)
                 CALL RHFTCE (BBB, B, EP, JDTYP, JTM, ONE, JMN)
C
                 DO JJQ = 1, ITM
                    IF (AAA(JJQ) .EQ. 0.0) GO TO 232
                    IL = IBTAND(IBTSHR(LMN(JJQ),20), 2**10 - 1)
                    IM = IBTAND(IBTSHR(LMN(JJQ),10), 2**10 - 1)
                    IN = IBTAND(LMN(JJQ), 2**10 - 1)
C
                    DO KKQ = 1, JTM
                       DD = AAA(JJQ)*BBB(KKQ)
                       IF (ABS(DD) .LT. THR2)  GO TO 233
                       JL = IBTAND(IBTSHR(JMN(KKQ),20), 2**10 - 1)
     &                      + IL
                       JM = IBTAND(IBTSHR(JMN(KKQ),10), 2**10 - 1)
     &                      + IM
                       JN = IBTAND(JMN(KKQ), 2**10 - 1) + IN
C
                       CALL ATTRACT(JL, JM, JN, GAMAP, V, DISTQD)
C
                       TEMP(4, 1) = EFACTP*DD*V(1) + TEMP(4, 1)
  233                 CONTINUE
                    ENDDO
  232            CONTINUE
                 ENDDO
C
#ifdef _DEBUG_LVM1
                 Write(6,"(a,F10.5)")  "dz,na+1,nb+1: ", TEMP(4,1)
#endif
C
C
C Now evaluate the contribution to the Z-derivative
C
                 DV(3) = NA*NB*TEMP(1, 1)-TWO*NA*EXPS(JEX)*TEMP(2, 1) 
     &                   -TWO*NB*EXPS(IEX)*TEMP(3, 1)+FOUR*EXPS(JEX)*
     &                   EXPS(IEX)*TEMP(4, 1)
C
#ifdef _DEBUG_LVM1
                 Write(6,*)
                 Write(6,"(a,F10.5)")  "DV(3) = ", DV(3)
                 write(luout,*)
                 write(luout,*) 'Center1  ', i ,'     ',ITYP
                 write(luout,*) 'Center2  ', j ,'     ',JTYP
                 Write(luout,*)
                 Write (luout, 999) (DV(IIII), IIII = 1 ,3)
 999             format(5X,F12.5, 5X,F12.5,5X,F12.5) 
#endif
C 
C Note that pVp integral must have a negative sign (1/i x 1/i).
C
                 DO 30 MM = 1, 3
                    VS(IJIND, 1) = CHARGE*DV(MM) + VS(IJIND, 1)
 30              CONTINUE
#ifdef _DEBUG_LVM1
C                 write(luout,*)
C                 write(luout,*) 'Center1  ', ii ,'     ',IE
C                 write(luout,*) 'Center2  ', jj ,'     ',JE
C                 Write (luout, 999) (VS(IJIND, IIII), IIII = 1 ,1)
#endif
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
C
#ifdef _DEBUG_LVM1
              Write(6,*)
              Write(6,"(a)") "After the contraction"
              Write(6,*)
              Write (luout, 999) (VS(1,IIII), IIII = 1 ,3)
              Write(6,"(a)") "--------------------------"
              write(6,*)
#endif
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
                  IF(MTYP .EQ. 22) XP = -XP
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
C     
*******DEBUG
C         call output(aoint, 1, nsabf, 1, nsabf, nsabf, nsabf, 1)
C          Write(6, "(3(1x,E20.5))") (aoint(i), i=1, nu2)
C&&&&Debug
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
