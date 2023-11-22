C
      SUBROUTINE GENQSD(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF)
C
C Driver for Q(ALPHA, BETA) and Q(ALPHA, BETA) one-body and
C             AE                 IJ,AB
C two body contributions from I(ALPHA) and G(ALPHA) intermediates.
C The algebraic expressions are given in actual subroutine which
C evaluates particular contributions.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION ICORE(MAXCOR)
      INTEGER POP, VRT, DIRPRD
C
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT1AA,NT1BB,NF1AA,NF1BB,
     &             NF2AA,NF2BB
      COMMON /MACHSP/ IINTLN, IFLTLN, IINTFP, IALONE, IBITWD
      COMMON /SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
C
C Common blocks used in the quadratic term
C
      COMMON /QADINTI/ INTIMI, INTIAE, INTIME
      COMMON /QADINTG/ INGMNAA, INGMNBB, INGMNAB, INGABAA,
     &                 INGABBB, INGABAB, INGMCAA, INGMCBB,
     &                 INGMCABAB, INGMCBABA, INGMCBAAB,
     &                 INGMCABBA
      COMMON /APERTT1/ IAPRT1AA
      COMMON /BPERTT1/ IBPRT1AA
      COMMON /APERTT2/ IAPRT2AA1, IAPRT2BB1, IAPRT2AB1, IAPRT2AB2, 
     &                 IAPRT2AB3, IAPRT2AB4, IAPRT2AA2, IAPRT2BB2,
     &                 IAPRT2AB5 
C
C Genrate one-body and two-body contributions from the intermediates.
C
C
      IF (IFLAGS(1) .GE. 10) THEN
         CALL HEADER('Calculating one and two body contributions', 
     &                0, LUOUT)
         CALL TIMER(1)
      ENDIF
C
C Zero out the lists to which the contributions are being written
C
      DO 10 ISPIN = 3, 3 - 2*IUHF, -1
         CALL ZEROLIST(ICORE, MAXCOR, 60 + ISPIN)
         CALL ZEROLIST(ICORE, MAXCOR, 40 + ISPIN)
 10   CONTINUE
C
      IF (IUHF .EQ. 0) THEN
         CALL ZEROLIST(ICORE, MAXCOR, 42)
      ELSE
         CALL ZEROLIST(ICORE, MAXCOR, 40)
      ENDIF
C
      DO 20 ISPIN = 1, IUHF + 1
         NFVO = IRPDPD(1, 8 + ISPIN)
         CALL IZERO(ICORE, NFVO*IINTFP)
         CALL PUTLST(ICORE, 1, 1, 1, 2 + ISPIN, 90)
 20   CONTINUE
C
C Create - P_(AB) SUM T(ALP) [1/2(SUM T(BET)*Hbar(AM,EF)]
C                 M    B,M        EF   (EF,IJ)
C to the doubles equation.
C
C
      CALL GABEFIND(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF)
C
C      
C Create + 1/2 SUM T(BET) * G(ALP)  hole-hole ladder contribution to 
C              M,N  (MN,AB)  (MN,IJ)
C to the doubles equation.
C
C
      CALL GMNIJIND(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF)
C
C                                E        A
C Calculate - P_(IJ)P(_(AB) SUM T(ALP) * T(BET) * Hbar(MB,EJ)
C                           M,E  I        M
C to the doubles equation.
C
C
      CALL T1T1IND1(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF)
C
C
C Create + P_(AB)P_(IJ) SUM T(BET) * G(ALP)  ring contribution to 
C                       M,E  (IM,AE)  (MB,EJ)
C to the doubles equation.
C
C
      IF (IUHF .EQ. 1) THEN 
         IBEGIN = 1
         IEND   = 3
      ELSE IF (IUHF. EQ. 0) THEN
         IBEGIN = 3
         IEND   = 3
      ENDIF
C     
      DO 30 ISPIN = IBEGIN, IEND
         CALL GMBEJIND(ICORE, MAXCOR, IRREPX, IOFFSET, ISPIN, IUHF)
 30   CONTINUE
C
C    
C Create P_(AB) SUM T(BET)  * I(ALP)  contribution to the doubles equation.
C               E   (IJ,AE)   (BE)
C In addition to that this routine also compute the SUM T(BET) * I(ALP)
C                                                   E    (I,E)    (AE)
C to the singles equation.
C
      CALL IAEINSD(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF)
C
C
C Create - P_(IJ) SUM T(BET)  * I(ALP)  contribution to the doubles equation.
C                 M   (IM,AB)   MJ
C In addition to that this routine also compute the - SUM T(BET) * I(ALP)
C                                                     M   (M,A)    (MI)
C to the singles equation.
C
C
      CALL IMIINSD(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF)
C
C
C Create SUM T(BET)  * I(ALP)  contribution to the singles equation.
C        M,E (IM,AE)   (ME)
C Note that there are no I(ALP) contribution to the doubles equation
C                        (ME)
C in the quadratic term.
C
C
      IF (IUHF .NE. 0) THEN
         CALL IMEINS(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF, 1)
         CALL IMEINS(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF, 2)
      ELSE IF (IUHF .EQ. 0) THEN
         CALL IMEINS(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF, 1)   
      ENDIF
C
C                          A        B
C Create + 1/2 P_(AB) SUM T(ALP) * T(BET) * Hbar(MN, IJ) and 
C                     M,N  M        N 
C                   E        F
C + 1/2 P_(AB) SUM T(ALP) * T(BET) * Hbar(AB, EF)  to the
C              E,F  I        J
C doubles equation.  The option to use AO ladder on 05/2019
C
C     
      IF (IFLAGS(93) .EQ. 2) THEN
         CALL T1T1IND2(ICORE, MAXCOR, IRREPX, IOFFSET, 1, IUHF)
         CALL T1T1IND2_AO(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF)
      ELSE
         DO 40 ITYPE = 1, 6, 5
            CALL T1T1IND2(ICORE, MAXCOR, IRREPX, IOFFSET, ITYPE, IUHF)
 40      CONTINUE
      ENDIF 
C 
C                       E        F        AB
C Create - P_(IJ) SUM  T(ALP) * T(BET) * T  * Hbar(MN,EF)
C                 M,E   M        I        NJ
C                 N,F
C                   E        A        FB
C and - P_(AB) SUM T(ALP) * T(BET) * T  * Hbar(MN,EF)
C              M,E  M        N        IJ
C              N,F
C Three body contributions to the doubles equation
C
      CALL QTHRBDY1(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF)
      CALL QTHRBDY2(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF)
C
C This will finish the evaluation of contribution to singles 
C and doubles equations from G(ALPHA) and I(ALPHA) intermediates.
C Write the time taken if requested.
C
      IF (IFLAGS(1) .GE. 10) THEN
         CALL TIMER (1)
         WRITE(LUOUT, *)
         WRITE(LUOUT, 101) TIMENEW
 101     FORMAT(T3, '@GENQSD-I, One body and two body terms required',
     &          F9.3, ' seconds.')
      ENDIF
C
      RETURN
      END
