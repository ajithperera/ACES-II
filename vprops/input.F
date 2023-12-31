      SUBROUTINE INPUT(EXPS,COEF,IT,NCONTC,ITRM,ITRN,CTRN,NOC,
     1 NCONT,INU,IOUTU,VLIST,CENTR,CENTS,CHARGS,IMOL,IPFLG,
     2 NPRIM, NCONTR, INDX, LABEL, IDONE, EXPNNT, CONTRC, NGROUP,
     3 NOADC, JOBIPH, INTERF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER*8 EXCLAM,NUCDAT,BLANKS
#ifdef _DOUG_KROLL
      DIMENSION VLIST(NOC,11)
#else
      DIMENSION VLIST(NOC,10)
#endif
      DIMENSION CONTRC(1),IT(1),CHARGS(1),NCONTC(1)
c The hardcoded limit 256 on ITRN and CTRN corresponds to the number of
c primitives in one contracted basis function. I think 256 is safe since
c it is doubtful that we will have that many primitives in one contraction.
c Ajith Perera, 05/2005.
      DIMENSION ITRAN(256), CTRAN(256)
      DIMENSION EXPNNT(1), EXPS(1), COEF(1)
      DIMENSION NPRIM(1), NCONTR(1), INDX(1), IDONE(1), LABEL(1)
      DIMENSION ITRM(1)
      DIMENSION ITRN(256,1),CTRN(256,1)
      COMMON /COMINF/ TOTMSS, COM(3), RMOI(3), PATRAN(3,3)
      COMMON /VARS/ FNU(20),DP(20),FN(15),FD(15),PI,DFTR(11)
      COMMON /HIGHL/ LMNVAL(3,84), ANORM(84)
      COMMON /PRTDAT/ IDOSPH, MAXAQN, NSABF
      DIMENSION MAP(84)
      DATA MAP/1,3*2,6*3,10*4,15*5,21*6,28*7/
C SG 8/25/98 CENTR changed to be character*4
      CHARACTER*4 ITXX, TYPE(84), IF3, IIF, IFX
      CHARACTER*4 CENTR(NOC), CENTS(NOC)
      DATA (TYPE(I),I = 1,56)
     O         /'S   ', 'X   ', 'Y   ', 'Z   ',
     1          'XX  ', 'XY  ', 'XZ  ', 'YY  ',
     2          'YZ  ', 'ZZ  ', 'F300', 'F210',
     3          'F201', 'F120', 'F111', 'F102',
     4          'F030', 'F021', 'F012', 'F003',
     5          'G400', 'G310', 'G301', 'G220',
     6          'G211', 'G202', 'G130', 'G121',
     7          'G112', 'G103', 'G040', 'G031',
     8          'G022', 'G013', 'G004', 'H500',
     9          'H410', 'H401', 'H320', 'H311',
     O          'H302', 'H230', 'H221', 'H212',
     1          'H203', 'H140', 'H131', 'H122',
     2          'H113', 'H104', 'H050', 'H041',
     3          'H032', 'H023', 'H014', 'H005'/
      DATA (TYPE(I),I = 57,84) /
     1 'I600', 'I510', 'I501', 'I420', 'I411', 'I402',
     2 'I330', 'I321', 'I312', 'I303', 'I240', 'I231',
     3 'I222', 'I213', 'I204', 'I150', 'I141', 'I132',
     4 'I123', 'I114', 'I105', 'I060', 'I051', 'I042',
     5 'I033', 'I024', 'I015', 'I006'/
      DATA EXCLAM /'!!!!!!!!'/, BLANKS /'        '/, 
     &     NUCDAT /'NUCLEI  '/
      IBTAND(I,J) = IAND(I,J)
      IBTOR(I,J)  = IOR(I,J)
      IBTXOR(I,J) = IEOR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTNOT(I)   = NOT(I)
      READ(INU ) (NCONTC(I),I=1,NCONT)
C      WRITE(6,110)
C110   FORMAT(1X,'GAUSSIAN  FUNCTION  COMPONENT   CENTER   TYPE        E
C     1PONENT     COEFFICIENT')
      IEX=0
      ICN=0
      DO 21 I=1,NCONT
      ICN=ICN+1
      IEXS=IEX
      IEX=IEX+1
      READ(INU ) IFX,K,ITXX, EXPNNT(IEX),CONTRC(IEX)
      CALL PUTLN(IFX,K)
C      IF(ABS(CONTRC(IEX)).LT.1.E-09) CONTRC(IEX) = 1.0
C      IF(IPFLG.NE.0)WRITE(6,111) IEX,I,1,IFX,ITXX,
C     1 EXPNNT(IEX),CONTRC(IEX)
      DO 4 J=1,NOC
      ICC=J
      IF(IFX.EQ.CENTR(J))GO TO 5
4     CONTINUE
      WRITE(6,300) I, IFX
300   FORMAT(' INVALID CENTER   FUNCTION ',I4,'  CENTER ',A8)
      CALL ERREX
5     CONTINUE
      DO 7 J=1,84
      MMM = J
      IF(TYPE(J).EQ.ITXX)GO TO 8
7     CONTINUE
      WRITE(6,301) I, ITXX
301   FORMAT(' INVALID TYPE     FUNCTION ',I4,'  CENTER ',A8)
      CALL ERREX
8     CONTINUE
#ifdef _UNICOS
c      IF (MMM .EQ. 1) CALL AUGCNT(VLIST(ICC,9), VLIST(ICC,10), 1)
c      IF (MMM .EQ. 2) CALL AUGCNT(VLIST(ICC,9), VLIST(ICC,10), 2)
c      IF (MMM .EQ. 5) CALL AUGCNT(VLIST(ICC,9), VLIST(ICC,10), 3)
c      IF (MMM .EQ.11) CALL AUGCNT(VLIST(ICC,9), VLIST(ICC,10), 4)
c      IF (MMM .EQ.21) CALL AUGCNT(VLIST(ICC,9), VLIST(ICC,10), 5)
c      IF (MMM .EQ.36) CALL AUGCNT(VLIST(ICC,9), VLIST(ICC,10), 6)
c      IF (MMM .EQ.57) CALL AUGCNT(VLIST(ICC,9), VLIST(ICC,10), 7)
#else
      IF (MMM .EQ. 1) CALL AUGCNT(VLIST(ICC,9), VLIST(ICC,10), 1)
      IF (MMM .EQ. 2) CALL AUGCNT(VLIST(ICC,9), VLIST(ICC,10), 2)
      IF (MMM .EQ. 5) CALL AUGCNT(VLIST(ICC,9), VLIST(ICC,10), 3)
      IF (MMM .EQ.11) CALL AUGCNT(VLIST(ICC,9), VLIST(ICC,10), 4)
      IF (MMM .EQ.21) CALL AUGCNT(VLIST(ICC,9), VLIST(ICC,10), 5)
      IF (MMM .EQ.36) CALL AUGCNT(VLIST(ICC,9), VLIST(ICC,10), 6)
      IF (MMM .EQ.57) CALL AUGCNT(VLIST(ICC,9), VLIST(ICC,10), 7)
#endif /* _UNICOS */
c      IT(ICN)=IBTOR(IBTSHL(ICC,20),MMM)
      IT(ICN)=IBTOR(IBTSHL(ICC,10),MMM)
      JE=NCONTC(I)
      IF(JE.EQ.1)GO TO 3
      DO 34 J=2,JE
      IEX=IEX+1
      READ(INU ) IIF,IF2,IF3, EXPNNT(IEX),CONTRC(IEX)
C      IF(ABS(CONTRC(IEX)).LT.1.E-09) CONTRC(IEX) = 1.0
C      IF(IPFLG.NE.0)WRITE(6,111) IEX,I,J,IFX,ITXX,
C     1 EXPNNT(IEX),CONTRC(IEX)
34    CONTINUE
3     CONTINUE
21    CONTINUE
C
C....    NOW RESTRUCTURE BASIS TO ALLOW GENERAL CONTRACTION TO BE
C....    PERFORMED WITHOUT DUPLICATING PRIMITIVES.
C
      ICNTR = 0
      IEXOLD = 0
      ICOEF = 0
      IEX = 0
      IGRP = 0
      DO 600 I = 1,NCONT
      IF (IDONE(I) .NE. 0) GOTO 610
C
C....    THIS FUNCTION HAS NOT YET BEEN PROCESSED
C
      IDONE(I) = 1
      ICNTR = ICNTR + 1
      INDX(ICNTR) = I
      IGRP = IGRP + 1
      LABEL(IGRP) = IT(I)
      NPRIMI = NCONTC(I)
      NPRIM(IGRP) = NPRIMI
      NCONTI = 1
      DO 615 J = 1, NPRIMI
      EXPS(IEX+J) = EXPNNT(IEXOLD+J)
      COEF(ICOEF+J) = CONTRC(IEXOLD+J)
615   CONTINUE
      ICOEF = ICOEF + NPRIMI
C
C....    NOW LOOK FOR OTHER CONTRACTED FUNCTIONS BUILT FROM THESE PRIMIT
C
      IF (I .EQ. NCONT) GOTO 625
      IEXOLE = IEXOLD + NPRIMI
      DO 620 J = I+1, NCONT
      IF (LABEL(IGRP) .NE. IT(J) .OR.
     1    NPRIMI .NE. NCONTC(J) .OR.
     2    ABS(EXPS(IEX+1) - EXPNNT(IEXOLE+1)) .GT. 0.0001 .OR.
     3    IDONE(J) .NE. 0) GOTO 630
C
C....    MATCH
C
      IDONE(J) = 1
      ICNTR = ICNTR + 1
      INDX(ICNTR) = J
      NCONTI = NCONTI + 1
      DO 640 K = 1, NPRIMI
      COEF(ICOEF+K) = CONTRC(IEXOLE+K)
640   CONTINUE
      ICOEF = ICOEF + NPRIMI
630   CONTINUE
      IEXOLE = IEXOLE + NCONTC(J)
620   CONTINUE
625   CONTINUE
      NCONTR(IGRP) = NCONTI
      IEX = IEX + NPRIMI
610   CONTINUE
      IEXOLD = IEXOLD + NCONTC(I)
600   CONTINUE
      NGROUP = IGRP
C
C....    PRINT CENTRE INFORMATION
C
c      WRITE(6,106)
106   FORMAT(///,1X,'CENTRE',16X,'CO-ORDINATES',15X,'CHARGE',
     1       7X,'MASS',4X,'SPIN',4X,'Q',8X,'BASIS',//)
      NONNY = NOC - NOADC
      IF (JOBIPH .GT. 0) THEN
        WRITE(INTERF,4790) EXCLAM, BLANKS, BLANKS, NUCDAT
        WRITE(INTERF,4797) NOC
      ENDIF
4790  FORMAT(4A8)
4797  FORMAT(I8)
      DO 900 I = 1,NOC
c      WRITE(6,107) CENTR(I), (VLIST(I,J),J=1,10)
107   FORMAT(1X,A6,4X,3F12.8,4X,F4.1,2X,A2,2X,F8.4,2X,F4.1,2X,F6.3,2X,
     1       2A8)
      IF (JOBIPH .GT. 0)
     1  WRITE(INTERF,107) CENTR(I), (VLIST(I,J),J=1,10)
      IF (I .GT. NONNY) GOTO 900
c      WRITE(IOUTU) CENTR(I), (VLIST(I,J),J=1,10)
900   CONTINUE
C
C....    NOW WRITE OUT CENTER OF MASS AND ROTATIONAL CONSTANTS
C
c      WRITE(6,8840) TOTMSS, COM
c      WRITE(6,8841) RMOI
c      WRITE(6,8842) ((PATRAN(I,J),I = 1,3), J = 1,3)
8840  FORMAT(///,6X,' TOTAL MOLECULAR MASS IS ',F10.4,' DALTON',/,
     1           6X,' CENTER OF MASS IS ',3F12.8)
8841  FORMAT(/,' ROTATIONAL CONSTANTS (GHZ) ',3(F8.4,4X))
8842  FORMAT(/,' TRANSFORMATION TO PRINCIPAL AXES OF INERTIA ',/,
     1         (8X,F12.6,'x',F12.6,'y',F12.6,'z',/))
      IF (JOBIPH .GT. 0)
     1  WRITE(INTERF,4830) COM, RMOI
4830  FORMAT(3(4X,F12.8),/,3(4X,F8.4))
C
C....    PRINT BASIS INFORMATION IF DESIRED
C
      IF (IPFLG .EQ. 0) GOTO 750
c      WRITE(6,10000)
10000 FORMAT(1H1,20X,' GAUSSIAN BASIS INFORMATION',///,
     1       '  CENTER   TYPE               EXPONENT',
     2       '          CONTRACTION COEFFICIENTS',
     3       ' (WITH SEQUENCE NUMBERS)',/)
      IEX = 0
      ICOEF = 0
      ICNTR = 0
      DO 700 IGRP = 1, NGROUP
      NCONTS = NCONTR(IGRP)
      NPRIMS = NPRIM(IGRP)
c      ITIGRP = IBTAND(LABEL(IGRP),2**12 - 1)
      ITIGRP = IBTAND(LABEL(IGRP),2**10 - 1)
c      ICIGRP = IBTSHR(LABEL(IGRP),20)
      ICIGRP = IBTSHR(LABEL(IGRP),10)
c      WRITE(6,1144) (INDX(ICNTR+I),I = 1,NCONTS)
1144  FORMAT('0',40X,7(I9,3X),/,(43X,7(I9,3X)))
c      WRITE(6,1143) CENTR(ICIGRP), TYPE(ITIGRP)
1143  FORMAT(2X,A6,3X,A5)
      DO 710 I = 1,NPRIMS
c      IOVERP = 1H
c      IF (I .EQ. 1) IOVERP = 1H+
c      WRITE(6,1145) IOVERP, IEX+I, EXPS(IEX+I),
c     1              (COEF(ICOEF+NPRIMS*(J-1)+I),J = 1,NCONTS)
1145  FORMAT(A1,16X,I4,5X,F15.6,7F12.6,/,(43X,7F12.6))
710   CONTINUE
      IEX = IEX + NPRIMS
      ICOEF = ICOEF + NCONTS*NPRIMS
      ICNTR = ICNTR + NCONTS
700   CONTINUE
750   CONTINUE
C
C....    NOW RENORMALIZE CONTRACTIONS
C
      IEX = 0
      ICOEF = 0
      DO 650 IGRP = 1,NGROUP
      MMM = IBTAND(LABEL(IGRP),2**10 - 1)
      JE = NPRIM(IGRP)
      IE = NCONTR(IGRP)
      DO 660 I = 1,IE
      DTYPE=1.0
      IF(IMOL.NE.0) GO TO 70
      DTYPE=SQRT(ANORM(MMM))
70    CONTINUE
      GII = 0.0
      DO 14 J=1,JE
      EE=EXPS(IEX+J)*2.0
      IJ=MAP(MMM)
      GJJ = ANORM(MMM)*(SQRT(PI/EE)**3)*((0.5/EE)**(IJ-1))
      COEF(ICOEF+J) = COEF(ICOEF+J)/SQRT(GJJ)
      DO 14 K=1,J
      EE=EXPS(IEX+J)+EXPS(IEX+K)
      GI = ANORM(MMM)*(SQRT(PI/EE)**3)*((0.5/EE)**(IJ-1))
      GI=GI*COEF(ICOEF+J)*COEF(ICOEF+K)
      IF(J.EQ.K) GI = GI*0.5
14    GII = GII + GI
      GII = 1.0/SQRT(GII+GII)
      GII=GII*DTYPE
      DO 18 J =1,JE
      COEF(ICOEF+J) = COEF(ICOEF+J)*GII
18    CONTINUE
      ICOEF = ICOEF + JE
660   CONTINUE
      IEX = IEX + JE
650   CONTINUE
C
C....    PRINT RENORMALIZED BASIS INFORMATION IF DESIRED
C
      IF (IPFLG .LE. 1) GOTO 850
c      WRITE(6,11000)
11000 FORMAT(1H1,20X,' RENORMALIZED GAUSSIAN BASIS INFORMATION',///,
     1       '  CENTER   TYPE               EXPONENT',
     2       '          CONTRACTION COEFFICIENTS',
     3       ' (WITH SEQUENCE NUMBERS)',/)
      IEX = 0
      ICOEF = 0
      ICNTR = 0
      DO 800 IGRP = 1, NGROUP
      NCONTS = NCONTR(IGRP)
      NPRIMS = NPRIM(IGRP)
c      ITIGRP = IBTAND(LABEL(IGRP),2**12 - 1)
c      ICIGRP = IBTSHR(LABEL(IGRP),20)
      ITIGRP = IBTAND(LABEL(IGRP),2**10 - 1)
      ICIGRP = IBTSHR(LABEL(IGRP),10)
c      WRITE(6,1144) (INDX(ICNTR+I),I = 1,NCONTS)
c      WRITE(6,1143) CENTR(ICIGRP), TYPE(ITIGRP)
      DO 810 I = 1,NPRIMS
c      IOVERP = 1H
c      IF (I .EQ. 1) IOVERP = 1H+
c      WRITE(6,1145) IOVERP, IEX+I, EXPS(IEX+I),
C     1              (COEF(ICOEF+NPRIMS*(J-1)+I),J = 1,NCONTS)
810   CONTINUE
      IEX = IEX + NPRIMS
      ICOEF = ICOEF + NCONTS*NPRIMS
      ICNTR = ICNTR + NCONTS
800   CONTINUE
850   CONTINUE
111   FORMAT(1X,I5,2I10,8X,A6,2X,A6,3X,2(F13.7,2X))
53    CONTINUE
      REWIND INU
      CALL LOCATE(INU,'MULPOPDA')
      READ(INU)NCENTS, NCONR
c      WRITE(IOUTU) NCONR
      READ(INU) NCENTS,( CENTS(I),CHARGS(I),I=1,NCENTS)
c      WRITE(IOUTU) NCENTS,( CENTS(I),CHARGS(I),I=1,NCENTS)
      DO 50 I=1,NCONR
      READ(INU) IFK,ITRM(I),ITRN(I,1)
c      IF(IPFLG.GT.1)WRITE(6,102)I,ITRM(I),ITRN(I,1)
102   FORMAT(' POPULATION INFO ',3I5)
50    CONTINUE
c      WRITE(IOUTU)(ITRN(I,1),I=1,NCONR)
c      WRITE(IOUTU)(ITRM(I),I=1,NCONR)
C     READ IN TRANSFORMATION MATRIX
      REWIND INU
      CALL LOCATE(INU,'SYMTRANS')
      READ(INU) NSABF
      DO 52 I=1,NCONT
      ITRM(I)=0
52    CONTINUE
      DO 55 I=1,NSABF
      READ(INU) K,(ITRAN(L),CTRAN(L),L=1,K)
c      IF(IPFLG.NE.0)WRITE(6,88)I,K,(ITRAN(L),CTRAN(L),L=1,K)
88    FORMAT(' FUNCTION ',I4,'-',I4,(8(I3,F7.3)))
      DO 56 L=1,K
      IR=ITRAN(L)
      ITRM(IR)=ITRM(IR)+1
      IPUT=ITRM(IR)
      ITRN(IPUT,IR)=I
      CTRN(IPUT,IR)=CTRAN(L)
56    CONTINUE
55    CONTINUE
      RETURN
      END
