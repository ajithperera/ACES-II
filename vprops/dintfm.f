      SUBROUTINE DINTFM(OPA,NUC,VLIST,NOC,C,EXPS,CONT,LABEL,
     1 NPRIM,NCONTR,ITRM,ITRN,CTRN,NCONT,BUF,IBUF,AOINT,LBUF,
     2 INDX,NGROUP,VS,SCR,MXPRM2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND
      EXTERNAL OPA,NUC
      character*8 duml1,labeld(3),LABELMOM(6),LABELJA(10,25),junklab
C
CJDW 1/8/98. For f90 on crunch, make this array a character.
      CHARACTER*8 TITLE
      COMMON /PRTBUG/ IPRTFS
      COMMON /PRTDAT/ IDOSPH, MAXAQN, NSABF
      COMMON /NUCINDX/ ICT
      COMMON /PARM/ISL2,ISL3,ISL4,islx,THRESH,MTYP,ILNMCM,INT,INUI,
     & NCENTS,
     1 ID20,NSYM,NAO(12),NMO(12),junk,PTOT(3,10), INTERF, JOBIPH,
     2 RUNMOS
      DIMENSION VLIST(NOC,10),EXPS(1),CONT(1),LABEL(1),ITRM(1),
     1 ITRN(256,1),CTRN(256,1),NCONTR(1),NPRIM(1),INDX(1)
      DIMENSION IHOLD(10),IADH(10), VS(MXPRM2,10), CKE(10)
     1 ,IKE(3,10),JKE(3,10),BKE(10)
      DIMENSION BUF(LBUF,10),IBUF(IINTFP*LBUF,10),AOINT(1)
      DIMENSION VAL(10),TITLE(10,25)
      DIMENSION C(1), SCR(1),gbuf(600),igndx(600)
      COMMON /HIGHL/ LMNVAL(3,84), ANORM(84)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      DIMENSION LMN(200),JMN(200)
      COMMON /VARS/ FNU(20),DP(20),FN(15),FD(15),PI,DFTR(11)
      DIMENSION Q(3),B(3),E(3),D(3),AA(200),BB(200),V(10)
      DATA LABELD /'DIPOLE_X','DIPOLE_Y','DIPOLE_Z'/
C
      DATA (TITLE(I,1),I = 1,10)
     1/'   1/R  ',9*'        '/
C     DATA (TITLE(I,1),I = 1,10)
C    1/6H   1/R,9*6H      /
      DATA (LABELJA(I,1),I = 1,10)
     1/'ONEOVERR',9*'        '/   
C
      DATA (TITLE(I,2),I = 1,10)
     2/'    EX  ','    EY  ','    EX  ',7*'        '/
C     DATA (TITLE(I,2),I = 1,10)
C    2/6H    EX,6H    EY,6H    EX,7*6H      /
      DATA (LABELJA(I,2),I = 1,10)
     1/'EFIELD_X','EFIELD_Y','EFIELD_Z',7*'        '/   
C
      DATA (TITLE(I,3),I = 1,10)
     3/'   FXX  ','   FYY  ','   FZZ  ','   FXY  ','   FXZ  ',
     3 '   FYZ  ',4*'        '/
C     DATA (TITLE(I,3),I = 1,10)
C    3/6H   FXX,6H   FYY,6H   FZZ,6H   FXY,6H   FXZ,6H   FYZ,4*6H      /
      DATA (LABELJA(I,3),I = 1,10)
     1/'FGRDX000','FGRDY000','FGRDZ000',7*'        '/
C
      DATA (TITLE(I,4),I = 1,10)
     4/'     X  ','     Y  ','     Z  ',7*'        '/
C     DATA (TITLE(I,4),I = 1,10)
C    4/6H     X,6H     Y,6H     Z,7*6H      /
      DATA (LABELJA(I,4),I = 1,10)
     1/'DIPOLE_X','DIPOLE_Y','DIPOLE_Z',7*'        '/
C
      DATA (TITLE(I,5),I = 1,10)
     5/'   QXX  ','   QYY  ','   QZZ  ','   QXY  ','   QXZ  ',
     5 '   QYZ  ','  R**2  ','    XX  ','    YY  ','    ZZ  '/
C     DATA (TITLE(I,5),I = 1,10)
C    5/6H   QXX,6H   QYY,6H   QZZ,6H   QXY,6H   QXZ,6H   QYZ,6H  R**2,
C    5 6H    XX,6H    YY,6H    ZZ/
      DATA (LABELJA(I,5),I = 1,10)
     6/'QUAD_XX ','QUAD_YY ','QUAD_ZZ ','QUAD_XY ','QUAD_XZ ',
     & 'QUAD_YZ ','QUAD_R2 ','MYSTERY ','MYSTERY ','MYSTERY '/
C
      DATA (TITLE(I,6),I = 1,10)
     6/'  OXXX  ','  OYYY  ','  OZZZ  ','  OXXY  ','  OXXZ  ',
     6 '  OXYY  ','  OYYZ  ','  OXZZ  ','  OYZZ  ','  OXYZ  '/
C     DATA (TITLE(I,6),I = 1,10)
C    6/6H  OXXX,6H  OYYY,6H  OZZZ,6H  OXXY,6H  OXXZ,6H  OXYY,6H  OYYZ,
C    6 6H  OXZZ,6H  OYZZ,6H  OXYZ/
      DATA (LABELJA(I,6),I=1,10)
     &/'OCTUPXXX','OCTUPYYY','OCTUPZZZ','OCTUPXXY','OCTUPXXZ',
     & 'OCTUPXYY','OCTUPYYZ','OCTUPXZZ','OCTUPYZZ','OCTUPXYZ'/
C
      DATA (TITLE(I,7),I = 1,10)
     7/'    DX  ','    DY  ','    DZ  ',7*'        '/
C     DATA (TITLE(I,7),I = 1,10)
C    7/6H    DX,6H    DY,6H    DZ,7*6H      /
      DATA (LABELJA(I,7),I=1,10)
     &/'PLNDEN_X','PLNDEN_Y','PLNDEN_Z',7*'        '/
C
      DATA (TITLE(I,8),I = 1,10)
     8/'   DXY  ','   DXZ  ','   DYZ  ',7*'        '/
C     DATA (TITLE(I,8),I = 1,10)
C    8/6H   DXY,6H   DXZ,6H   DYZ,7*6H      /
      DATA (LABELJA(I,8),I=1,10)
     &/'LINDEN_Z','LINDEN_Y','PLNDEN_X',7*'        '/
C
      DATA (TITLE(I,9),I = 1,10)
     9/'   DEN  ',9*'        '/
      DATA (LABELJA(I,9),I=1,10)
     &/'RHOSQ000',9*'        '/
C
      DATA (TITLE(I,10),I = 1,10)
     O/'   OVL  ',9*'        '/
C     DATA (TITLE(I,10),I = 1,10)
C    O/6H   OVL,9*6H      /
C
      DATA (TITLE(I,11),I = 1,10)
     1/'CHIDXX  ','CHIDYY  ','CHIDZZ  ','CHIDXY  ','CHIDXZ  ',
     1 'CHIDYZ  ',4*'        '/
C     DATA (TITLE(I,11),I = 1,10)
C    1/6HCHIDXX,6HCHIDYY,6HCHIDZZ,6HCHIDXY,6HCHIDXZ,6HCHIDYZ,4*6H      /
C
      DATA (TITLE(I,12),I = 1,10)
     2/'   P4   ','DARWIN  ','TOTAL   ',7*'        '/
C     DATA (TITLE(I,12),I = 1,10)
C    2/6H   P4 ,6HDARWIN,6HTOTAL ,7*6H      /
      DATA (LABELJA(I,12),I = 1,10)
     2/'P**4    ','DARWIN  ','TOTAL   ',7*'        '/
C
      DATA (TITLE(I,13),I = 1,10)
     3/'   HFS  ',9*'        '/
C     DATA (TITLE(I,13),I = 1,10)
C    3/6H   HFS,9*6H      /
C
      DATA (TITLE(I,14),I = 1,10)
     4/'   2XX  ','   2YY  ','   2ZZ  ','   2XY  ','   2XZ  ',
     4 '   2YZ  ',4*'        '/
C     DATA (TITLE(I,14),I = 1,10)
C    4/6H   2XX,6H   2YY,6H   2ZZ,6H   2XY,6H   2XZ,6H   2YZ,4*6H      /
      DATA (LABELJA(I,14),I = 1,10)
     & /'2NDMO_XX','2NDMO_YY','2NDMO_ZZ',
     & '2NDMO_XY','2NDMO_XZ','2NDMO_YZ',4*'        '/
C
      DATA (TITLE(I,15),I = 1,10)
     5/'   XXX  ','   YYY  ','   ZZZ  ','   XXY  ','   XXZ  ',
     5 '   XYY  ','   YYZ  ','   XZZ  ','   YZZ  ','   XYZ  '/
C     DATA (TITLE(I,15),I = 1,10)
C    5/6H   XXX,6H   YYY,6H   ZZZ,6H   XXY,6H   XXZ,6H   XYY,6H   YYZ,
C    5 6H   XZZ,6H   YZZ,6H   XYZ/
C
      DATA (TITLE(I,16),I = 1,10)
     &/'VELOCI  ',9*'        '/
      DATA (LABELJA(I,16),I=1,10) / ' VELOCI ', 9*'        '/
C
C     DATA (TITLE(I,16),I = 1,10)
C    &/6HVELOCI,9*6H      /
C
      DATA ((TITLE(I,J),I = 1,10), J = 17,22) /60*'        '/
C     DATA ((TITLE(I,J),I = 1,10), J = 17,22) /60*6H      /
C
      DATA (TITLE(I,23),I = 1,10)
     3/' FX''X''  ',' FY''Y''  ',' FZ''Z''  ',' FX''Y''  ',' FX''Z''  ',
     3 ' FY''Z''  ',4*'        '/
C     DATA (TITLE(I,23),I = 1,10)
C    3/6H FX'X',6H FY'Y',6H FZ'Z',6H FX'Y',6H FX'Z',6H FY'Z',4*6H      /
C
      DATA (TITLE(I,24),I = 1,10)
     4/' QX''X''  ',' QY''Y''  ',' QZ''Z''  ',' QX''Y''  ',' QX''Z''  ',
     4 ' QY''Z''  ',4*'        '/
C     DATA (TITLE(I,24),I = 1,10)
C    4/6H QX'X',6H QY'Y',6H QZ'Z',6H QX'Y',6H QX'Z',6H QY'Z',4*6H      /
C
      DATA (TITLE(I,25),I = 1,10)
     3/'  X''X''  ','  Y''Y''  ','  Z''Z''  ','  X''Y''  ','  X''Z''  ',
     3 '  Y''Z''  ',4*'        '/
C     DATA (TITLE(I,25),I = 1,10)
C    3/6H  X'X',6H  Y'Y',6H  Z'Z',6H  X'Y',6H  X'Z',6H  Y'Z',4*6H      /
C
      IBTAND(I,J) =IAND(I,J)
      IBTOR(I,J)  =IOR(I,J)
      IBTXOR(I,J) =IEOR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTNOT(I)   = NOT(I)
      UNITY = 1.0
c      write(6,*)' ict is ',ict
      LXQ=0
      RELFAC=137.03604**2*8.0
      RELFAC=-1.0/RELFAC
      RF2=PI/(137.03604**2*2.)
      THR2=THRESH/100000.
      CALL NUC(NOC,VLIST,C,NOC,NT)
      IF(MTYP.EQ.12)NT=3
      CALL SETRHF
      IWPB=(LBUF-2)/2
c      IWPB=(LBUF-2)
      lenbuf=iintfp*lbuf
      DO 100 I=1,NT
      IHOLD(I)=0
      IBUF(1,I)=0
CRAY 1
      IBUF(2,I)=-1
100   CONTINUE
      IEXH=0
      ICONT = 0
      ICOFF = 0
      DO 50 I=1,NGROUP
CRAY 1
c      ICNT=IBTAND(IBTSHR(LABEL(I),20),2**12 -1)
      ICNT=AND(IBTSHR(LABEL(I),10),2**10-1)
c      ITYP=IBTAND(LABEL(I),2**15 - 1)
      ITYP=IBTAND(LABEL(I),2**10-1)
      DO 410 II=1,3
      Q(II)=VLIST(ICNT,II)
410   CONTINUE
      IE=NPRIM(I)
      ICE = NCONTR(I)
      JEXH=0
      JCONT = 0
      JCOFF = 0
      DO 49 J=1,I
CRAY 1
      JCNT=IBTAND(IBTSHR(LABEL(J),10),2**10 - 1)
c      JTYP=IBTAND(LABEL(J),2**15 - 1)
      JTYP=IBTAND(LABEL(J),2**10-1)
      DO 411 JJ=1,3
      B(JJ)=VLIST(JCNT,JJ)
411   CONTINUE
      JE=NPRIM(J)
      JCE = NCONTR(J)
      DO 2  M=1,NT
      DO 28 MM = 1,MXPRM2
      VS(MM,M)=0.0
28    CONTINUE
2     CONTINUE
      IEX=IEXH
      DO 40 II=1,IE
      IEX=IEX+1
      JEX=JEXH
      DO 40 JJ=1,JE
      JEX=JEX+1
      IJIND = IE*(JJ-1) + II
      IF(MTYP.NE.12)GO TO 300
C
C RELATIVISTIC INTEGRALS
C
C     P4  TWO AND TWO
C
      CALL CIVPT(Q,EXPS(IEX),B,EXPS(JEX),E,GAMA,EFACT)
      LOOPI=4
      IL1=LMNVAL(1,ITYP)
      IM1=LMNVAL(2,ITYP)
      IN1=LMNVAL(3,ITYP)
      CKE(1)=-2.*EXPS(IEX)*(2*IL1+2*IM1+2*IN1+3)
      IKE(1,1)=IL1+1
      IKE(2,1)=IM1+1
      IKE(3,1)=IN1+1
      CKE(2)=4.*EXPS(IEX)**2
      IKE(1,2)=IL1+3
      IKE(2,2)=IM1+1
      IKE(3,2)=IN1+1
      CKE(3)=4.*EXPS(IEX)**2
      IKE(1,3)=IL1+1
      IKE(2,3)=IM1+3
      IKE(3,3)=IN1+1
      CKE(4)=4.*EXPS(IEX)**2
      IKE(1,4)=IL1+1
      IKE(2,4)=IM1+1
      IKE(3,4)=IN1+3
      IF(IL1.LT.2)GO TO 341
      LOOPI=LOOPI+1
      CKE(5)=2.*DFLOAT(IL1*(IL1-1)/2)
      IKE(1,5)=IL1-1
      IKE(2,5)=IM1+1
      IKE(3,5)=IN1+1
341   CONTINUE
      IF(IM1.LT.2)GO TO 342
      LOOPI=LOOPI+1
      CKE(LOOPI)=2.*DFLOAT(IM1*(IM1-1)/2)
      IKE(1,LOOPI)=IL1+1
      IKE(2,LOOPI)=IM1-1
      IKE(3,LOOPI)=IN1+1
342   CONTINUE
      IF(IN1.LT.2)GO TO 343
      LOOPI=LOOPI+1
      CKE(LOOPI)=2.*DFLOAT(IN1*(IN1-1)/2)
      IKE(1,LOOPI)=IL1+1
      IKE(2,LOOPI)=IM1+1
      IKE(3,LOOPI)=IN1-1
343   CONTINUE
      LOOPJ=4
      JL1=LMNVAL(1,JTYP)
      JM1=LMNVAL(2,JTYP)
      JN1=LMNVAL(3,JTYP)
      BKE(1)=-2.*EXPS(JEX)*(2*JL1+2*JM1+2*JN1+3)
      JKE(1,1)=JL1+1
      JKE(2,1)=JM1+1
      JKE(3,1)=JN1+1
      BKE(2)=4.*EXPS(JEX)**2
      JKE(1,2)=JL1+3
      JKE(2,2)=JM1+1
      JKE(3,2)=JN1+1
      BKE(3)=4.*EXPS(JEX)**2
      JKE(1,3)=JL1+1
      JKE(2,3)=JM1+3
      JKE(3,3)=JN1+1
      BKE(4)=4.*EXPS(JEX)**2
      JKE(1,4)=JL1+1
      JKE(2,4)=JM1+1
      JKE(3,4)=JN1+3
      IF(JL1.LT.2)GO TO 351
      LOOPJ=LOOPJ+1
      BKE(5)=2.*DFLOAT(JL1*(JL1-1)/2)
      JKE(1,5)=JL1-1
      JKE(2,5)=JM1+1
      JKE(3,5)=JN1+1
351   CONTINUE
      IF(JM1.LT.2)GO TO 352
      LOOPJ=LOOPJ+1
      BKE(LOOPJ)=2.*DFLOAT(JM1*(JM1-1)/2)
      JKE(1,LOOPJ)=JL1+1
      JKE(2,LOOPJ)=JM1-1
      JKE(3,LOOPJ)=JN1+1
352   CONTINUE
      IF(JN1.LT.2)GO TO 353
      LOOPJ=LOOPJ+1
      BKE(LOOPJ)=2.*dFLOAT(JN1*(JN1-1)/2)
      JKE(1,LOOPJ)=JL1+1
      JKE(2,LOOPJ)=JM1+1
      JKE(3,LOOPJ)=JN1-1
353   CONTINUE
C
C     DARWIN TERMS
C
      XX=0.0
      DO 420 IC=1,NOC
      XC1=VLIST(IC,1)-Q(1)
      YC1=VLIST(IC,2)-Q(2)
      ZC1=VLIST(IC,3)-Q(3)
      XC2=VLIST(IC,1)-B(1)
      YC2=VLIST(IC,2)-B(2)
      ZC2=VLIST(IC,3)-B(3)
      CHRG=VLIST(IC,4)
      EXPI=-EXPS(IEX)*(XC1*XC1+YC1*YC1+ZC1*ZC1)
      EXPJ=-EXPS(JEX)*(XC2*XC2+YC2*YC2+ZC2*ZC2)
      EXPI=EXP(EXPI)
      EXPJ=EXP(EXPJ)
      ANG=1.0
      IF(IL1.NE.0)ANG=XC1**IL1
      IF(IM1.NE.0)ANG=ANG*YC1**IM1
      IF(IN1.NE.0)ANG=ANG*ZC1**IN1
      IF(JL1.NE.0)ANG=ANG*XC2**JL1
      IF(JM1.NE.0)ANG=ANG*YC2**JM1
      IF(JN1.NE.0)ANG=ANG*ZC2**JN1
      ANG=ANG*EXPI*EXPJ
      XX=XX+ANG*CHRG
420   CONTINUE
      VS(IJIND,2)=VS(IJIND,2)+XX*RF2
      IL1=IL1+1
      IM1=IM1+1
      IN1=IN1+1
      JL1=JL1+1
      JM1=JM1+1
      JN1=JN1+1
      CALL OVPP(OX,Q(1),B(1),E(1),JL1,IL1,GAMA)
      CALL OVPP(OY,Q(2),B(2),E(2),JM1,IM1,GAMA)
      CALL OVPP(OZ,Q(3),B(3),E(3),JN1,IN1,GAMA)
      VALX=0.0
      DO 346 LL=1,LOOPJ
      DO 347 KK=1,LOOPI
      IF(JKE(1,LL).EQ.JL1.AND.IKE(1,KK).EQ.IL1)GO TO 500
      CALL OVPP(X,Q(1),B(1),E(1),JKE(1,LL),IKE(1,KK),GAMA)
      GO TO 501
500   CONTINUE
      X=OX
501   CONTINUE
      IF(X.EQ.0.0)GO TO 347
      IF(JKE(2,LL).EQ.JM1.AND.IKE(2,KK).EQ.IM1)GO TO 510
      CALL OVPP(Y,Q(2),B(2),E(2),JKE(2,LL),IKE(2,KK),GAMA)
      GO TO 511
510   CONTINUE
      Y=OY
511   CONTINUE
      IF(Y.EQ.0.0)GO TO 347
      IF(JKE(3,LL).EQ.JN1.AND.IKE(3,KK).EQ.IN1)GO TO 520
      CALL OVPP(Z,Q(3),B(3),E(3),JKE(3,LL),IKE(3,KK),GAMA)
      GO TO 521
520   CONTINUE
      Z=OZ
521   CONTINUE
      OLP=X*Y*Z
      VALX=CKE(KK)*BKE(LL)*OLP+VALX
347   CONTINUE
346   CONTINUE
      VS(IJIND,1)=(SQRT(PI/GAMA))**3*EFACT*VALX*RELFAC+VS(IJIND,1)
      VS(IJIND,3)=VS(IJIND,1)+VS(IJIND,2)
      GO TO 40
300   CONTINUE
      IF(MTYP.NE.16)GO TO 301
C
C     VELOCITY
C
      CALL CIVPT(Q,EXPS(IEX),B,EXPS(JEX),E,GAMA,EFACT)
      LOOPI=3
      IL1=LMNVAL(1,ITYP)
      IM1=LMNVAL(2,ITYP)
      IN1=LMNVAL(3,ITYP)
      CKE(1)=-2.*EXPS(IEX)
      IKE(1,1)=IL1+2
      IKE(2,1)=IM1+1
      IKE(3,1)=IN1+1
      CKE(2)=-2.*EXPS(IEX)
      IKE(1,2)=IL1+1
      IKE(2,2)=IM1+2
      IKE(3,2)=IN1+1
      CKE(3)=-2.*EXPS(IEX)
      IKE(1,3)=IL1+1
      IKE(2,3)=IM1+1
      IKE(3,3)=IN1+2
      IF(IL1.LT.1)GO TO 641
      LOOPI=LOOPI+1
      CKE(LOOPI)=DFLOAT(IL1)
      IKE(1,LOOPI)=IL1
      IKE(2,LOOPI)=IM1+1
      IKE(3,LOOPI)=IN1+1
641   CONTINUE
      IF(IM1.LT.1)GO TO 642
      LOOPI=LOOPI+1
      CKE(LOOPI)=DFLOAT(IM1)
      IKE(1,LOOPI)=IL1+1
      IKE(2,LOOPI)=IM1
      IKE(3,LOOPI)=IN1+1
642   CONTINUE
      IF(IN1.LT.1)GO TO 643
      LOOPI=LOOPI+1
      CKE(LOOPI)=DFLOAT(IN1)
      IKE(1,LOOPI)=IL1+1
      IKE(2,LOOPI)=IM1+1
      IKE(3,LOOPI)=IN1
643   CONTINUE
      JKE(1,1)=LMNVAL(1,JTYP)+1
      JKE(2,1)=LMNVAL(2,JTYP)+1
      JKE(3,1)=LMNVAL(3,JTYP)+1
      VALX=0.0
      DO 647 KK=1,LOOPI
      CALL OVPP(X,Q(1),B(1),E(1),JKE(1,1),IKE(1,KK),GAMA)
      CALL OVPP(Y,Q(2),B(2),E(2),JKE(2,1),IKE(2,KK),GAMA)
      CALL OVPP(Z,Q(3),B(3),E(3),JKE(3,1),IKE(3,KK),GAMA)
      OLP=X*Y*Z
      VALX=CKE(KK)*OLP+VALX
647   CONTINUE
      VS(IJIND,1)=(SQRT(PI/GAMA))**3*EFACT*VALX+VS(IJIND,1)
      GO TO 40
301   CONTINUE
C
C ALL OTHER PROPERTIES - RELATIVISTIC PROPS ARE ABOVE
C
      CALL CIVPT(Q,EXPS(IEX),B,EXPS(JEX),E,GAMA,EFACT)
      CALL RHFTCE(AA,Q,E,ITYP,ITM,UNITY,LMN)
      CALL RHFTCE(BB,B,E,JTYP,JTM,UNITY,JMN)
      DO 31 K=1,NT
      VAL(K)=0.
31    CONTINUE
      DO 25 K=1,3
      D(K)=E(K)-C(K)
25    CONTINUE
      DO 210 IIQ=1,ITM
      IF(AA(IIQ).EQ.0.0)GO TO 210
CRAY 2
c      IL=IBTAND(IBTSHR(LMN(IIQ),40),2**9 - 1)
c      IM=IBTAND(IBTSHR(LMN(IIQ),20),2**9 - 1)
c      IN=IBTAND(LMN(IIQ),2**9 - 1)
      IL=IBTAND(IBTSHR(LMN(IIQ),20),2**10 - 1)
      IM=IBTAND(IBTSHR(LMN(IIQ),10),2**10 - 1)
      IN=IBTAND(LMN(IIQ),2**10 - 1)
      DO 211 JJQ=1,JTM
      DD=AA(IIQ)*BB(JJQ)
      IF(ABS(DD).LT.THR2)GO TO 211
CRAY 2
c      JL=IBTAND(IBTSHR(JMN(JJQ),40),2**9 - 1) + IL
c      JM=IBTAND(IBTSHR(JMN(JJQ),20),2**9 - 1) + IM
c      JN=IBTAND(JMN(JJQ),2**9 - 1) + IN
      JL=IBTAND(IBTSHR(JMN(JJQ),20),2**10 - 1) + IL
      JM=IBTAND(IBTSHR(JMN(JJQ),10),2**10 - 1) + IM
      JN=IBTAND(JMN(JJQ),2**10 - 1) + IN
      CALL OPA(JL,JM,JN,GAMA,V,NT,D)
      DO 51 L=1,NT
      VAL(L)=DD*V(L)+VAL(L)
51    CONTINUE
211   CONTINUE
210   CONTINUE
C&&&&&&&DEBUG&&&&
      if (mtyp .eq. 2) then
      write(6,*) 
      write(6,*) 'Center1  ', i ,'     ',II
      write(6,*) 'Center2  ', j ,'     ',JJ
      Write (6, 999) (EFACT*VAL(IIII), IIII = 1 ,3)
 999  format(3X, 6(5X,F12.5)) 
      endif
C&&&&&&&&&&&&&&&&
      DO 30 M=1,NT
30    VS(IJIND,M) = EFACT*VAL(M) + VS(IJIND,M)
40    CONTINUE

      DO 101 M=1,NT
C
C....    PERFORM GENERAL CONTRACTION, SPECIAL CODE FOR 1X1 CASE
C
      IF (IE*JE .EQ. 1) THEN
        VS(1,M) = VS(1,M)*CONT(ICONT+1)*CONT(JCONT+1)
      ELSE
        CALL MXM(VS(1,M),IE,CONT(JCONT+1),JE,SCR,JCE)
        CALL MXMT(CONT(ICONT+1),ICE,SCR,IE,VS(1,M),JCE)
      ENDIF
      DO 111 II = 1, ICE
      IIND = INDX(ICOFF+II)
      JTOP = JCE
      IF (I.EQ.J) JTOP = II
      DO 112 JJ = 1, JTOP
      JJND = INDX(JCOFF+JJ)
      VSINT = VS(ICE*(JJ-1)+II,M)
      IF(ABS(VSINT).LT.THRESH)GO TO 112
      IHOLD(M)=IHOLD(M)+1
      IPUT=IHOLD(M)+2
      IF(JJND.GT.IIND.AND.MTYP.EQ.16)VSINT=-VSINT
      BUF(IPUT,M)=VSINT
C&&&&&&&DEBUG
C      Write (luout, 999) BUF(IPUT,M)
C 999  format(5X,F12.5, 5X,F12.5,5X,F12.5) 
C&&&&&&&DEBUG
      ioff=iput+iwpb
      IBUF(ioff*iintfp,M) =
     1    IBTOR(IBTSHL(MAX0(IIND,JJND),16),MIN0(IIND,JJND))
      
c      IBUF(IPUT+IWPB,M) =
c     1    IBTOR(IBTSHL(MAX0(IIND,JJND),15),MIN0(IIND,JJND))
c     1    IBTOR(IBTSHL(MAX0(IIND,JJND),16),MIN0(IIND,JJND))
      IF(IHOLD(M).LT.IWPB)GO TO 112
      IBUF(1,M)=IWPB
      CALL OUTAPD(IBUF(1,M),LenBUF,ID20,LXQ)
      IBUF(2,M)=LXQ
      IHOLD(M)=0
112   CONTINUE
111   CONTINUE
101   CONTINUE
      JEXH=JEXH+JE
      JCONT = JCONT + JE*JCE
      JCOFF = JCOFF + JCE
C&&&&&&DEBUG
C      Write(luout, 999) (BUF(j, m), m = 1, 3)
C&&&&&&&DEBUG
49    CONTINUE
      IEXH=IEXH+IE
      ICONT = ICONT + IE*ICE
      ICOFF = ICOFF + ICE
50    CONTINUE
C     CLOSE BINS
      DO 102 I=1,NT
      IF(IHOLD(I).EQ.0)GO TO 103
      IBUF(1,I)=IHOLD(I)
      CALL OUTAPD(IBUF(1,I),LenBUF,ID20,LXQ)
      IADH(I)=LXQ
      GO TO 102
103   CONTINUE
      IADH(I)=IBUF(2,I)
102   CONTINUE
C     NOW TRANSFORM TO SYMMETRY INTS
      NU2=NSABF*(NSABF+1)/2
      NU=NU2+3
c      WRITE(IOUTU) NT,NU
      DO 600 NX=1,NT
      DO 70 I=1,NU2
      AOINT(I)=0.0
70    CONTINUE
      LL=IADH(NX)
75    CONTINUE
      IF(LL.EQ.-1)GO TO 76
      CALL INTAPD(IBUF,LenBUF,ID20,LL)
      LOOP=IBUF(1,1)
      LL=IBUF(2,1)
      DO 74 IQ=1,LOOP
      X=BUF(IQ+2,1)
      ioff=iq+2+iwpb
      I=IBTAND(IBTSHR(IBUF(ioff*iintfp,1),16),2**16 - 1)
      J=IBTAND(IBUF(ioff*iintfp,1),2**16 - 1)
      
c      I=IBTAND(IBTSHR(IBUF(IQ+2+IWPB,1),16),2**16 - 1)
c      J=IBTAND(IBUF(IQ+2+IWPB,1),2**16 - 1)
      IE=ITRM(I)
      DO 702 II=1,IE
      CC=CTRN(II,I)
      IS=ITRN(II,I)
      JE=ITRM(J)
      DO 704 JJ=1,JE
      JS=ITRN(JJ,J)
      IF(IS.LT.JS)GO TO 704
      XP=CC*CTRN(JJ,J)*X
      IJF=IS*(IS-1)/2+JS
      AOINT(IJF)=XP+AOINT(IJF)
704   CONTINUE
702   CONTINUE
      IF(I.EQ.J)GO TO 74
      DO 722 II=1,JE
      CC=CTRN(II,J)
      IS=ITRN(II,J)
      DO 724 JJ=1,IE
      JS=ITRN(JJ,I)
      IF(IS.LT.JS)GO TO 724
      XP=CC*CTRN(JJ,I)*X
      IF(MTYP.EQ.16)XP=-XP
      IJF=IS*(IS-1)/2+JS
      AOINT(IJF)=XP+AOINT(IJF)
724   CONTINUE
722   CONTINUE
74    CONTINUE
      GO TO 75
76    CONTINUE
c      WRITE(IOUTU) (PTOT(I,NX),I=1,3),(AOINT(J),J=1,NU2)
C
C CRAPS "INTERFACE"
C
      ndprop=35
      lcbuf=600
      duml1='********'
C SG 12/6/96
      DUMLAB = 0.0D0
C
      write(ndprop)duml1,dumlab,dumlab,title(nx,mtyp),ptot(1,nx)
      icount=0
C
CSSSS         If (mtyp .eq. 2 .or. mtyp .eq. 24 .or. mtyp .eq.14 .or.
         If (mtyp .eq. 1) Then 
         Write(6,*) "Triangularly Packed AO integrals"
         Write(6,"(5F12.8)") (aoint(i), i=1, nu2)
         Endif

      if(mtyp.eq.4)then
       call putrec(20,'JOBARC',LABELD(NX),NU2*IINTFP,AOINT)
      elseif(mtyp.eq.3.or.mtyp.eq.9)then
       if(ict.le.9)then
        junklab=labelja(nx,mtyp)
        WRITE(JUNKLAB(8:8),'(I1)')ICT
       elseif(ict.gt.9.and.ict.lt.100)then
        WRITE(JUNKLAB(7:8),'(I2)')ICT
       else
        WRITE(JUNKLAB(6:8),'(I3)')ICT
       endif
       call putrec(20,'JOBARC',JUNKLAB,NU2*IINTFP,AOINT)
      else
       call putrec(20,'JOBARC',LABELJA(NX,MTYP),NU2*IINTFP,AOINT)
      endif
      do 750 ijt=1,nu2
       if(icount.lt.lcbuf)goto 751
       limit=lcbuf
       write(ndprop) gbuf,igndx,limit
       icount=0
751    icount=icount+1
       igndx(icount)=ijt
       gbuf(icount)=aoint(ijt)
750   continue
      if(icount.gt.0)then
       limit=icount
       write(ndprop)gbuf,igndx,limit
      endif
      limit=-1
      write(ndprop)gbuf,igndx,limit
C
C       
C  END OF "INTERFACE" 
C      
      IF(ISL2.EQ.0)GO TO 600
c      CALL PRINTM(AOINT,NCONT)
600   CONTINUE
c      close(26)
      RETURN
      END
