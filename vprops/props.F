C
      SUBROUTINE PROPS(S,MAXCOR,IS,IMAXCOR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C Modified 08/93 Ajith
      PARAMETER (NUMPROP = 26)
      PARAMETER (NGAUSS = 26)
      PARAMETER (VALUEM = 0.8D+00)      
C
C-----------------------------------------------------------------------
C Modifications for f90. 
C
C     IPRPCT is initialized (should be done anyway).
C
C     See routines INCENT, GETCNT.
C
C     Make ILNMCM an integer throughout this main program. It is not
C     clear if ILNMCM has any use, but various logic depends on its
C     value, so we'll try to do something.
C
C     Initialize it to zero at beginning.
C     Usually it is set to values in IS (S(IB)), which are set in
C     INCENT. Of course, ILNMCM is then an integer representation of
C     a character*4 quantity. Also, the equivalencing between S and IS
C     probably does not work on 32-bit integer machines. As for explicit
C     character stuff, this is what we are doing.
C
C     Previous    New
C     6H          999
C     6HC-O-M  999999
C
C     The original hollerith code for ILNMCM is commented.
C-----------------------------------------------------------------------
C
      Logical ISTHER,IOPN
C
C NATMMX is the maximum number of atoms allowed.  To increase the size
C of the molecule allowed, increase NATMMX.  NVARMX is the proper
C size of VARM and VARIC for a given NATMMX.
C SG 10/94
C
      PARAMETER (NATMMX = 100,NVARMX = (9+9*NATMMX+1))
      INTEGER VARM(NVARMX),VARIC(NVARMX)
      CHARACTER*8 STARS,BLANKS,PROPSS,EXCLAM,HEAD,EODATA
      CHARACTER*8 ATOMSYM(NATMMX)
      CHARACTER*20 PRPTYP(NUMPROP)
      COMMON /PRTBUG/ IPRTFS
      COMMON /PRTDAT/ IDOSPH, MAXAQN, NSABF
      COMMON /COMINF/ TOTMSS, COM(3), RMOI(3), PATRAN(3,3)
      COMMON /NUCINDX/ ICT
 
#include "icdacc.com"

      COMMON/FILES/LUOUT, MOINTS
      COMMON/FLAGS/IFLAGS(100)
      COMMON/LABS/ ILAB(12),ILBL(12),A(20,3),TITLE(10,25),ICOD(12)
      COMMON /PARM/ISL2,ISL3,ISL4,islx,THRESH,M,ILNMCM,INT,INUI,
     & NCENTS,
     1 ID20,NSYM,NAO(12),NMO(12),junk,PTOT(3,10), INTERF, JOBIPH,
     2 RUNMOS
      COMMON /VARS/ FNU(20),DP(20),FN(15),FD(15),PI,DFTR(11)
      EXTERNAL OPAA0,OPAA1,OPAA2,OPAB1,OPAB2,OPAB3,OPAC1,OPAC2, OPAC3,
     1         OPAD1,OPAB4, OPAB5, OPAB6
      EXTERNAL POTNUC,FLDNUC,GRDNUC,DIPNUC,QUDNUC,OCTNUC,NONUC,MAGNUC,
     1         SECNUC, THINUC
      DIMENSION S(MAXCOR)
      DIMENSION IS(IMAXCOR),IPART(3)
      COMMON /CWBP/THRSMP,ISEPOP,NPART(25),IEXPOP,junque,PRCONV(20),
     &             UCONV(20)
      DIMENSION C(3), Centr1(3), Centr2(3), LAB(10)
      DIMENSION IDOCOM(22), IDOPRP(2,20)
      DIMENSION IRDMSS(25), XMASS(25), XSPIN(25), XNAME(25), XQ(25)
      DIMENSION ISETIN(10)
C Added 08/93 Ajith
      DIMENSION ROOT(NGAUSS), WEIGHT(NGAUSS)
C
      EQUIVALENCE (IPMOI,ISL3), (IPAOI,ISL2),(IDM,ISL4)
      NAMELIST /INPUTP/ IPMOI, IPAOI, INT, IPROP, ISETIN,INU,
     1 THRESH,IALL,IDM,IBUF,IPORB, IPBF, THRSMP,ISEPOP,IOLDV,
     2 NPART,IPART,IEXPOP, PRCONV, IDOCOM, IAUTO, IRDMSS, XMASS,
     3 JOBIPH, IACC
      DATA  IALL/0/, IBUF/512/
      DATA ISETIN/10*1/, IPORB/0/,IPBF/0/,IPART/3*0/
      DATA IAUTO/0/
      DATA INUCMX/25/
      DATA IRDMSS/25*0/
      DATA IPASS/0/
      DATA IACC/0/
C Added 08/93 Ajith
      DATA PRPTYP /'Potential           ','Electric Field      ',
     &             'Field gradient      ','Dipole moment       ',
     &             'Quadrupole moment   ','Octopole moment     ',
     &             'Planar density      ','Linear density      ',
     &             'Electron density    ','Overlap             ',
     &             'Diamagnetic suscep. ','Relativisitic corr. ',
     &             'Hyperfine splitting ','Second moment       ',
     &             'Third moment        ','Velocity            ',
     &             'Orbital diamagnetic ','Orbital paramagnetic',
     &             'Spin-dipole         ','Diamagnetic shield. ',
     &             'Second moment       ','Dipole velocity     ',
     &             'Grad. potential     ','Angular Momentum    ',
     &             'Blank               ','Mag. Quadrupole     '/
C
      DATA INU /16/
      DATA IDOCOM /22*-1/
      DATA IDOPRP / 10,0, 12,0, 4,0, 5,0, 6,0, 3,-1, 11,0, 0,0, 24*0/
      DATA STARS /'********'/, BLANKS/'        '/,PROPSS/'PROPS   '/
      DATA EXCLAM /'!!!!!!!!'/, HEAD/'HEADER  '/, EODATA/'ENDFILE '/
      CALL ZERO(S,MAXCOR)
C
C INITIALIZE SOME VARIABLES.  MOVED FROM DATA STATEMENTS WHICH DID
C  NOT FOLLOW ANSI F77
C
      MAXAQN=0
      OneM = -1.00D+00
      One  = 1.00D+00
      THRSMP=0.005
      ISEPOP=1
      IEXPOP=0
      THRESH=1.E-12
      IDM=0
      IPMOI=1
      IPAOI=0
      INT=0
      INUI=5
      PI=DACOS(-1.0D0)
      DFTR(1)=1.0
      DO 104 I=2,11
       DFTR(I)=DFTR(I-1)*(1.0+2.0*(DFLOAT(I)-2.0))
104   CONTINUE
C
C     IPRPCT needs to be initialized.
C
      IPRPCT=0
C
      ID20=60
      INTERF=74
      JOBIPH=-1
      WRITE(6,151)
      WRITE(6,*)
      WRITE(6,152)
      WRITE(6,153)
      WRITE(6,154)
      WRITE(6,152)
151   FORMAT(T4,'Property integrals will be calculated.')
152   FORMAT(4X,56('-'))
153   FORMAT(T33,'Property evaluated at ')
154   FORMAT(T5,'Property type',T32,'X',T44,'Y',T56,'Z')
C
      OPEN(UNIT=INU,FILE='IIII',FORM='UNFORMATTED',STATUS='OLD')
      OPEN(UNIT=35,FILE='VPOUT',FORM='UNFORMATTED',STATUS='UNKNOWN')
      MXPRM = 25
      MXPRM2 = MXPRM*MXPRM
C
C....    SET FACTORIAL ARRAYS.  THE DIMENSION 15 WILL HANDLE UP TO I TYP
C....    FUNCTIONS.
C
      FN(1) = 1.
      FD(1) = 1.
      DO 7700 I = 2,15
         FN(I) = (I - 1)*FN(I-1)
         FD(I) = 1./FN(I)
7700  CONTINUE
CJFS
      INQUIRE(FILE='VPINP',EXIST=ISTHER)
      IF (ISTHER) THEN
         OPEN (UNIT=INUI, FILE='VPINP', FORM='FORMATTED', STATUS='OLD')
      ENDIF
CJFS
      IF(ISTHER)THEN
c       READ(5,101) (ILBL(I), I=1,12)
      ELSE
       IPBF=0
       IPAOI=0
       INT=1
      ENDIF
c      REWIND IOUTU
      IF(ISTHER)READ(INUI,INPUTP)
c      IF(INT.LT.0)IPROP(1)=0
C
C....    SET DEFAULTS FOR DOING PROPERTIES IN CENTER OF MASS AND FOR
C....    OUTPUT INTERFACE
C
      IF (INT .GE. 0) THEN
        IF (IDOCOM(5) .EQ. -1) IDOCOM(5) = 1
        IF (IDOCOM(6) .EQ. -1) IDOCOM(6) = 1
        IF (IDOCOM(11) .EQ. -1) IDOCOM(11) = 1
      ENDIF
c      WRITE(6,INPUTP)
c YAU - I am not sure if this should be 2*idaccm or iintfp*idaccm.
      OPEN(ID20,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=IBUF*2*idaccm)
4790  FORMAT(4A8)
4795  FORMAT(12A6)
4799  FORMAT(9A8)
C
C....    AFTER THIS JOBIPH MUST GE .GE. 0
C
      
      CALL GETSPH(INU)
      CALL SETLMN
      IF(INT.EQ.0) GOTO 200
      REWIND INU
      CALL LOCATE(INU,'POLYAINP')
      READ(INU)NON,NOADC,IKK,IMOL
      NOC=NON+NOADC
      IA=1
      IB=IA+NOC*10
      IC=IB+NOC
cSSS CRAY1
      IF(IC.GT.MAXCOR)CALL INSMEM('VPROPS',IC,MAXCOR)
CDC 3
C      IPP=0
C      NWORD=LOCF(S(10))+IC
C      CALL GETCOR(NWORD,IPP)
C
C....    PRESET MASSES BEFORE READING NUCLEAR DATA
C
      DO 77840 I = 1,NOC
      S(IA+5*NOC-1+I) = -1.
77840 CONTINUE
      DO 77850 I = 1,INUCMX
      IF (IRDMSS(I) .LE. 0) GOTO 77860
      INDI = IA + 5*NOC - 1 + IRDMSS(I)
      S(INDI) = XMASS(I)
77850 CONTINUE
77860 CONTINUE
#ifdef _UNICOS
      CALL INCENT(S(IA),ATOMSYM,NOC,NON,NOADC,INU,S(IB),IOUTU)
#else
      CALL INCENT(S(IA),NOC,NON,NOADC,INU,S(IB),IOUTU)
#endif /* _UNICOS */
      READ(INU)NUNC,NCONT
      ID=IC+NUNC
      IE=ID+NUNC
      IF=IE+NCONT
      IG=IF+NCONT
      IH=IG+NCONT
      II=IH+NCONT*256
      IJ=II+NCONT*256
      IK=IJ+NOC
      IK1=IK+NOC
      IK2=IK1+NCONT
      IK3=IK2+NCONT
      IK4=IK3+NCONT
      IK5=IK4+NCONT
      IK6=IK5+10*MXPRM2
      IKS1=IK6+MXPRM2
      IKS2=IKS1+NCONT
      IKS3=IKS2+NUNC
      IL2=IKS3+NUNC
C     CORE FOR INT SORT
      IKN=IKS1+10*IBUF
      IN=IKS1+IBUF
      IM=MAX0(IN+NCONT*(NCONT+1)/2,IKN,IL2)
cSSS CRAY 1
      IF(IM.GT.MAXCOR)CALL INSMEM('VPROPS',IM,MAXCOR)
CDC 3
C      IPP=0
C      NWORD=LOCF(S(10))+IM
C      CALL GETCOR(NWORD,IPP)
      CALL INPUT(S(IC),S(ID),S(IE),S(IF),S(IG),S(IH),S(II),NOC,
     1 NCONT,INU,IOUTU,S(IA),S(IB),IS(IJ),S(IK),IMOL,IPBF,
     2 S(IK1), S(IK2), S(IK3), S(IK4), S(IKS1), S(IKS2), S(IKS3),
     3 NGROUP, NOADC, JOBIPH, INTERF)
C
C Write atomic ordering to the 'JOBARC' file
      CALL DUMPJ (S(IA), NOC)
C
C SET UP FOR AUTOMATICALLY DOING SOME PROPERTIES
C
      IDOAUT = IAUTO
      IACNT = 0
      if(.not.isther)then
       call figprp(varm,varic,NVARMX,noc)
      endif
      ICENTER = 0
25    CONTINUE
      M = 0
      IF (IDOAUT .NE. 0) THEN
        IACNT = IACNT + 1
        M = IDOPRP(1,IACNT)
        IF (M .EQ. 0) THEN
          IDOAUT = 0
          GOTO 25
        ENDIF
        ICT = IDOPRP(2,IACNT)
        C(1) = 0.
        C(2) = 0.
        C(3) = 0.
COLD    ILNMCM = 6H
        ILNMCM = 999
      ELSE
       if(isther)then
        READ(INUI,102,END=199) M,ILNMCM,ICT,(C(K),K=1,3)
       else
        iprpct=iprpct+1
        m=varm(iprpct)
        ict=varic(iprpct)
       endif
      ENDIF
C 
C This is a brain damaged algorithm. This whole program has to
C be rewritten, so I don't try to make it nicer. Ajith 08/93
C
      IF (M .EQ. 0) GOTO 199
      IF(ICT.LE.0) GO TO 126
      IX=IA+ICT-1
C
      DO I=1,3
      C(I)=S(IX)
      IX=IX+NOC
      ENDDO
C
 126  CONTINUE
COLD  IF (ICT .GE. 0 .AND. ILNMCM .EQ. 6H      )
      IF (ICT .GE. 0 .AND. ILNMCM .EQ. 999     )
     1  CALL GETCNT(ILNMCM, S(IA), NOC, S(IB), C)
C
C....    IF ICT IS -VE, THIS IS FOR EF OR EFG AND ALL CENTRES ARE DESIRE
C
      NPASSX = 1
      IF (ICT .LT. 0) NPASSX = NOC
300   FORMAT(T5,A,T27,F10.7,T39,F10.7,T51,F10.7)
302   FORMAT(T27,F10.7,T39,F10.7,T51,F10.7)
301   FORMAT(T5,A,T27,'   ----   ',T39,'   ----   ',T51,'   ----   ')
c300   FORMAT(/' PROPERTY TO BE COMPUTED ',3A6/
c     1 '     AT ',A6,'  ',3F15.8)
      IPASSY = 0
C
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
     &       20,21,22,23,24,25,26,27), M
1     CONTINUE
      WRITE(6,300)PRPTYP(M),(C(J),J=1,3)
c      WRITE(IOUTU) M,ILNMCM,0,(C(K),K=1,3)
      CALL DINTFM(OPAA0,POTNUC,S(IA),NOC,C,S(IC),S(ID),S(IK4),
     1 S(IK1),S(IK2),S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),
     2 S(IN),IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
      GO TO 25
2     CONTINUE
      IF (NPASSX .EQ. 1) THEN
      WRITE(6,300)PRPTYP(M),(C(J),J=1,3)
c        WRITE(IOUTU) M,ILNMCM,0,(C(K),K=1,3)
        CALL DINTFM(OPAA1,FLDNUC,S(IA),NOC,C,S(IC),S(ID),S(IK4),
     1   S(IK1),S(IK2),S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),
     2   S(IN),IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
      ELSE
        DO 42 ICENT = 1,NPASSX
          C(1) = S(IA+ICENT-1)
          C(2) = S(IA+ICENT+NOC-1)
          C(3) = S(IA+ICENT+2*NOC-1)
          ILNMCM = IS(IB+ICENT-1)
c          WRITE(IOUTU ) M,ILNMCM,ICENT,(C(K),K=1,3)
      WRITE(6,300)PRPTYP(M),(C(J),J=1,3)
          CALL DINTFM(OPAA1,FLDNUC,S(IA),NOC,C,S(IC),S(ID),S(IK4),
     1    S(IK1),S(IK2),S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),
     2    S(IN),IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
42      CONTINUE
      ENDIF
      GO TO 25
3     CONTINUE
      IF (NPASSX .EQ. 1) THEN
      WRITE(6,300)PRPTYP(M),(C(J),J=1,3)
c        WRITE(IOUTU) M,ILNMCM,0,(C(K),K=1,3)
        CALL DINTFM(OPAA2,GRDNUC,S(IA),NOC,C,S(IC),S(ID),S(IK4),
     1   S(IK1),S(IK2),S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),
     2   S(IN),IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
      ELSE
         DO 43 ICENT = 1,NPASSX
          C(1) = S(IA+ICENT-1)
          C(2) = S(IA+ICENT+NOC-1)
          C(3) = S(IA+ICENT+2*NOC-1)
          ILNMCM = IS(IB+ICENT-1)
c          WRITE(IOUTU ) M,ILNMCM,ICENT,(C(K),K=1,3)
      WRITE(6,300)PRPTYP(M),(C(J),J=1,3)
          CALL DINTFM(OPAA2,GRDNUC,S(IA),NOC,C,S(IC),S(ID),S(IK4),
     1     S(IK1),S(IK2),S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),
     2     S(IN),IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
43      CONTINUE
      ENDIF
      GO TO 25
4     CONTINUE
      WRITE(6,301)PRPTYP(M)
CSSS      WRITE(IOUTU) M,ILNMCM,0,(C(K),K=1,3)
      CALL DINTFM(OPAB1,DIPNUC,S(IA),NOC,C,S(IC),S(ID),S(IK4),
     1 S(IK1),S(IK2),S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),
     2 S(IN),IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
      IF (IPASSY .EQ. 1 .OR. IDOCOM(M) .LE. 0) GOTO 25
      IPASSY = 1
      C(1) = COM(1)
      C(2) = COM(2)
      C(3) = COM(3)
COLD  ILNMCM = 6HC-O-M
      ILNMCM = 999999
      IDOCOM(M) = 0
      GOTO 4
5     CONTINUE
      IF(IDOCOM(M).LE.0)WRITE(6,301)PRPTYP(M)
c      WRITE(IOUTU) M,ILNMCM,0,(C(K),K=1,3)
      CALL DINTFM(OPAB2,QUDNUC,S(IA),NOC,C,S(IC),S(ID),S(IK4),
     1 S(IK1),S(IK2),S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),
     2 S(IN),IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
      IF (IPASSY .EQ. 1 .OR. IDOCOM(M) .LE. 0) GOTO 25
      IPASSY = 1
      C(1) = COM(1)
      C(2) = COM(2)
      C(3) = COM(3)
COLD  ILNMCM = 6HC-O-M
      ILNMCM = 999999
      IDOCOM(M) = 0
      GOTO 5
6     CONTINUE
      IF(IDOCOM(M).LE.0)WRITE(6,301)PRPTYP(M)
c      WRITE(IOUTU) M,ILNMCM,0,(C(K),K=1,3)
      CALL DINTFM(OPAB3,OCTNUC,S(IA),NOC,C,S(IC),S(ID),S(IK4),
     1 S(IK1),S(IK2),S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),
     2 S(IN),IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
      IF (IPASSY .EQ. 1 .OR. IDOCOM(M) .LE. 0) GOTO 25
      IPASSY = 1
      C(1) = COM(1)
      C(2) = COM(2)
      C(3) = COM(3)
COLD  ILNMCM = 6HC-O-M
      ILNMCM = 999999
      IDOCOM(M) = 0
      GOTO 6
7     CONTINUE
      WRITE(6,300)PRPTYP(M),(C(J),J=1,3)
c      WRITE(IOUTU) M,ILNMCM,0,(C(K),K=1,3)
      CALL DINTFM(OPAC1, NONUC,S(IA),NOC,C,S(IC),S(ID),S(IK4),
     1 S(IK1),S(IK2),S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),
     2 S(IN),IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
      GO TO 25
8     CONTINUE
      WRITE(6,300)PRPTYP(M),(C(J),J=1,3)
c      WRITE(IOUTU) M,ILNMCM,0,(C(K),K=1,3)
      CALL DINTFM(OPAC2, NONUC,S(IA),NOC,C,S(IC),S(ID),S(IK4),
     1 S(IK1),S(IK2),S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),
     2 S(IN),IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
      GO TO 25
9     CONTINUE
      WRITE(6,300)PRPTYP(M),(C(J),J=1,3)
c      WRITE(IOUTU) M,ILNMCM,0,(C(K),K=1,3)
      CALL DINTFM(OPAC3, NONUC,S(IA),NOC,C,S(IC),S(ID),S(IK4),
     1 S(IK1),S(IK2),S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),
     2 S(IN),IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
      GO TO 25
10    CONTINUE
      WRITE(6,300)PRPTYP(M),(C(J),J=1,3)
c      WRITE(IOUTU) M,ILNMCM,0,(C(K),K=1,3)
      CALL DINTFM(OPAD1, NONUC,S(IA),NOC,C,S(IC),S(ID),S(IK4),
     1 S(IK1),S(IK2),S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),
     2 S(IN),IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
      GO TO 25
11    CONTINUE
      WRITE(6,300)PRPTYP(M),(C(J),J=1,3)
c      WRITE(IOUTU) M,ILNMCM,0,(C(K),K=1,3)
      CALL DINTFM(OPAB4,MAGNUC,S(IA),NOC,C,S(IC),S(ID),S(IK4),
     1 S(IK1),S(IK2),S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),
     2 S(IN),IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
      IF (IPASSY .EQ. 1 .OR. IDOCOM(M) .LE. 0) GOTO 25
      IPASSY = 1
      C(1) = COM(1)
      C(2) = COM(2)
      C(3) = COM(3)
COLD  ILNMCM = 6HC-O-M
      ILNMCM = 999999
      IDOCOM(M) = 0
      GOTO 11
12    CONTINUE
      WRITE(6,301)PRPTYP(M)
c      WRITE(IOUTU) M,ILNMCM,0,(C(K),K=1,3)
      CALL DINTFM(OPAD1, NONUC,S(IA),NOC,C,S(IC),S(ID),S(IK4),
     1 S(IK1),S(IK2),S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),
     2 S(IN),IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
      GO TO 25
13    CONTINUE
      GOTO 25
14    CONTINUE
      WRITE(6,301)PRPTYP(M)
c      WRITE(6,*) M,ILNMCM,0,(C(K),K=1,3)
      CALL DINTFM(OPAB5,SECNUC,S(IA),NOC,C,S(IC),S(ID),S(IK4),
     1 S(IK1),S(IK2),S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),
     2 S(IN),IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
      IF (IPASSY .EQ. 1 .OR. IDOCOM(M) .LE. 0) GOTO 25
      IPASSY = 1
      C(1) = COM(1)
      C(2) = COM(2)
      C(3) = COM(3)
COLD  ILNMCM = 6HC-O-M
      ILNMCM = 999999
      IDOCOM(M) = 0
      GOTO 14
15    CONTINUE
      WRITE(6,300)PRPTYP(M),(C(J),J=1,3)
c      WRITE(IOUTU) M,ILNMCM,0,(C(K),K=1,3)
      CALL DINTFM(OPAB6,THINUC,S(IA),NOC,C,S(IC),S(ID),S(IK4),
     1 S(IK1),S(IK2),S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),
     2 S(IN),IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
      IF (IPASSY .EQ. 1 .OR. IDOCOM(M) .LE. 0) GOTO 25
      IPASSY = 1
      C(1) = COM(1)
      C(2) = COM(2)
      C(3) = COM(3)
COLD  ILNMCM = 6HC-O-M
      ILNMCM = 999999
      IDOCOM(M) = 0
      GOTO 15
16    CONTINUE
      WRITE(6,301)PRPTYP(M)
c      WRITE(IOUTU) M,ILNMCM,0,(C(K),K=1,3)
      CALL DINTFM(OPAD1, NONUC,S(IA),NOC,C,S(IC),S(ID),S(IK4),
     1 S(IK1),S(IK2),S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),
     2 S(IN),IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
      GOTO 25
17    CONTINUE
C
C *****************************************************************
C * Orbital diamagnetic contribution to the NMR coupling constant *
C * O. Matsuoka and T. Aoyama J. Chem. Phys. 73(11), 5719, 1980.  *
C * Implemented by Ajith 08/93.                                   *
C *****************************************************************
C  
C  Set up to do Gauss-Legendre Numerical integration
C  At Present number of integration points is fixed to
C  24 (defined as a parameter NGauss at the begining).
C  Get the roots and weights and store them in arrays Root
C  and Weight also defined at the begining.
C      
      Call GaussL (OneM, One, Root, Weight, NGauss)
      DO 37 I = 1, 3
         Centr1(I) = C(I)
 37   CONTINUE
C     
      IFirst = ICT + 1
C
C Comply with the condition in <|X_C*X_D/R(3)_C*R(3)_D|> C Not equal D.
C Loop over distinct pair of centers
C
      DO 1005 ICENT = IFirst, NOC
         Centr2(1) = S(IA + ICENT - 1)
         Centr2(2) = S(IA + ICENT + NOC - 1)
         Centr2(3) = S(IA + ICENT + 2*NOC - 1) 
C
C******DEBUG
C      Centr1(1) = -0.4
C      Centr1(2) = -0.3
C      Centr1(3) =  0.5
C
C      Centr2(1) =  0.3
C      Centr2(2) =  0.5
C      Centr2(3) = -0.6
C
C******DEBUG
C      Centr2(1) = -0.4
C      Centr2(2) = -0.3
C      Centr2(3) =  0.5
C
C      Centr1(1) =  0.3
C      Centr1(2) =  0.5
C      Centr1(3) = -0.6
C
         WRITE(6,300) PRPTYP(M), (CENTR1(J), J = 1, 3)
         WRITE(6,302) (CENTR2(J), J = 1, 3)
C
         Call OrbDia(S(IA),NOC,Centr1,Centr2,S(IC),S(ID),S(IK4),S(IK1),
     &               S(IK2),S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),
     &               S(IN),IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2,
     &               Root,Weight,ValueM,NGauss)
         

 1005 CONTINUE
C      
      GOTO 25
C     
 18   CONTINUE
C     
      IF (NPASSX .EQ. 1) THEN
         WRITE(6,300) PRPTYP(M), (C(J), J = 1, 3)
C
         Call OrbPara(S(IA),NOC,C,S(IC),S(ID),S(IK4),S(IK1),S(IK2),
     &               S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),S(IN),
     &               IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
      ELSE
         DO 44 ICENT = 1, NPASSX
            C(1) = S(IA + ICENT - 1)
            C(2) = S(IA + ICENT + NOC - 1)
            C(3) = S(IA + ICENT + 2*NOC - 1)
            ILNMCM = IS(IB + ICENT - 1)
            WRITE(6,300) PRPTYP(M),(C(J), J=1, 3)
C
            Call OrbPara(S(IA),NOC,C,S(IC),S(ID),S(IK4),S(IK1),S(IK2),
     &                   S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),S(IN),
     &                   IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
 
 44      CONTINUE
      ENDIF
C
      GOTO 25
C
 19   CONTINUE
C     
      IF (NPASSX .EQ. 1) THEN
         WRITE(6,300) PRPTYP(M), (C(J), J = 1, 3)
C
         Call SpinDip(S(IA),NOC,C,S(IC),S(ID),S(IK4),S(IK1),S(IK2),
     &               S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),S(IN),
     &               IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
      ELSE
         DO 66 ICENT = 1, NPASSX
            C(1) = S(IA + ICENT - 1)
            C(2) = S(IA + ICENT + NOC - 1)
            C(3) = S(IA + ICENT + 2*NOC - 1)
            ILNMCM = IS(IB + ICENT - 1)
            WRITE(6,300) PRPTYP(M),(C(J), J=1, 3)
C
            Call SpinDip(S(IA),NOC,C,S(IC),S(ID),S(IK4),S(IK1),S(IK2),
     &                   S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),S(IN),
     &                   IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
 
 66      CONTINUE
      ENDIF
C
      GOTO 25
C
 20   CONTINUE
C     
CSSS      IF (NPASSX .EQ. 1) THEN
CSSS         WRITE(6,300) PRPTYP(M), (C(J), J = 1, 3)
C
CSSS         Call DiaShil(S(IA),NOC,C,S(IC),S(ID),S(IK4),S(IK1),S(IK2),
CSSS     &                S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),S(IN),
CSSS     &                IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
CSSS      ELSE
         NPASSX = 3
         DO 67 ICENT = 1, NPASSX
            C(1) = S(IA + ICENT - 1)
            C(2) = S(IA + ICENT + NOC - 1)
            C(3) = S(IA + ICENT + 2*NOC - 1)
            ILNMCM = IS(IB + ICENT - 1)
            WRITE(6,300) PRPTYP(M),(C(J), J=1, 3)
C
            Call DiaShil(S(IA),NOC,C,S(IC),S(ID),S(IK4),S(IK1),S(IK2),
     &                  S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),S(IN),
     &                  IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)

 67      CONTINUE
CSSSS      ENDIF
C
      GO TO 25
C
 21   CONTINUE
C
      IF (NPASSX .EQ. 1) THEN
         WRITE(6,300) PRPTYP(M), (C(J), J = 1, 3)
C
      Call SecMOM(S(IA),NOC,C,S(IC),S(ID),S(IK4),S(IK1),S(IK2),
     &             S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),S(IN),
     &             IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
      ELSE
         DO 68 ICENT = 1, NPASSX
            C(1) = S(IA + ICENT - 1)
            C(2) = S(IA + ICENT + NOC - 1)
            C(3) = S(IA + ICENT + 2*NOC - 1)
            ILNMCM = IS(IB + ICENT - 1)
            WRITE(6,300) PRPTYP(M),(C(J), J=1, 3)
C
            Call SECMOM(S(IA),NOC,C,S(IC),S(ID),S(IK4),S(IK1),S(IK2),
     &                  S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),S(IN),
     &                  IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
 
 68      CONTINUE
      ENDIF
C
      GO TO 25
C
 22   CONTINUE
C     
      WRITE(6,301) PRPTYP(M)
C     
      Call DipVelo(S(IA),NOC,C,S(IC),S(ID),S(IK4),S(IK1),S(IK2),
     &             S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),S(IN),
     &             IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
C 
      GO TO 25
C
 23   CONTINUE
C
      IF (NPASSX .EQ. 1) THEN
         WRITE(6,300) PRPTYP(M), (C(J), J = 1, 3)
         ICENTER =  ICENTER + 1 
         CHARGE = S(3*NOC + ICENTER)

         Call GraDPotnl(S(IA),NOC,C,S(IC),S(ID),S(IK4),S(IK1),
     &                 S(IK2),S(IG),S(IH),S(II),NCONT,S(IKS1),
     &                 S(IKS1),S(IN),IBUF,S(IK3),NGROUP,S(IK5),
     &                 S(IK6),CHARGE,MXPRM2)
      ELSE
         DO ICENT = 1, NPASSX
            C(1) = S(IA + ICENT - 1)
            C(2) = S(IA + ICENT + NOC - 1)
            C(3) = S(IA + ICENT + 2*NOC - 1)
            LNMCM = IS(IB + ICENT - 1)
            WRITE(6,300) PRPTYP(M),(C(J), J=1, 3)
            ICENTER =  ICENTER + 1
            CHARGE = S(3*NOC + ICENTER )

            Call GraDPotnl(S(IA),NOC,C,S(IC),S(ID),S(IK4),S(IK1),
     &                     S(IK2),S(IG),S(IH),S(II),NCONT,S(IKS1),
     &                     S(IKS1),S(IN),IBUF,S(IK3),NGROUP,S(IK5),
     &                     S(IK6),CHARGE,MXPRM2)
         ENDDO 
      ENDIF
C
      GO TO 25
C     
 24   CONTINUE

      IF (NPASSX .EQ. 1) THEN
         WRITE(6,301) PRPTYP(M)
C    
         Call AngMom(S(IA),NOC,C,S(IC),S(ID),S(IK4),S(IK1),S(IK2),
     &               S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),S(IN),
     &               IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
      ELSE  
         DO ICENT = 1, NPASSX
            C(1) = S(IA + ICENT - 1)
            C(2) = S(IA + ICENT + NOC - 1)
            C(3) = S(IA + ICENT + 2*NOC - 1)
            ILNMCM = IS(IB + ICENT - 1)
            WRITE(6,300) PRPTYP(M),(C(J), J=1, 3)
C    
            Call AngMom(S(IA),NOC,C,S(IC),S(ID),S(IK4),S(IK1),S(IK2),
     &                  S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),S(IN),
     &                  IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
      
         ENDDO
      ENDIF 
C
      GOTO 25

 26   CONTINUE 

      IF (NPASSX .EQ. 1) THEN
         WRITE(6,301) PRPTYP(M)
C   
         Call Magqmom(S(IA),NOC,C,S(IC),S(ID),S(IK4),S(IK1),S(IK2),
     &                S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),S(IN),
     &                IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)
      ELSE
         DO ICENT = 1, NPASSX
            C(1) = S(IA + ICENT - 1)
            C(2) = S(IA + ICENT + NOC - 1)
            C(3) = S(IA + ICENT + 2*NOC - 1)
            ILNMCM = IS(IB + ICENT - 1)
            WRITE(6,300) PRPTYP(M),(C(J), J=1, 3)
C   
            Call Magqmom(S(IA),NOC,C,S(IC),S(ID),S(IK4),S(IK1),S(IK2),
     &                   S(IG),S(IH),S(II),NCONT,S(IKS1),S(IKS1),S(IN),
     &                   IBUF,S(IK3),NGROUP,S(IK5),S(IK6),MXPRM2)

         ENDDO
      ENDIF
C
      GOTO 25

 27   CONTINUE
C
C This take care of the requirement of the VPINP file. 
C But this is ugly way of doing it, but the whole program
C is dirty.
C
      GO TO 25
C
199   CONTINUE
C
C      WRITE(IOUTU) M,ILNMCM,0,(C(K),K=1,3)
 200  CONTINUE
C
      WRITE(6,152)
C
      close(unit=id20,status='delete')
      call aces_fin
C
      stop

100   FORMAT(6I5)
101   FORMAT(12A6)
1001  FORMAT(///,20X,'MOLECULAR PROPERTIES PROGRAM',//,
     1           40X,'Cartesian or spherical harmonic basis',/
     2           40X,'Up to i-type functions',
     3           50X,'January  1989',//)
102   FORMAT(I5,4X,A6,I5,3F20.12)
103   FORMAT(//)
      end
