c--------------------------------------------------------------------------
c The modules below are extracted from the GAMESS source with the following
c modifications:
c 1. 'cwc' string used to comment out unnecessary lines
c 2. 'cwc1' and 'cwc0' denote begin and end of added fragments
c 3. all external names have 'x' added in front
c------------------------------------------------------------------------
C*MODULE INT2A   *DECK DEBUT
      SUBROUTINE xDEBUT(INTG76,DIRSCF,BUFP,BUFK,IX,NINTMX,NEED,DIRTRF)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL DIRSCF,DIRTRF
      LOGICAL DIR,GOPARR,DSKWRK,MASWRK,PK,PANDK,BLOCK,OUT
C
      DIMENSION BUFP(NINTMX),BUFK(NINTMX),IX(*)
C
#include "gamess.par"
C
      COMMON /xB     / CO(MXSH,3)
      COMMON /xIJPAIR/ IA(MXAO)
      COMMON /xINFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM)
      COMMON /xINTPR / QINT(3),VALINT(3),JCINT(16)
      COMMON /xIOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /xMACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     2                CF(MXGTOT),CG(MXGTOT),
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /xOUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /xPAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /xPCKLAB/ LABSIZ
      COMMON /xPKFIL / PK,PANDK,BLOCK
      COMMON /xRESTAR/ TIMLIM,IREST,NREC,INTLOC,IST,JST,KST,LST
      COMMON /xRUNOPT/ RUNTYP,EXETYP,NEVALS
      COMMON /xSHLT  / TOL,CUTOFF,OUT,ICOUNT
C
      PARAMETER (ONE=1.0D+00, TEN=10.D+00, RLN10=2.30258D+00)
C
      DATA CHECK/8HCHECK   /
C
C     ----- INITIALIZE TWO ELECTRON INTEGRAL CALCULATION -----
C
      TOL= ITOL*RLN10
      CUTOFF= ONE/(TEN**ICUT)
      DIR= DIRSCF .OR. DIRTRF
C
      DO 100 I=1,NUM
         IA(I)=(I*I-I)/2
  100 CONTINUE
C
C     ----- READ ERROR FUNCTION INTERPOLATION TABLE -----
C
      IF(INTG76.EQ.1 .AND. EXETYP.NE.CHECK) CALL xGAMGEN(2)
C
      IF((NPRINT.NE.-5) .AND. MASWRK) THEN
         IF(DIR) THEN
            IF(DIRSCF) WRITE(IW,9000)
cwc            IF(DIRSCF) WRITE(IW,9040)
cwc            IF(DIRTRF) WRITE(IW,9045)
         ELSE
            WRITE(IW,9000)
            IF(PK) THEN
               IF(     PANDK .AND. NPRINT.NE.-5) WRITE(IW,9010)
               IF(.NOT.PANDK) WRITE(IW,9020)
            ELSE
               WRITE(IW,9030)
            END IF
            NBYTES = 8 + (LABSIZ*8)/NWDVAR
            IF(PANDK) NBYTES = NBYTES+8
            WRITE(IW,9050) NINTMX,NBYTES
            WRITE(IW,9060) NEED
         END IF
      END IF
C
      OUT = NPRINT.EQ.4 .AND. MASWRK
      JCINT(1) = 0
C
      DO 200 I=1,NSHELL
         ICC = KATOM(I)
         CO(I,1)= C(1,ICC)
         CO(I,2)= C(2,ICC)
         CO(I,3)= C(3,ICC)
  200 CONTINUE
C
      IF(DIRSCF .OR. DIRTRF .OR. EXETYP.EQ.CHECK) GO TO 400
C
cwc      CALL SEQREW(IS)
C
      IF(IREST.LT.1 .OR. NREC.LE.1 .OR. INTLOC.LE.1) GO TO 400
  300 CONTINUE
C
C     ----- POSITION THE INTEGRAL FILE FOR A RESTART JOB -----
C
      ICOUNT = INTLOC
      N = NREC-1
      IF(MASWRK) THEN
         DO 310 I=1,N
            READ(IS)
  310    CONTINUE
      END IF
cwc      IF(.NOT.PANDK) CALL PREAD (IS,BUFP,IX,NXX,NINTMX)
cwc      IF(     PANDK) CALL PKREAD(IS,BUFP,BUFK,IX,NXX,NINTMX)
cwc      CALL SEQREW(IS)
      IF(MASWRK) THEN
         DO 320 I=1,N
            READ(IS)
  320    CONTINUE
      END IF
      GO TO 999
  400 CONTINUE
C
C     ----- NORMAL START -----
C
      IF(IST.LT.1) IST = 1
      IF(JST.LT.1) JST = 1
      IF(KST.LT.1) KST = 1
      IF(LST.LT.1) LST = 1
      IF(IST.GT.NSHELL) GO TO 999
      IF(IST.NE.1 .OR. JST.NE.1 .OR. KST.NE.1 .OR. LST.NE.1) GO TO 300
      NREC   = 1
      INTLOC = 1
      ICOUNT = 1
  999 CONTINUE
      RETURN
 9000 FORMAT(/10X,20(1H-)/10X,'2 ELECTRON INTEGRALS'/10X,20(1H-)/)
 9010 FORMAT(' THE -PK- OPTION IS ON, CREATING -P- AND -K-',
     *       ' SUPERMATRICES.')
 9020 FORMAT(' THE -PK- OPTION IS ON, CREATING A -P- SUPERMATRIX.')
 9030 FORMAT(' THE -PK- OPTION IS OFF, THE INTEGRALS ARE NOT IN',
     *       ' SUPERMATRIX FORM.')
 9040 FORMAT(' DIRECT SCF METHOD SKIPS INTEGRAL STORAGE ON DISK.')
 9045 FORMAT(' DIRECT TRANSFORMATION SKIPS AO INTEGRAL STORAGE',
     *       ' ON DISK.')
 9050 FORMAT(' STORING',I8,' INTEGRALS/RECORD ON DISK, USING',I3,
     *       ' BYTES/INTEGRAL.')
 9060 FORMAT(' TWO ELECTRON INTEGRAL EVALUATION REQUIRES',I8,
     *       ' WORDS OF MEMORY.')
      END
c------------------------------------------------------------------------
C*MODULE INT2A   *DECK SHELLS
cwc      SUBROUTINE xSHELLS(NELEC,ISH,JSH,KSH,LSH,FLIP)
      SUBROUTINE xSHELLS(NELEC,ISH,JSH,KSH,LSH,FLIP,nosame)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL FLIP
      LOGICAL IANDJ,KANDL,SAME,FIRST
C
      DIMENSION IX(35),IY(35),IZ(35),
     *          JX(35),JY(35),JZ(35),
     *          KX(35),KY(35),KZ(35),
     *          LX(35),LY(35),LZ(35)
C
#include "gamess.par"
C
      COMMON /xGOUT  / GPOPLE(768),NORG
      COMMON /xINTDEX/ IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784),
     *                KLGT(784),KLX(784),KLY(784),KLZ(784)
      COMMON /xINFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM)
      COMMON /xMISC  / IANDJ,KANDL,SAME
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /xROOT  / XX,U(13),W(13),NROOTS
      COMMON /xSHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     +                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     +                NIJ,IJ,KL,IJKL
      COMMON /xSHLINF/ GA(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH),
     *                CFA(MXGSH),CGA(MXGSH),
     *                GB(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH),
     *                CFB(MXGSH),CGB(MXGSH),
     *                GC(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH),
     *                CFC(MXGSH),CGC(MXGSH),
     *                GD(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH),
     *                CFD(MXGSH),CGD(MXGSH),
     *                AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD,
     *                NGA,NGB,NGC,NGD
cwc1
      common /shtrue/ inu,jnu,knu,lnu
      logical nosame
cwc0
C
      SAVE FIRST,IGT,JGT,KGT,LGT
C
      DATA LX /   0,  1,  0,  0,  2,  0,  0,  1,  1,  0,
     *            3,  0,  0,  2,  2,  1,  0,  1,  0,  1,
     *            4,  0,  0,  3,  3,  1,  0,  1,  0,  2,
     *            2,  0,  2,  1,  1/
      DATA KX /   0,  5,  0,  0, 10,  0,  0,  5,  5,  0,
     *           15,  0,  0, 10, 10,  5,  0,  5,  0,  5,
     *           20,  0,  0, 15, 15,  5,  0,  5,  0, 10,
     *           10,  0, 10,  5,  5/
      DATA JX /   0, 25,  0,  0, 50,  0,  0, 25, 25,  0,
     *           75,  0,  0, 50, 50, 25,  0, 25,  0, 25,
     *          100,  0,  0, 75, 75, 25,  0, 25,  0, 50,
     *           50,  0, 50, 25, 25/
      DATA IX /   1,126,  1,  1,251,  1,  1,126,126,  1,
     *          376,  1,  1,251,251,126,  1,126,  1,126,
     *          501,  1,  1,376,376,126,  1,126,  1,251,
     *          251,  1,251,126,126/
      DATA LY /   0,  0,  1,  0,  0,  2,  0,  1,  0,  1,
     *            0,  3,  0,  1,  0,  2,  2,  0,  1,  1,
     *            0,  4,  0,  1,  0,  3,  3,  0,  1,  2,
     *            0,  2,  1,  2,  1/
      DATA KY /   0,  0,  5,  0,  0, 10,  0,  5,  0,  5,
     *            0, 15,  0,  5,  0, 10, 10,  0,  5,  5,
     *            0, 20,  0,  5,  0, 15, 15,  0,  5, 10,
     *            0, 10,  5, 10,  5/
      DATA JY /   0,  0, 25,  0,  0, 50,  0, 25,  0, 25,
     *            0, 75,  0, 25,  0, 50, 50,  0, 25, 25,
     *            0,100,  0, 25,  0, 75, 75,  0, 25, 50,
     *            0, 50, 25, 50, 25/
      DATA IY /   1,  1,126,  1,  1,251,  1,126,  1,126,
     *            1,376,  1,126,  1,251,251,  1,126,126,
     *            1,501,  1,126,  1,376,376,  1,126,251,
     *            1,251,126,251,126/
      DATA LZ /   0,  0,  0,  1,  0,  0,  2,  0,  1,  1,
     *            0,  0,  3,  0,  1,  0,  1,  2,  2,  1,
     *            0,  0,  4,  0,  1,  0,  1,  3,  3,  0,
     *            2,  2,  1,  1,  2/
      DATA KZ /   0,  0,  0,  5,  0,  0, 10,  0,  5,  5,
     *            0,  0, 15,  0,  5,  0,  5, 10, 10,  5,
     *            0,  0, 20,  0,  5,  0,  5, 15, 15,  0,
     *           10, 10,  5,  5, 10/
      DATA JZ /   0,  0,  0, 25,  0,  0, 50,  0, 25, 25,
     *            0,  0, 75,  0, 25,  0, 25, 50, 50, 25,
     *            0,  0,100,  0, 25,  0, 25, 75, 75,  0,
     *           50, 50, 25, 25, 50/
      DATA IZ /   1,  1,  1,126,  1,  1,251,  1,126,126,
     *            1,  1,376,  1,126,  1,126,251,251,126,
     *            1,  1,501,  1,126,  1,126,376,376,  1,
     *          251,251,126,126,251/
      DATA FIRST/.TRUE./
C
      IF(FIRST) THEN
         FIRST=.FALSE.
         CALL xBASCHK(LMAX)
                       NANGM =  4
         IF(LMAX.EQ.2) NANGM =  6
         IF(LMAX.EQ.3) NANGM = 10
         IF(LMAX.EQ.4) NANGM = 15
         LGT = 1
         KGT = LGT * NANGM
         JGT = KGT * NANGM
         IGT = JGT * NANGM
      END IF
C
C     PREPARE SHELL INFORMATION/FOR HONDO INTEGRATION
C
      IF(NELEC.EQ.2) GO TO 200
C
C     ----- PERMUTE ISH AND JSH SHELLS, FOR THEIR TYPE
C     THIS IS DONE FOR SPEED REASONS.  THE CODE GETS THE RIGHT ANSWER
C     WITHOUT THE ANGULAR MOMENTUM FLIPPING, AND THEREFORE A CALLING
C     ARGUMENT ALLOWS ONE DO EXACTLY THE INTEGRAL BLOCK AS SPECIFIED,
C     SHOULD THAT BE DESIRED.
C
      IANDJ = ISH .EQ. JSH
      IF (KTYPE(ISH) .LT. KTYPE(JSH)  .AND.  FLIP) THEN
         INU = JSH
         JNU = ISH
         NGTI = JGT
         NGTJ = IGT
      ELSE
         INU = ISH
         JNU = JSH
         NGTI = IGT
         NGTJ = JGT
      END IF
C
C     ----- ISHELL
C
      I = KATOM(INU)
      AX = C(1,I)
      AY = C(2,I)
      AZ = C(3,I)
      I1 = KSTART(INU)
      I2 = I1+KNG(INU)-1
      LIT = KTYPE(INU)
      MINI = KMIN(INU)
      MAXI = KMAX(INU)
      LOCI = KLOC(INU)-MINI
      NGA = 0
      DO 140 I = I1,I2
         NGA = NGA+1
         GA(NGA) = EX(I)
         CSA(NGA) = CS(I)
         CPA(NGA) = CP(I)
         CDA(NGA) = CD(I)
         CFA(NGA) = CF(I)
         CGA(NGA) = CG(I)
  140 CONTINUE
C
C     ----- JSHELL
C
      J = KATOM(JNU)
      BX = C(1,J)
      BY = C(2,J)
      BZ = C(3,J)
      J1 = KSTART(JNU)
      J2 = J1+KNG(JNU)-1
      LJT = KTYPE(JNU)
      MINJ = KMIN(JNU)
      MAXJ = KMAX(JNU)
      LOCJ = KLOC(JNU)-MINJ
      NGB = 0
      DO 160 J = J1,J2
         NGB = NGB+1
         GB(NGB) = EX(J)
         CSB(NGB) = CS(J)
         CPB(NGB) = CP(J)
         CDB(NGB) = CD(J)
         CFB(NGB) = CF(J)
         CGB(NGB) = CG(J)
  160 CONTINUE
      RAB = ((AX-BX)*(AX-BX) + (AY-BY)*(AY-BY) + (AZ-BZ)*(AZ-BZ))
C
C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS
C
      IJ = 0
      JMAX = MAXJ
      DO 190 I = MINI,MAXI
         NX = IX(I)
         NY = IY(I)
         NZ = IZ(I)
         IF (IANDJ) JMAX = I
         DO 180 J = MINJ,JMAX
            IJ = IJ+1
            IJX(IJ) = NX+JX(J)
            IJY(IJ) = NY+JY(J)
            IJZ(IJ) = NZ+JZ(J)
            IJGT(IJ) = NGTI*(I-MINI)+NGTJ*(J-MINJ)+1+NORG
  180    CONTINUE
  190 CONTINUE
      RETURN
C     ******
C
C        K AND L SHELL
C
  200 CONTINUE
      KANDL = KSH .EQ. LSH
      SAME = ISH .EQ. KSH .AND. JSH .EQ. LSH
cwc1
      if (nosame) same=.false.
cwc0

C
C     ----- PERMUTE KSH AND LSH SHELLS, FOR THEIR TYPE
C
      IF (KTYPE(KSH) .LT. KTYPE(LSH)  .AND.  FLIP) THEN
         KNU = LSH
         LNU = KSH
         NGTK = LGT
         NGTL = KGT
      ELSE
         KNU = KSH
         LNU = LSH
         NGTK = KGT
         NGTL = LGT
      END IF
C
C     ----- K SHELL
C
      K = KATOM(KNU)
      CX = C(1,K)
      CY = C(2,K)
      CZ = C(3,K)
      K1 = KSTART(KNU)
      K2 = K1+KNG(KNU)-1
      LKT = KTYPE(KNU)
      MINK = KMIN(KNU)
      MAXK = KMAX(KNU)
      LOCK = KLOC(KNU)-MINK
      NGC = 0
      DO 260 K = K1,K2
         NGC = NGC+1
         GC(NGC) = EX(K)
         CSC(NGC) = CS(K)
         CPC(NGC) = CP(K)
         CDC(NGC) = CD(K)
         CFC(NGC) = CF(K)
         CGC(NGC) = CG(K)
  260 CONTINUE
C
C     ----- LSHELL
C
      L = KATOM(LNU)
      DX = C(1,L)
      DY = C(2,L)
      DZ = C(3,L)
      L1 = KSTART(LNU)
      L2 = L1+KNG(LNU)-1
      LLT = KTYPE(LNU)
      MINL = KMIN(LNU)
      MAXL = KMAX(LNU)
      LOCL = KLOC(LNU)-MINL
      NGD = 0
      DO 280 L = L1,L2
         NGD = NGD+1
         GD(NGD) = EX(L)
         CSD(NGD) = CS(L)
         CPD(NGD) = CP(L)
         CDD(NGD) = CD(L)
         CFD(NGD) = CF(L)
         CGD(NGD) = CG(L)
  280 CONTINUE
      NROOTS = (LIT+LJT+LKT+LLT-2)/2
      RCD = ((CX-DX)*(CX-DX) + (CY-DY)*(CY-DY) + (CZ-DZ)*(CZ-DZ))
C
C     ----- PREPARE INDICES FOR PAIRS OF (K,L) FUNCTIONS
C
      KL = 0
      LMAX = MAXL
      DO 310 K = MINK,MAXK
         NX = KX(K)
         NY = KY(K)
         NZ = KZ(K)
         IF (KANDL) LMAX = K
         DO 300 L = MINL,LMAX
            KL = KL+1
            KLX(KL) = NX+LX(L)
            KLY(KL) = NY+LY(L)
            KLZ(KL) = NZ+LZ(L)
            KLGT(KL) = NGTK*(K-MINK)+NGTL*(L-MINL)
  300    CONTINUE
  310 CONTINUE
      MAX = KL
      DO 320 I = 1,IJ
      IF (SAME) MAX = I
  320 IK(I) = MAX
      IJKL = IJ*KL
      IF (SAME) IJKL = IJ*(IJ+1)/2
      RETURN
      END
c-------------------------------------------------------------------
C*MODULE INT2A   *DECK IJPRIM
      SUBROUTINE xIJPRIM(DDIJ)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL IANDJ,KANDL,SAME,OUT,NORM
C
#include "gamess.par"
C
      DIMENSION DDIJ(16*MXG2)
C
      COMMON /xIJGNRL/ A(MXG2),R(MXG2),X1(MXG2),Y1(MXG2),Z1(MXG2),
     *                IJD(784)
      COMMON /xMISC  / IANDJ,KANDL,SAME
      COMMON /xOUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /xSHLINF/ AG(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH),
     *                CFA(MXGSH),CGA(MXGSH),
     *                BG(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH),
     *                CFB(MXGSH),CGB(MXGSH),
     *                CG(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH),
     *                CFC(MXGSH),CGC(MXGSH),
     *                DG(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH),
     *                CFD(MXGSH),CGD(MXGSH),
     *                XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK,
     *                NGA,NGB,NGC,NGD
      COMMON /xSHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     *                NIJ,IJ,KL,IJKL
      COMMON /xSHLT  / TOL,CUTOFF,OUT,ICOUNT
C
      PARAMETER (SQRT3=1.73205080756888D+00,
     *           SQRT5=2.23606797749979D+00,
     *           SQRT7=2.64575131106459D+00,
     *           ZERO=0.0D+00, ONE=1.0D+00)
C

c      write (99,*) 'ijprim:',loci,locj,tol
      NORM = NORMF .NE. 1 .OR. NORMP .NE. 1
      MAX = MAXJ
      N = 0
      NN = 0
      NM = -2**20
      DO 180 I = MINI,MAXI
         GO TO (100,100,120,120,100,120,120,100,120,120,
     1          100,120,120,100,120,120,120,120,120,100,
     1          100,120,120,100,120,120,120,120,120,100,
     1          120,120,100,120,120),I
  100    NM = NN
  120    NN = NM
         IF (IANDJ) MAX = I
         DO 170 J = MINJ,MAX
            GO TO (140,140,160,160,140,160,160,140,160,160,
     1             140,160,160,140,160,160,160,160,160,140,
     1             140,160,160,140,160,160,160,160,160,140,
     1             160,160,140,160,160),J
  140       NN = NN+1
  160       N = N+1
            IJD(N) = NN
  170    CONTINUE
  180 CONTINUE
C
C     ----- I PRIMITIVE
C
      NIJ = 0
      JBMAX = NGB
      DO 540 IA = 1,NGA
         AI = AG(IA)
         ARRI = AI*RRI
         AXI = AI*XI
         AYI = AI*YI
         AZI = AI*ZI
         CSI = CSA(IA)
         CPI = CPA(IA)
         CDI = CDA(IA)
         CFI = CFA(IA)
         CGI = CGA(IA)
C
C        ----- J PRIMITIVE
C
         IF (IANDJ) JBMAX = IA
         DO 520 JB = 1,JBMAX
            AJ = BG(JB)
            AA = AI+AJ
            AAINV = ONE/AA
            DUM = AJ*ARRI*AAINV
c            write (99,*) ia,jb,dum,tol,nij

            IF (DUM .GT. TOL) GO TO 520
            CSJ = CSB(JB)
            CPJ = CPB(JB)
            CDJ = CDB(JB)
            CFJ = CFB(JB)
            CGJ = CGB(JB)
            NM = 16*NIJ
            NN = NM
            NIJ = NIJ+1
            R(NIJ) = DUM
            A(NIJ) = AA
            X1(NIJ) = (AXI+AJ*XJ)*AAINV
            Y1(NIJ) = (AYI+AJ*YJ)*AAINV
            Z1(NIJ) = (AZI+AJ*ZJ)*AAINV
C
C           ----- DENSITY FACTOR
C
            DUM1 = ZERO
            DUM2 = ZERO
            DO 420 I = MINI,MAXI
               GO TO (200,220,420,420,240,420,420,260,420,420,
     1                261,420,420,262,420,420,420,420,420,263,
     1                264,420,420,265,420,420,420,420,420,266,
     1                420,420,267,420,420),I
  200          DUM1 = CSI*AAINV
               GO TO 280
  220          DUM1 = CPI*AAINV
               GO TO 280
  240          DUM1 = CDI*AAINV
               GO TO 280
  260          IF (NORM) DUM1 = DUM1*SQRT3
               GO TO 280
  261          DUM1 = CFI*AAINV
               GO TO 280
  262          IF (NORM) DUM1 = DUM1*SQRT5
               GO TO 280
  263          IF (NORM) DUM1 = DUM1*SQRT3
               GO TO 280
  264          DUM1 = CGI*AAINV
               GO TO 280
  265          IF (NORM) DUM1 = DUM1*SQRT7
               GO TO 280
  266          IF (NORM) DUM1 = DUM1*SQRT5/SQRT3
               GO TO 280
  267          IF (NORM) DUM1 = DUM1*SQRT3
  280          IF (IANDJ) MAX = I
               DO 400 J = MINJ,MAX
                  GO TO (300,320,400,400,340,400,400,360,400,400,
     1                   361,400,400,362,400,400,400,400,400,363,
     1                   364,400,400,365,400,400,400,400,400,366,
     1                   400,400,367,400,400),J
  300             DUM2 = DUM1*CSJ
                  GO TO 380
  320             DUM2 = DUM1*CPJ
                  GO TO 380
  340             DUM2 = DUM1*CDJ
                  GO TO 380
  360             IF (NORM) DUM2 = DUM2*SQRT3
                  GO TO 380
  361             DUM2 = DUM1*CFJ
                  GO TO 380
  362             IF (NORM) DUM2 = DUM2*SQRT5
                  GO TO 380
  363             IF (NORM) DUM2 = DUM2*SQRT3
                  GO TO 380
  364             DUM2 = DUM1*CGJ
                  GO TO 380
  365             IF (NORM) DUM2 = DUM2*SQRT7
                  GO TO 380
  366             IF (NORM) DUM2 = DUM2*SQRT5/SQRT3
                  GO TO 380
  367             IF (NORM) DUM2 = DUM2*SQRT3
  380             NN = NN+1
                  DDIJ(NN) = DUM2
  400          CONTINUE
  420       CONTINUE
            IF ( .NOT. IANDJ) GO TO 520
            IF (IA .EQ. JB) GO TO 520
            GO TO (500,440,460,455,450),LIT
  440       IF (MINI .EQ. 2) GO TO 500
            DDIJ(NM+2) = DDIJ(NM+2)+CSI*CPJ*AAINV
            GO TO 480
  450       DDIJ(NM+10) = DDIJ(NM+10)+DDIJ(NM+10)
            DDIJ(NM+9) = DDIJ(NM+9)+DDIJ(NM+9)
            DDIJ(NM+8) = DDIJ(NM+8)+DDIJ(NM+8)
            DDIJ(NM+7) = DDIJ(NM+7)+DDIJ(NM+7)
  455       DDIJ(NM+6) = DDIJ(NM+6)+DDIJ(NM+6)
            DDIJ(NM+5) = DDIJ(NM+5)+DDIJ(NM+5)
            DDIJ(NM+4) = DDIJ(NM+4)+DDIJ(NM+4)
  460       DDIJ(NM+2) = DDIJ(NM+2)+DDIJ(NM+2)
  480       DDIJ(NM+3) = DDIJ(NM+3)+DDIJ(NM+3)
  500       DDIJ(NM+1) = DDIJ(NM+1)+DDIJ(NM+1)
  520    CONTINUE
  540 CONTINUE
c      write (99,*) 'ijprim end'
      RETURN
      END
c------------------------------------------------------------------
C*MODULE INT2A   *DECK GENRAL
      SUBROUTINE xGENRAL(GHONDO,DDIJ)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION GHONDO(*),DDIJ(*)
C
      LOGICAL IANDJ,KANDL,SAME,OUT,NORM,DOUBLE
C
#include "gamess.par"
C
      COMMON /xDENS  / DKL(784),DIJ(784)
      COMMON /xIJGNRL/ AA(MXG2),R(MXG2),X1(MXG2),Y1(MXG2),Z1(MXG2),
     *                IJD(784)
      COMMON /xMISC  / IANDJ,KANDL,SAME
      COMMON /xOUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /xROOT  / XX,U(13),W(13),NROOTS
      COMMON /xSETINT/ IN(13),KN(13),NI,NJ,NK,NL,NMAX,MMAX,
     +                BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00,
     +                DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL
      COMMON /xSHLINF/ AG(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH),
     *                CFA(MXGSH),CGA(MXGSH),
     *                BG(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH),
     *                CFB(MXGSH),CGB(MXGSH),
     *                CG(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH),
     *                CFC(MXGSH),CGC(MXGSH),
     *                DG(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH),
     *                CFD(MXGSH),CGD(MXGSH),
     *                XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK,
     *                NGA,NGB,NGC,NGD
      COMMON /xSHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     +                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     +                NIJ,IJ,KL,IJKL
      COMMON /xSHLT  / TOL,CUTOFF,OUT,ICOUNT
C
      DIMENSION IN1(9)
C
      PARAMETER (SQRT3=1.73205080756888D+00, SQRT5=2.23606797749979D+00,
     *           SQRT7=2.64575131106459D+00, PI252=34.986836655250D+00,
     *           ZERO=0.0D+00, HALF=0.5D+00, ONE=1.0D+00)
C
C     GENERAL INTEGRAL ROUTINE FOR SPD FUNCTIONS
C
      FACTOR = PI252*QQ4
      NORM = NORMF .NE. 1 .OR. NORMP .NE. 1
      NI = LIT-1
      NJ = LJT-1
      NK = LKT-1
      NL = LLT-1
      DXIJ = XI-XJ
      DYIJ = YI-YJ
      DZIJ = ZI-ZJ
      DXKL = XK-XL
      DYKL = YK-YL
      DZKL = ZK-ZL
      NMAX = NI+NJ
      MMAX = NK+NL
      MAX = NMAX+1
      DO 100 I = 1,MAX
         N = I-1
         IF (N .LE. NI) IN1(I) = 125*N+1
         IF (N .GT. NI) IN1(I) = 125*NI+25*(N-NI)+1
  100 CONTINUE
      MAX = MMAX+1
      DO 120 K = 1,MAX
         N = K-1
         IF (N .LE. NK) KN(K) = 5*N
         IF (N .GT. NK) KN(K) = 5*NK+N-NK
  120 CONTINUE
C
C     ----- K PRIMITIVE
C
      LGMAX = NGD
      DO 480 KG = 1,NGC
         AK = CG(KG)
         BRRK = AK*RRK
         AKXK = AK*XK
         AKYK = AK*YK
         AKZK = AK*ZK
         CSK = CSC(KG)*FACTOR
         CPK = CPC(KG)*FACTOR
         CDK = CDC(KG)*FACTOR
         CFK = CFC(KG)*FACTOR
         CGK = CGC(KG)*FACTOR
C
C        ----- L PRIMITIVE
C
         IF (KANDL) LGMAX = KG
         DO 460 LG = 1,LGMAX
            AL = DG(LG)
            B = AK+AL
            BINV = ONE/B
            BBRRK = AL*BRRK*BINV
            IF (BBRRK .GT. TOL) GO TO 460
            CSL = CSD(LG)
            CPL = CPD(LG)
            CDL = CDD(LG)
            CFL = CFD(LG)
            CGL = CGD(LG)
            XB = (AKXK+AL*XL)*BINV
            YB = (AKYK+AL*YL)*BINV
            ZB = (AKZK+AL*ZL)*BINV
            BXBK = B*(XB-XK)
            BYBK = B*(YB-YK)
            BZBK = B*(ZB-ZK)
            BXBI = B*(XB-XI)
            BYBI = B*(YB-YI)
            BZBI = B*(ZB-ZI)
C
C           ----- DENSITY FACTOR
C
            DOUBLE=KANDL.AND.KG.NE.LG
            N = 0
            MAX = MAXL
            DUM1 = ZERO
            DUM2 = ZERO
            DO 370 K = MINK,MAXK
               GO TO (140,160,220,220,180,220,220,200,220,220,
     1                201,220,220,202,220,220,220,220,220,203,
     1                204,220,220,205,220,220,220,220,220,206,
     1                220,220,207,220,220),K
  140          DUM1 = CSK*BINV
               GO TO 220
  160          DUM1 = CPK*BINV
               GO TO 220
  180          DUM1 = CDK*BINV
               GO TO 220
  200          IF (NORM) DUM1 = DUM1*SQRT3
               GO TO 220
  201          DUM1 = CFK*BINV
               GO TO 220
  202          IF (NORM) DUM1 = DUM1*SQRT5
               GO TO 220
  203          IF (NORM) DUM1 = DUM1*SQRT3
               GO TO 220
  204          DUM1 = CGK*BINV
               GO TO 220
  205          IF (NORM) DUM1 = DUM1*SQRT7
               GO TO 220
  206          IF (NORM) DUM1 = DUM1*SQRT5/SQRT3
               GO TO 220
  207          IF (NORM) DUM1 = DUM1*SQRT3
  220          IF (KANDL) MAX = K
               DO 360 L = MINL,MAX
                  GO TO (240,280,340,340,300,340,340,320,340,340,
     1                   321,340,340,322,340,340,340,340,340,323,
     1                   324,340,340,325,340,340,340,340,340,326,
     1                   340,340,327,340,340),L
  240             DUM2 = DUM1*CSL
                  IF ( .NOT. DOUBLE) GO TO 340
                  IF (K .GT. 1) GO TO 260
                  DUM2 = DUM2+DUM2
                  GO TO 340
  260             DUM2 = DUM2+CSK*CPL*BINV
                  GO TO 340
  280             DUM2 = DUM1*CPL
                  IF (DOUBLE) DUM2 = DUM2+DUM2
                  GO TO 340
  300             DUM2 = DUM1*CDL
                  IF (DOUBLE) DUM2 = DUM2+DUM2
                  GO TO 340
  320             IF (NORM) DUM2 = DUM2*SQRT3
                  GO TO 340
  321             DUM2 = DUM1*CFL
                  IF (DOUBLE) DUM2 = DUM2+DUM2
                  GO TO 340
  322             IF (NORM) DUM2 = DUM2*SQRT5
                  GO TO 340
  323             IF (NORM) DUM2 = DUM2*SQRT3
                  GO TO 340
  324             DUM2 = DUM1*CGL
                  IF (DOUBLE) DUM2 = DUM2+DUM2
                  GO TO 340
  325             IF (NORM) DUM2 = DUM2*SQRT7
                  GO TO 340
  326             IF (NORM) DUM2 = DUM2*SQRT5/SQRT3
                  GO TO 340
  327             IF (NORM) DUM2 = DUM2*SQRT3
  340             N = N+1
                  DKL(N) = DUM2
  360          CONTINUE
  370       CONTINUE
C
C           ----- PAIR OF I,J PRIMITIVES
C
            NN = 0
            DO 440 N = 1,NIJ
               DUM = BBRRK+R(N)
               IF (DUM .GT. TOL) GO TO 440
               DO 380 I = 1,IJ
                  DIJ(I) = DDIJ(IJD(I)+NN)
  380          CONTINUE
               A = AA(N)
               AB = A*B
               AANDB = A+B
               EXPE = EXP(-DUM)/SQRT(AANDB)
               RHO = AB/AANDB
               XA = X1(N)
               YA = Y1(N)
               ZA = Z1(N)
               XX = RHO*((XA-XB)*(XA-XB) + (YA-YB)*(YA-YB)
     *                                   + (ZA-ZB)*(ZA-ZB))
               AXAK = A*(XA-XK)
               AYAK = A*(YA-YK)
               AZAK = A*(ZA-ZK)
               AXAI = A*(XA-XI)
               AYAI = A*(YA-YI)
               AZAI = A*(ZA-ZI)
               C1X = BXBK+AXAK
               C2X = A*BXBK
               C3X = BXBI+AXAI
               C4X = B*AXAI
               C1Y = BYBK+AYAK
               C2Y = A*BYBK
               C3Y = BYBI+AYAI
               C4Y = B*AYAI
               C1Z = BZBK+AZAK
               C2Z = A*BZBK
               C3Z = BZBI+AZAI
               C4Z = B*AZAI
C
C              ----- ROOTS AND WEIGHTS FOR QUADRATURE
C
               IF (NROOTS .LE. 3) CALL xRT123
               IF (NROOTS .EQ. 4) CALL xROOT4
               IF (NROOTS .EQ. 5) CALL xROOT5
               IF (NROOTS .GE. 6) CALL xROOT6
               MM = 0
               MAX = NMAX+1
C
C              COMPUTE TWO-ELECTRON INTEGRALS FOR EACH ROOT
C
               DO 420 M = 1,NROOTS
                  U2 = U(M)*RHO
                  F00 = EXPE*W(M)
                  DO 400 I = 1,MAX
                     IN(I) = IN1(I)+MM
  400             CONTINUE
                  DUMINV = ONE/(AB+U2*AANDB)
                  DM2INV = HALF*DUMINV
                  BP01 = (A+U2)*DM2INV
                  B00 = U2*DM2INV
                  B10 = (B+U2)*DM2INV
                  XCP00 = (U2*C1X+C2X)*DUMINV
                  XC00 = (U2*C3X+C4X)*DUMINV
                  YCP00 = (U2*C1Y+C2Y)*DUMINV
                  YC00 = (U2*C3Y+C4Y)*DUMINV
                  ZCP00 = (U2*C1Z+C2Z)*DUMINV
                  ZC00 = (U2*C3Z+C4Z)*DUMINV
                  CALL xXYZINT
                  MM = MM+625
  420          CONTINUE
C
C              ----- FORM (I,J//K,L) INTEGRALS OVER FUNCTIONS
C
               CALL xFORMS(GHONDO)
  440       NN = NN+16
  460    CONTINUE
  480 CONTINUE
C
      RETURN
      END
c-------------------------------------------------------------------
C*MODULE INT2A   *DECK BASCHK
      SUBROUTINE xBASCHK(LMAX)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
#include "gamess.par"
C
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
C
C     RETURN THE HIGHEST ANGULAR MOMENTUM PRESENT IN THE BASIS.
C     NOTE THAT KTYPE=1,2,3,4,5 MEANS S, P(L), D, F, G FUNCTION.
C
      KANG = 0
      DO 100 N=1,NSHELL
          IF(KTYPE(N).GT.KANG) KANG = KTYPE(N)
  100 CONTINUE
      LMAX = KANG-1
      RETURN
      END
c-------------------------------------------------------------------
C*MODULE RYSPOL  *DECK ROOT4
      SUBROUTINE xROOT4
CDIR$ INLINEALWAYS XROOT4
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C          *****   VERSION FEBRUARY 16,1975   *****
      COMMON /xROOT  / X,U(13),W(13),NROOTS
C
      EQUIVALENCE (U(1),RT1),(U(2),RT2),(U(3),RT3),(U(4),RT4),(U(5),RT5)
      EQUIVALENCE (W(1),WW1),(W(2),WW2),(W(3),WW3),(W(4),WW4),(W(5),WW5)
C
      DATA R14,PIE4/1.45303521503316D-01, 7.85398163397448D-01/
      DATA R24,W24/ 1.33909728812636D+00, 2.34479815323517D-01/
      DATA R34,W34/ 3.92696350135829D+00, 1.92704402415764D-02/
      DATA R44,W44/ 8.58863568901199D+00, 2.25229076750736D-04/
C
      IF (X .GT. 15.0D+00) GO TO 180
      IF (X .GT. 5.0D+00) GO TO 140
      IF (X .GT. 1.0D+00) GO TO 120
      IF (X .GT. 3.0D-07) GO TO 100
C     X IS APPROXIMATELY ZERO.                   NROOTS = 4
      RT1 = 3.48198973061471D-02 -4.09645850660395D-03 *X
      RT2 = 3.81567185080042D-01 -4.48902570656719D-02 *X
      RT3 = 1.73730726945891D+00 -2.04389090547327D-01 *X
      RT4 = 1.18463056481549D+01 -1.39368301742312D+00 *X
      WW1 = 3.62683783378362D-01 -3.13844305713928D-02 *X
      WW2 = 3.13706645877886D-01 -8.98046242557724D-02 *X
      WW3 = 2.22381034453372D-01 -1.29314370958973D-01 *X
      WW4 = 1.01228536290376D-01 -8.28299075414321D-02 *X
      RETURN
C
C     X=0.0 TO 1.0                               NROOTS = 4
  100 RT1 = ((((((-1.95309614628539D-10*X+5.19765728707592D-09)*X-
     +     1.01756452250573D-07 )*X+1.72365935872131D-06 )*X-
     +     2.61203523522184D-05 )*X+3.52921308769880D-04 )*X-
     +     4.09645850658433D-03 )*X+3.48198973061469D-02
      RT2 = (((((-1.89554881382342D-08*X+3.07583114342365D-07)*X+
     +     1.270981734393D-06)*X-1.417298563884D-04)*X+
     +     3.226979163176D-03)*X-4.48902570678178D-02 )*X+
     +     3.81567185080039D-01
      RT3 = (((((( 1.77280535300416D-09*X+3.36524958870615D-08)*X-
     +     2.58341529013893D-07 )*X-1.13644895662320D-05 )*X-
     +     7.91549618884063D-05 )*X+1.03825827346828D-02 )*X-
     +     2.04389090525137D-01 )*X+1.73730726945889D+00
      RT4 = (((((-5.61188882415248D-08*X-2.49480733072460D-07)*X+
     +     3.428685057114D-06)*X+1.679007454539D-04)*X+
     +     4.722855585715D-02)*X-1.39368301737828D+00 )*X+
     +     1.18463056481543D+01
      WW1 = ((((((-1.14649303201279D-08*X+1.88015570196787D-07)*X-
     +     2.33305875372323D-06 )*X+2.68880044371597D-05 )*X-
     +     2.94268428977387D-04 )*X+3.06548909776613D-03 )*X-
     +     3.13844305680096D-02 )*X+3.62683783378335D-01
      WW2 = ((((((((-4.11720483772634D-09*X+6.54963481852134D-08)*X-
     +     7.20045285129626D-07 )*X+6.93779646721723D-06 )*X-
     +     6.05367572016373D-05 )*X+4.74241566251899D-04 )*X-
     +     3.26956188125316D-03 )*X+1.91883866626681D-02 )*X-
     +     8.98046242565811D-02 )*X+3.13706645877886D-01
      WW3 = ((((((((-3.41688436990215D-08*X+5.07238960340773D-07)*X-
     +     5.01675628408220D-06 )*X+4.20363420922845D-05 )*X-
     +     3.08040221166823D-04 )*X+1.94431864731239D-03 )*X-
     +     1.02477820460278D-02 )*X+4.28670143840073D-02 )*X-
     +     1.29314370962569D-01 )*X+2.22381034453369D-01
      WW4 = ((((((((( 4.99660550769508D-09*X-7.94585963310120D-08)*X+
     +     8.359072409485D-07)*X-7.422369210610D-06)*X+
     +     5.763374308160D-05)*X-3.86645606718233D-04 )*X+
     +     2.18417516259781D-03 )*X-9.99791027771119D-03 )*X+
     +     3.48791097377370D-02 )*X-8.28299075413889D-02 )*X+
     +     1.01228536290376D-01
      RETURN
C
C     X= 1.0 TO 5.0                              NROOTS = 4
  120 Y = X-3.0D+00
      RT1 = (((((((((-1.48570633747284D-15*Y-1.33273068108777D-13)*Y+
     +     4.068543696670D-12)*Y-9.163164161821D-11)*Y+
     +     2.046819017845D-09)*Y-4.03076426299031D-08 )*Y+
     +     7.29407420660149D-07 )*Y-1.23118059980833D-05 )*Y+
     +     1.88796581246938D-04 )*Y-2.53262912046853D-03 )*Y+
     +     2.51198234505021D-02
      RT2 = ((((((((( 1.35830583483312D-13*Y-2.29772605964836D-12)*Y-
     +     3.821500128045D-12)*Y+6.844424214735D-10)*Y-
     +     1.048063352259D-08)*Y+1.50083186233363D-08 )*Y+
     +     3.48848942324454D-06 )*Y-1.08694174399193D-04 )*Y+
     +     2.08048885251999D-03 )*Y-2.91205805373793D-02 )*Y+
     +     2.72276489515713D-01
      RT3 = ((((((((( 5.02799392850289D-13*Y+1.07461812944084D-11)*Y-
     +     1.482277886411D-10)*Y-2.153585661215D-09)*Y+
     +     3.654087802817D-08)*Y+5.15929575830120D-07 )*Y-
     +     9.52388379435709D-06 )*Y-2.16552440036426D-04 )*Y+
     +     9.03551469568320D-03 )*Y-1.45505469175613D-01 )*Y+
     +     1.21449092319186D+00
      RT4 = (((((((((-1.08510370291979D-12*Y+6.41492397277798D-11)*Y+
     +     7.542387436125D-10)*Y-2.213111836647D-09)*Y-
     +     1.448228963549D-07)*Y-1.95670833237101D-06 )*Y-
     +     1.07481314670844D-05 )*Y+1.49335941252765D-04 )*Y+
     +     4.87791531990593D-02 )*Y-1.10559909038653D+00 )*Y+
     +     8.09502028611780D+00
      WW1 = ((((((((((-4.65801912689961D-14*Y+7.58669507106800D-13)*Y-
     +     1.186387548048D-11)*Y+1.862334710665D-10)*Y-
     +     2.799399389539D-09)*Y+4.148972684255D-08)*Y-
     +     5.933568079600D-07)*Y+8.168349266115D-06)*Y-
     +     1.08989176177409D-04 )*Y+1.41357961729531D-03 )*Y-
     +     1.87588361833659D-02 )*Y+2.89898651436026D-01
      WW2 = ((((((((((((-1.46345073267549D-14*Y+2.25644205432182D-13)*Y-
     +     3.116258693847D-12)*Y+4.321908756610D-11)*Y-
     +     5.673270062669D-10)*Y+7.006295962960D-09)*Y-
     +     8.120186517000D-08)*Y+8.775294645770D-07)*Y-
     +     8.77829235749024D-06 )*Y+8.04372147732379D-05 )*Y-
     +     6.64149238804153D-04 )*Y+4.81181506827225D-03 )*Y-
     +     2.88982669486183D-02 )*Y+1.56247249979288D-01
      WW3 = ((((((((((((( 9.06812118895365D-15*Y-1.40541322766087D-13)*
     +     Y+1.919270015269D-12)*Y-2.605135739010D-11)*Y+
     +     3.299685839012D-10)*Y-3.86354139348735D-09 )*Y+
     +     4.16265847927498D-08 )*Y-4.09462835471470D-07 )*Y+
     +     3.64018881086111D-06 )*Y-2.88665153269386D-05 )*Y+
     +     2.00515819789028D-04 )*Y-1.18791896897934D-03 )*Y+
     +     5.75223633388589D-03 )*Y-2.09400418772687D-02 )*Y+
     +     4.85368861938873D-02
      WW4 = ((((((((((((((-9.74835552342257D-16*Y+1.57857099317175D-14)*
     +     Y-2.249993780112D-13)*Y+3.173422008953D-12)*Y-
     +     4.161159459680D-11)*Y+5.021343560166D-10)*Y-
     +     5.545047534808D-09)*Y+5.554146993491D-08)*Y-
     +     4.99048696190133D-07 )*Y+3.96650392371311D-06 )*Y-
     +     2.73816413291214D-05 )*Y+1.60106988333186D-04 )*Y-
     +     7.64560567879592D-04 )*Y+2.81330044426892D-03 )*Y-
     +     7.16227030134947D-03 )*Y+9.66077262223353D-03
      RETURN
C
  140 IF (X .GT. 10.0D+00) GO TO 160
C     X=5.0 TO 10.0                              NROOTS = 4
      Y = X-7.5D+00
      RT1 = ((((((((( 4.64217329776215D-15*Y-6.27892383644164D-15)*Y+
     +     3.462236347446D-13)*Y-2.927229355350D-11)*Y+
     +     5.090355371676D-10)*Y-9.97272656345253D-09 )*Y+
     +     2.37835295639281D-07 )*Y-4.60301761310921D-06 )*Y+
     +     8.42824204233222D-05 )*Y-1.37983082233081D-03 )*Y+
     +     1.66630865869375D-02
      RT2 = ((((((((( 2.93981127919047D-14*Y+8.47635639065744D-13)*Y-
     +     1.446314544774D-11)*Y-6.149155555753D-12)*Y+
     +     8.484275604612D-10)*Y-6.10898827887652D-08 )*Y+
     +     2.39156093611106D-06 )*Y-5.35837089462592D-05 )*Y+
     +     1.00967602595557D-03 )*Y-1.57769317127372D-02 )*Y+
     +     1.74853819464285D-01
      RT3 = (((((((((( 2.93523563363000D-14*Y-6.40041776667020D-14)*Y-
     +     2.695740446312D-12)*Y+1.027082960169D-10)*Y-
     +     5.822038656780D-10)*Y-3.159991002539D-08)*Y+
     +     4.327249251331D-07)*Y+4.856768455119D-06)*Y-
     +     2.54617989427762D-04 )*Y+5.54843378106589D-03 )*Y-
     +     7.95013029486684D-02 )*Y+7.20206142703162D-01
      RT4 = (((((((((((-1.62212382394553D-14*Y+7.68943641360593D-13)*Y+
     +     5.764015756615D-12)*Y-1.380635298784D-10)*Y-
     +     1.476849808675D-09)*Y+1.84347052385605D-08 )*Y+
     +     3.34382940759405D-07 )*Y-1.39428366421645D-06 )*Y-
     +     7.50249313713996D-05 )*Y-6.26495899187507D-04 )*Y+
     +     4.69716410901162D-02 )*Y-6.66871297428209D-01 )*Y+
     +     4.11207530217806D+00
      WW1 = ((((((((((-1.65995045235997D-15*Y+6.91838935879598D-14)*Y-
     +     9.131223418888D-13)*Y+1.403341829454D-11)*Y-
     +     3.672235069444D-10)*Y+6.366962546990D-09)*Y-
     +     1.039220021671D-07)*Y+1.959098751715D-06)*Y-
     +     3.33474893152939D-05 )*Y+5.72164211151013D-04 )*Y-
     +     1.05583210553392D-02 )*Y+2.26696066029591D-01
      WW2 = ((((((((((((-3.57248951192047D-16*Y+6.25708409149331D-15)*Y-
     +     9.657033089714D-14)*Y+1.507864898748D-12)*Y-
     +     2.332522256110D-11)*Y+3.428545616603D-10)*Y-
     +     4.698730937661D-09)*Y+6.219977635130D-08)*Y-
     +     7.83008889613661D-07 )*Y+9.08621687041567D-06 )*Y-
     +     9.86368311253873D-05 )*Y+9.69632496710088D-04 )*Y-
     +     8.14594214284187D-03 )*Y+8.50218447733457D-02
      WW3 = ((((((((((((( 1.64742458534277D-16*Y-2.68512265928410D-15)*
     +     Y+3.788890667676D-14)*Y-5.508918529823D-13)*Y+
     +     7.555896810069D-12)*Y-9.69039768312637D-11 )*Y+
     +     1.16034263529672D-09 )*Y-1.28771698573873D-08 )*Y+
     +     1.31949431805798D-07 )*Y-1.23673915616005D-06 )*Y+
     +     1.04189803544936D-05 )*Y-7.79566003744742D-05 )*Y+
     +     5.03162624754434D-04 )*Y-2.55138844587555D-03 )*Y+
     +     1.13250730954014D-02
      WW4 = ((((((((((((((-1.55714130075679D-17*Y+2.57193722698891D-16)*
     +     Y-3.626606654097D-15)*Y+5.234734676175D-14)*Y-
     +     7.067105402134D-13)*Y+8.793512664890D-12)*Y-
     +     1.006088923498D-10)*Y+1.050565098393D-09)*Y-
     +     9.91517881772662D-09 )*Y+8.35835975882941D-08 )*Y-
     +     6.19785782240693D-07 )*Y+3.95841149373135D-06 )*Y-
     +     2.11366761402403D-05 )*Y+9.00474771229507D-05 )*Y-
     +     2.78777909813289D-04 )*Y+5.26543779837487D-04
      RETURN
C
C     X=10.0 TO 15.0                             NROOTS = 4
  160 Y = X-12.5D+00
      RT1 = ((((((((((( 4.94869622744119D-17*Y+8.03568805739160D-16)*Y-
     +     5.599125915431D-15)*Y-1.378685560217D-13)*Y+
     +     7.006511663249D-13)*Y+1.30391406991118D-11 )*Y+
     +     8.06987313467541D-11 )*Y-5.20644072732933D-09 )*Y+
     +     7.72794187755457D-08 )*Y-1.61512612564194D-06 )*Y+
     +     4.15083811185831D-05 )*Y-7.87855975560199D-04 )*Y+
     +     1.14189319050009D-02
      RT2 = ((((((((((( 4.89224285522336D-16*Y+1.06390248099712D-14)*Y-
     +     5.446260182933D-14)*Y-1.613630106295D-12)*Y+
     +     3.910179118937D-12)*Y+1.90712434258806D-10 )*Y+
     +     8.78470199094761D-10 )*Y-5.97332993206797D-08 )*Y+
     +     9.25750831481589D-07 )*Y-2.02362185197088D-05 )*Y+
     +     4.92341968336776D-04 )*Y-8.68438439874703D-03 )*Y+
     +     1.15825965127958D-01
      RT3 = (((((((((( 6.12419396208408D-14*Y+1.12328861406073D-13)*Y-
     +     9.051094103059D-12)*Y-4.781797525341D-11)*Y+
     +     1.660828868694D-09)*Y+4.499058798868D-10)*Y-
     +     2.519549641933D-07)*Y+4.977444040180D-06)*Y-
     +     1.25858350034589D-04 )*Y+2.70279176970044D-03 )*Y-
     +     3.99327850801083D-02 )*Y+4.33467200855434D-01
      RT4 = ((((((((((( 4.63414725924048D-14*Y-4.72757262693062D-14)*Y-
     +     1.001926833832D-11)*Y+6.074107718414D-11)*Y+
     +     1.576976911942D-09)*Y-2.01186401974027D-08 )*Y-
     +     1.84530195217118D-07 )*Y+5.02333087806827D-06 )*Y+
     +     9.66961790843006D-06 )*Y-1.58522208889528D-03 )*Y+
     +     2.80539673938339D-02 )*Y-2.78953904330072D-01 )*Y+
     +     1.82835655238235D+00
      WW4 = ((((((((((((( 2.90401781000996D-18*Y-4.63389683098251D-17)*
     +     Y+6.274018198326D-16)*Y-8.936002188168D-15)*Y+
     +     1.194719074934D-13)*Y-1.45501321259466D-12 )*Y+
     +     1.64090830181013D-11 )*Y-1.71987745310181D-10 )*Y+
     +     1.63738403295718D-09 )*Y-1.39237504892842D-08 )*Y+
     +     1.06527318142151D-07 )*Y-7.27634957230524D-07 )*Y+
     +     4.12159381310339D-06 )*Y-1.74648169719173D-05 )*Y+
     +     8.50290130067818D-05
      WW3 = ((((((((((((-4.19569145459480D-17*Y+5.94344180261644D-16)*Y-
     +     1.148797566469D-14)*Y+1.881303962576D-13)*Y-
     +     2.413554618391D-12)*Y+3.372127423047D-11)*Y-
     +     4.933988617784D-10)*Y+6.116545396281D-09)*Y-
     +     6.69965691739299D-08 )*Y+7.52380085447161D-07 )*Y-
     +     8.08708393262321D-06 )*Y+6.88603417296672D-05 )*Y-
     +     4.67067112993427D-04 )*Y+5.42313365864597D-03
      WW2 = ((((((((((-6.22272689880615D-15*Y+1.04126809657554D-13)*Y-
     +     6.842418230913D-13)*Y+1.576841731919D-11)*Y-
     +     4.203948834175D-10)*Y+6.287255934781D-09)*Y-
     +     8.307159819228D-08)*Y+1.356478091922D-06)*Y-
     +     2.08065576105639D-05 )*Y+2.52396730332340D-04 )*Y-
     +     2.94484050194539D-03 )*Y+6.01396183129168D-02
      WW1 = (((-1.8784686463512D-01/X+2.2991849164985D-01)/X -
     +     4.9893752514047D-01)/X-2.1916512131607D-05)*EXP(-X) +
     +     SQRT(PIE4/X)-WW4-WW3-WW2
      RETURN
C
  180 WW1 = SQRT(PIE4/X)
      IF (X .GT. 35.0D+00) GO TO 220
      IF (X .GT. 20.0D+00) GO TO 200
C     X=15.0 TO 20.0                             NROOTS = 4
      Y = X-17.5D+00
      RT1 = ((((((((((( 4.36701759531398D-17*Y-1.12860600219889D-16)*Y-
     +     6.149849164164D-15)*Y+5.820231579541D-14)*Y+
     +     4.396602872143D-13)*Y-1.24330365320172D-11 )*Y+
     +     6.71083474044549D-11 )*Y+2.43865205376067D-10 )*Y+
     +     1.67559587099969D-08 )*Y-9.32738632357572D-07 )*Y+
     +     2.39030487004977D-05 )*Y-4.68648206591515D-04 )*Y+
     +     8.34977776583956D-03
      RT2 = ((((((((((( 4.98913142288158D-16*Y-2.60732537093612D-16)*Y-
     +     7.775156445127D-14)*Y+5.766105220086D-13)*Y+
     +     6.432696729600D-12)*Y-1.39571683725792D-10 )*Y+
     +     5.95451479522191D-10 )*Y+2.42471442836205D-09 )*Y+
     +     2.47485710143120D-07 )*Y-1.14710398652091D-05 )*Y+
     +     2.71252453754519D-04 )*Y-4.96812745851408D-03 )*Y+
     +     8.26020602026780D-02
      RT3 = ((((((((((( 1.91498302509009D-15*Y+1.48840394311115D-14)*Y-
     +     4.316925145767D-13)*Y+1.186495793471D-12)*Y+
     +     4.615806713055D-11)*Y-5.54336148667141D-10 )*Y+
     +     3.48789978951367D-10 )*Y-2.79188977451042D-09 )*Y+
     +     2.09563208958551D-06 )*Y-6.76512715080324D-05 )*Y+
     +     1.32129867629062D-03 )*Y-2.05062147771513D-02 )*Y+
     +     2.88068671894324D-01
      RT4 = (((((((((((-5.43697691672942D-15*Y-1.12483395714468D-13)*Y+
     +     2.826607936174D-12)*Y-1.266734493280D-11)*Y-
     +     4.258722866437D-10)*Y+9.45486578503261D-09 )*Y-
     +     5.86635622821309D-08 )*Y-1.28835028104639D-06 )*Y+
     +     4.41413815691885D-05 )*Y-7.61738385590776D-04 )*Y+
     +     9.66090902985550D-03 )*Y-1.01410568057649D-01 )*Y+
     +     9.54714798156712D-01
      WW4 = ((((((((((((-7.56882223582704D-19*Y+7.53541779268175D-18)*Y-
     +     1.157318032236D-16)*Y+2.411195002314D-15)*Y-
     +     3.601794386996D-14)*Y+4.082150659615D-13)*Y-
     +     4.289542980767D-12)*Y+5.086829642731D-11)*Y-
     +     6.35435561050807D-10 )*Y+6.82309323251123D-09 )*Y-
     +     5.63374555753167D-08 )*Y+3.57005361100431D-07 )*Y-
     +     2.40050045173721D-06 )*Y+4.94171300536397D-05
      WW3 = (((((((((((-5.54451040921657D-17*Y+2.68748367250999D-16)*Y+
     +     1.349020069254D-14)*Y-2.507452792892D-13)*Y+
     +     1.944339743818D-12)*Y-1.29816917658823D-11 )*Y+
     +     3.49977768819641D-10 )*Y-8.67270669346398D-09 )*Y+
     +     1.31381116840118D-07 )*Y-1.36790720600822D-06 )*Y+
     +     1.19210697673160D-05 )*Y-1.42181943986587D-04 )*Y+
     +     4.12615396191829D-03
      WW2 = (((((((((((-1.86506057729700D-16*Y+1.16661114435809D-15)*Y+
     +     2.563712856363D-14)*Y-4.498350984631D-13)*Y+
     +     1.765194089338D-12)*Y+9.04483676345625D-12 )*Y+
     +     4.98930345609785D-10 )*Y-2.11964170928181D-08 )*Y+
     +     3.98295476005614D-07 )*Y-5.49390160829409D-06 )*Y+
     +     7.74065155353262D-05 )*Y-1.48201933009105D-03 )*Y+
     +     4.97836392625268D-02
      WW1 = (( 1.9623264149430D-01/X-4.9695241464490D-01)/X -
     +     6.0156581186481D-05)*EXP(-X)+WW1-WW2-WW3-WW4
      RETURN
C
C     X=20.0 TO 35.0                             NROOTS = 4
  200 E = EXP(-X)
      RT1 = ((((((-4.45711399441838D-05*X+1.27267770241379D-03)*X -
     +     2.36954961381262D-01)*X+1.54330657903756D+01)*X -
     +     5.22799159267808D+02)*X+1.05951216669313D+04)*X + (-
     +     2.51177235556236D+06/X+8.72975373557709D+05)/X -
     +     1.29194382386499D+05)*E + R14/(X-R14)
      RT2 = (((((-7.85617372254488D-02*X+6.35653573484868D+00)*X -
     +     3.38296938763990D+02)*X+1.25120495802096D+04)*X -
     +     3.16847570511637D+05)*X + ((-1.02427466127427D+09/X +
     +     3.70104713293016D+08)/X-5.87119005093822D+07)/X +
     +     5.38614211391604D+06)*E + R24/(X-R24)
      RT3 = (((((-2.37900485051067D-01*X+1.84122184400896D+01)*X -
     +     1.00200731304146D+03)*X+3.75151841595736D+04)*X -
     +     9.50626663390130D+05)*X + ((-2.88139014651985D+09/X +
     +     1.06625915044526D+09)/X-1.72465289687396D+08)/X +
     +     1.60419390230055D+07)*E + R34/(X-R34)
      RT4 = ((((((-6.00691586407385D-04*X-3.64479545338439D-01)*X +
     +     1.57496131755179D+01)*X-6.54944248734901D+02)*X +
     +     1.70830039597097D+04)*X-2.90517939780207D+05)*X + (+
     +     3.49059698304732D+07/X-1.64944522586065D+07)/X +
     +     2.96817940164703D+06)*E + R44/(X-R44)
      IF (X .LE. 25.0D+00) WW4 = ((((((( 2.33766206773151D-07*X-
     +     3.81542906607063D-05)*X +3.51416601267000D-03)*X-
     +     1.66538571864728D-01)*X +4.80006136831847D+00)*X-
     +     8.73165934223603D+01)*X +9.77683627474638D+02)*X +
     +     1.66000945117640D+04/X -6.14479071209961D+03)*E + W44*WW1
      IF (X .GT. 25.0D+00) WW4 = (((((( 5.74245945342286D-06*X-
     +     7.58735928102351D-05)*X +2.35072857922892D-04)*X-
     +     3.78812134013125D-03)*X +3.09871652785805D-01)*X-
     +     7.11108633061306D+00)*X +5.55297573149528D+01)*E + W44*WW1
      WW3 = (((((( 2.36392855180768D-04*X-9.16785337967013D-03)*X +
     +     4.62186525041313D-01)*X-1.96943786006540D+01)*X +
     +     4.99169195295559D+02)*X-6.21419845845090D+03)*X + ((+
     +     5.21445053212414D+07/X-1.34113464389309D+07)/X +
     +     1.13673298305631D+06)/X-2.81501182042707D+03)*E + W34*WW1
      WW2 = (((((( 7.29841848989391D-04*X-3.53899555749875D-02)*X +
     +     2.07797425718513D+00)*X-1.00464709786287D+02)*X +
     +     3.15206108877819D+03)*X-6.27054715090012D+04)*X + (+
     +     1.54721246264919D+07/X-5.26074391316381D+06)/X +
     +     7.67135400969617D+05)*E + W24*WW1
      WW1 = (( 1.9623264149430D-01/X-4.9695241464490D-01)/X -
     +     6.0156581186481D-05)*E + WW1-WW2-WW3-WW4
      RETURN
C
  220 IF (X .GT. 53.0D+00) GO TO 240
C     X=35.0 TO 53.0                             NROOTS = 4
      E = EXP(-X)*(X*X)**2
      RT4 = ((-2.19135070169653D-03*X-1.19108256987623D-01)*X -
     +     7.50238795695573D-01)*E + R44/(X-R44)
      RT3 = ((-9.65842534508637D-04*X-4.49822013469279D-02)*X +
     +     6.08784033347757D-01)*E + R34/(X-R34)
      RT2 = ((-3.62569791162153D-04*X-9.09231717268466D-03)*X +
     +     1.84336760556262D-01)*E + R24/(X-R24)
      RT1 = ((-4.07557525914600D-05*X-6.88846864931685D-04)*X +
     +     1.74725309199384D-02)*E + R14/(X-R14)
      WW4 = (( 5.76631982000990D-06*X-7.89187283804890D-05)*X +
     +     3.28297971853126D-04)*E + W44*WW1
      WW3 = (( 2.08294969857230D-04*X-3.77489954837361D-03)*X +
     +     2.09857151617436D-02)*E + W34*WW1
      WW2 = (( 6.16374517326469D-04*X-1.26711744680092D-02)*X +
     +     8.14504890732155D-02)*E + W24*WW1
      WW1 = WW1-WW2-WW3-WW4
      RETURN
C
C     X=47.0 TO INFINITY                         NROOTS = 4
  240 RT1 = R14/(X-R14)
      RT2 = R24/(X-R24)
      RT3 = R34/(X-R34)
      RT4 = R44/(X-R44)
      WW4 = W44*WW1
      WW3 = W34*WW1
      WW2 = W24*WW1
      WW1 = WW1-WW2-WW3-WW4
      RETURN
      END
c-------------------------------------------------------------------
C*MODULE RYSPOL  *DECK ROOT5
      SUBROUTINE xROOT5
CDIR$ INLINEALWAYS XROOT5
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C          *****   VERSION  FEBRUARY 27,1975   *****
C
      COMMON /xROOT  / X,U(13),W(13),NROOTS
C
      EQUIVALENCE (U(1),RT1),(U(2),RT2),(U(3),RT3),(U(4),RT4),(U(5),RT5)
      EQUIVALENCE (W(1),WW1),(W(2),WW2),(W(3),WW3),(W(4),WW4),(W(5),WW5)
C
      DATA R15,PIE4/1.17581320211778D-01, 7.85398163397448D-01/
      DATA R25,W25/ 1.07456201243690D+00, 2.70967405960535D-01/
      DATA R35,W35/ 3.08593744371754D+00, 3.82231610015404D-02/
      DATA R45,W45/ 6.41472973366203D+00, 1.51614186862443D-03/
      DATA R55,W55/ 1.18071894899717D+01, 8.62130526143657D-06/
C
      IF (X .GT. 15.0D+00) GO TO 180
      IF (X .GT. 5.0D+00) GO TO 140
      IF (X .GT. 1.0D+00) GO TO 120
      IF (X .GT. 3.0D-07) GO TO 100
C     X IS APPROXIMATELY ZERO.                   NROOTS = 5
      RT1 = 2.26659266316985D-02 -2.15865967920897D-03 *X
      RT2 = 2.31271692140903D-01 -2.20258754389745D-02 *X
      RT3 = 8.57346024118836D-01 -8.16520023025515D-02 *X
      RT4 = 2.97353038120346D+00 -2.83193369647137D-01 *X
      RT5 = 1.84151859759051D+01 -1.75382723579439D+00 *X
      WW1 = 2.95524224714752D-01 -1.96867576909777D-02 *X
      WW2 = 2.69266719309995D-01 -5.61737590184721D-02 *X
      WW3 = 2.19086362515981D-01 -9.71152726793658D-02 *X
      WW4 = 1.49451349150580D-01 -1.02979262193565D-01 *X
      WW5 = 6.66713443086877D-02 -5.73782817488315D-02 *X
      RETURN
C
C     X=0.0 TO 1.0                               NROOTS = 5
  100 RT1 = ((((((-4.46679165328413D-11*X+1.21879111988031D-09)*X-
     +     2.62975022612104D-08 )*X+5.15106194905897D-07 )*X-
     +     9.27933625824749D-06 )*X+1.51794097682482D-04 )*X-
     +     2.15865967920301D-03 )*X+2.26659266316985D-02
      RT2 = (((((( 1.93117331714174D-10*X-4.57267589660699D-09)*X+
     +     2.48339908218932D-08 )*X+1.50716729438474D-06 )*X-
     +     6.07268757707381D-05 )*X+1.37506939145643D-03 )*X-
     +     2.20258754419939D-02 )*X+2.31271692140905D-01
      RT3 = ((((( 4.84989776180094D-09*X+1.31538893944284D-07)*X-
     +     2.766753852879D-06)*X-7.651163510626D-05)*X+
     +     4.033058545972D-03)*X-8.16520022916145D-02 )*X+
     +     8.57346024118779D-01
      RT4 = ((((-2.48581772214623D-07*X-4.34482635782585D-06)*X-
     +     7.46018257987630D-07 )*X+1.01210776517279D-02 )*X-
     +     2.83193369640005D-01 )*X+2.97353038120345D+00
      RT5 = (((((-8.92432153868554D-09*X+1.77288899268988D-08)*X+
     +     3.040754680666D-06)*X+1.058229325071D-04)*X+
     +     4.596379534985D-02)*X-1.75382723579114D+00 )*X+
     +     1.84151859759049D+01
      WW1 = ((((((-2.03822632771791D-09*X+3.89110229133810D-08)*X-
     +     5.84914787904823D-07 )*X+8.30316168666696D-06 )*X-
     +     1.13218402310546D-04 )*X+1.49128888586790D-03 )*X-
     +     1.96867576904816D-02 )*X+2.95524224714749D-01
      WW2 = ((((((( 8.62848118397570D-09*X-1.38975551148989D-07)*X+
     +     1.602894068228D-06)*X-1.646364300836D-05)*X+
     +     1.538445806778D-04)*X-1.28848868034502D-03 )*X+
     +     9.38866933338584D-03 )*X-5.61737590178812D-02 )*X+
     +     2.69266719309991D-01
      WW3 = ((((((((-9.41953204205665D-09*X+1.47452251067755D-07)*X-
     +     1.57456991199322D-06 )*X+1.45098401798393D-05 )*X-
     +     1.18858834181513D-04 )*X+8.53697675984210D-04 )*X-
     +     5.22877807397165D-03 )*X+2.60854524809786D-02 )*X-
     +     9.71152726809059D-02 )*X+2.19086362515979D-01
      WW4 = ((((((((-3.84961617022042D-08*X+5.66595396544470D-07)*X-
     +     5.52351805403748D-06 )*X+4.53160377546073D-05 )*X-
     +     3.22542784865557D-04 )*X+1.95682017370967D-03 )*X-
     +     9.77232537679229D-03 )*X+3.79455945268632D-02 )*X-
     +     1.02979262192227D-01 )*X+1.49451349150573D-01
      WW5 = ((((((((( 4.09594812521430D-09*X-6.47097874264417D-08)*X+
     +     6.743541482689D-07)*X-5.917993920224D-06)*X+
     +     4.531969237381D-05)*X-2.99102856679638D-04 )*X+
     +     1.65695765202643D-03 )*X-7.40671222520653D-03 )*X+
     +     2.50889946832192D-02 )*X-5.73782817487958D-02 )*X+
     +     6.66713443086877D-02
      RETURN
C
C     X=1.0 TO 5.0                               NROOTS = 5
  120 Y = X-3.0D+00
      RT1 = ((((((((-2.58163897135138D-14*Y+8.14127461488273D-13)*Y-
     +     2.11414838976129D-11 )*Y+5.09822003260014D-10 )*Y-
     +     1.16002134438663D-08 )*Y+2.46810694414540D-07 )*Y-
     +     4.92556826124502D-06 )*Y+9.02580687971053D-05 )*Y-
     +     1.45190025120726D-03 )*Y+1.73416786387475D-02
      RT2 = ((((((((( 1.04525287289788D-14*Y+5.44611782010773D-14)*Y-
     +     4.831059411392D-12)*Y+1.136643908832D-10)*Y-
     +     1.104373076913D-09)*Y-2.35346740649916D-08 )*Y+
     +     1.43772622028764D-06 )*Y-4.23405023015273D-05 )*Y+
     +     9.12034574793379D-04 )*Y-1.52479441718739D-02 )*Y+
     +     1.76055265928744D-01
      RT3 = (((((((((-6.89693150857911D-14*Y+5.92064260918861D-13)*Y+
     +     1.847170956043D-11)*Y-3.390752744265D-10)*Y-
     +     2.995532064116D-09)*Y+1.57456141058535D-07 )*Y-
     +     3.95859409711346D-07 )*Y-9.58924580919747D-05 )*Y+
     +     3.23551502557785D-03 )*Y-5.97587007636479D-02 )*Y+
     +     6.46432853383057D-01
      RT4 = ((((((((-3.61293809667763D-12*Y-2.70803518291085D-11)*Y+
     +     8.83758848468769D-10 )*Y+1.59166632851267D-08 )*Y-
     +     1.32581997983422D-07 )*Y-7.60223407443995D-06 )*Y-
     +     7.41019244900952D-05 )*Y+9.81432631743423D-03 )*Y-
     +     2.23055570487771D-01 )*Y+2.21460798080643D+00
      RT5 = ((((((((( 7.12332088345321D-13*Y+3.16578501501894D-12)*Y-
     +     8.776668218053D-11)*Y-2.342817613343D-09)*Y-
     +     3.496962018025D-08)*Y-3.03172870136802D-07 )*Y+
     +     1.50511293969805D-06 )*Y+1.37704919387696D-04 )*Y+
     +     4.70723869619745D-02 )*Y-1.47486623003693D+00 )*Y+
     +     1.35704792175847D+01
      WW1 = ((((((((( 1.04348658616398D-13*Y-1.94147461891055D-12)*Y+
     +     3.485512360993D-11)*Y-6.277497362235D-10)*Y+
     +     1.100758247388D-08)*Y-1.88329804969573D-07 )*Y+
     +     3.12338120839468D-06 )*Y-5.04404167403568D-05 )*Y+
     +     8.00338056610995D-04 )*Y-1.30892406559521D-02 )*Y+
     +     2.47383140241103D-01
      WW2 = ((((((((((( 3.23496149760478D-14*Y-5.24314473469311D-13)*Y+
     +     7.743219385056D-12)*Y-1.146022750992D-10)*Y+
     +     1.615238462197D-09)*Y-2.15479017572233D-08 )*Y+
     +     2.70933462557631D-07 )*Y-3.18750295288531D-06 )*Y+
     +     3.47425221210099D-05 )*Y-3.45558237388223D-04 )*Y+
     +     3.05779768191621D-03 )*Y-2.29118251223003D-02 )*Y+
     +     1.59834227924213D-01
      WW3 = ((((((((((((-3.42790561802876D-14*Y+5.26475736681542D-13)*Y-
     +     7.184330797139D-12)*Y+9.763932908544D-11)*Y-
     +     1.244014559219D-09)*Y+1.472744068942D-08)*Y-
     +     1.611749975234D-07)*Y+1.616487851917D-06)*Y-
     +     1.46852359124154D-05 )*Y+1.18900349101069D-04 )*Y-
     +     8.37562373221756D-04 )*Y+4.93752683045845D-03 )*Y-
     +     2.25514728915673D-02 )*Y+6.95211812453929D-02
      WW4 = ((((((((((((( 1.04072340345039D-14*Y-1.60808044529211D-13)*
     +     Y+2.183534866798D-12)*Y-2.939403008391D-11)*Y+
     +     3.679254029085D-10)*Y-4.23775673047899D-09 )*Y+
     +     4.46559231067006D-08 )*Y-4.26488836563267D-07 )*Y+
     +     3.64721335274973D-06 )*Y-2.74868382777722D-05 )*Y+
     +     1.78586118867488D-04 )*Y-9.68428981886534D-04 )*Y+
     +     4.16002324339929D-03 )*Y-1.28290192663141D-02 )*Y+
     +     2.22353727685016D-02
      WW5 = ((((((((((((((-8.16770412525963D-16*Y+1.31376515047977D-14)*
     +     Y-1.856950818865D-13)*Y+2.596836515749D-12)*Y-
     +     3.372639523006D-11)*Y+4.025371849467D-10)*Y-
     +     4.389453269417D-09)*Y+4.332753856271D-08)*Y-
     +     3.82673275931962D-07 )*Y+2.98006900751543D-06 )*Y-
     +     2.00718990300052D-05 )*Y+1.13876001386361D-04 )*Y-
     +     5.23627942443563D-04 )*Y+1.83524565118203D-03 )*Y-
     +     4.37785737450783D-03 )*Y+5.36963805223095D-03
      RETURN
C
  140 IF (X .GT. 10.0D+00) GO TO 160
C     X=5.0 TO 10.0                              NROOTS = 5
      Y = X-7.5D+00
      RT1 = ((((((((-1.13825201010775D-14*Y+1.89737681670375D-13)*Y-
     +     4.81561201185876D-12 )*Y+1.56666512163407D-10 )*Y-
     +     3.73782213255083D-09 )*Y+9.15858355075147D-08 )*Y-
     +     2.13775073585629D-06 )*Y+4.56547356365536D-05 )*Y-
     +     8.68003909323740D-04 )*Y+1.22703754069176D-02
      RT2 = (((((((((-3.67160504428358D-15*Y+1.27876280158297D-14)*Y-
     +     1.296476623788D-12)*Y+1.477175434354D-11)*Y+
     +     5.464102147892D-10)*Y-2.42538340602723D-08 )*Y+
     +     8.20460740637617D-07 )*Y-2.20379304598661D-05 )*Y+
     +     4.90295372978785D-04 )*Y-9.14294111576119D-03 )*Y+
     +     1.22590403403690D-01
      RT3 = ((((((((( 1.39017367502123D-14*Y-6.96391385426890D-13)*Y+
     +     1.176946020731D-12)*Y+1.725627235645D-10)*Y-
     +     3.686383856300D-09)*Y+2.87495324207095D-08 )*Y+
     +     1.71307311000282D-06 )*Y-7.94273603184629D-05 )*Y+
     +     2.00938064965897D-03 )*Y-3.63329491677178D-02 )*Y+
     +     4.34393683888443D-01
      RT4 = ((((((((((-1.27815158195209D-14*Y+1.99910415869821D-14)*Y+
     +     3.753542914426D-12)*Y-2.708018219579D-11)*Y-
     +     1.190574776587D-09)*Y+1.106696436509D-08)*Y+
     +     3.954955671326D-07)*Y-4.398596059588D-06)*Y-
     +     2.01087998907735D-04 )*Y+7.89092425542937D-03 )*Y-
     +     1.42056749162695D-01 )*Y+1.39964149420683D+00
      RT5 = ((((((((((-1.19442341030461D-13*Y-2.34074833275956D-12)*Y+
     +     6.861649627426D-12)*Y+6.082671496226D-10)*Y+
     +     5.381160105420D-09)*Y-6.253297138700D-08)*Y-
     +     2.135966835050D-06)*Y-2.373394341886D-05)*Y+
     +     2.88711171412814D-06 )*Y+4.85221195290753D-02 )*Y-
     +     1.04346091985269D+00 )*Y+7.89901551676692D+00
      WW1 = ((((((((( 7.95526040108997D-15*Y-2.48593096128045D-13)*Y+
     +     4.761246208720D-12)*Y-9.535763686605D-11)*Y+
     +     2.225273630974D-09)*Y-4.49796778054865D-08 )*Y+
     +     9.17812870287386D-07 )*Y-1.86764236490502D-05 )*Y+
     +     3.76807779068053D-04 )*Y-8.10456360143408D-03 )*Y+
     +     2.01097936411496D-01
      WW2 = ((((((((((( 1.25678686624734D-15*Y-2.34266248891173D-14)*Y+
     +     3.973252415832D-13)*Y-6.830539401049D-12)*Y+
     +     1.140771033372D-10)*Y-1.82546185762009D-09 )*Y+
     +     2.77209637550134D-08 )*Y-4.01726946190383D-07 )*Y+
     +     5.48227244014763D-06 )*Y-6.95676245982121D-05 )*Y+
     +     8.05193921815776D-04 )*Y-8.15528438784469D-03 )*Y+
     +     9.71769901268114D-02
      WW3 = ((((((((((((-8.20929494859896D-16*Y+1.37356038393016D-14)*Y-
     +     2.022863065220D-13)*Y+3.058055403795D-12)*Y-
     +     4.387890955243D-11)*Y+5.923946274445D-10)*Y-
     +     7.503659964159D-09)*Y+8.851599803902D-08)*Y-
     +     9.65561998415038D-07 )*Y+9.60884622778092D-06 )*Y-
     +     8.56551787594404D-05 )*Y+6.66057194311179D-04 )*Y-
     +     4.17753183902198D-03 )*Y+2.25443826852447D-02
      WW4 = ((((((((((((((-1.08764612488790D-17*Y+1.85299909689937D-16)*
     +     Y-2.730195628655D-15)*Y+4.127368817265D-14)*Y-
     +     5.881379088074D-13)*Y+7.805245193391D-12)*Y-
     +     9.632707991704D-11)*Y+1.099047050624D-09)*Y-
     +     1.15042731790748D-08 )*Y+1.09415155268932D-07 )*Y-
     +     9.33687124875935D-07 )*Y+7.02338477986218D-06 )*Y-
     +     4.53759748787756D-05 )*Y+2.41722511389146D-04 )*Y-
     +     9.75935943447037D-04 )*Y+2.57520532789644D-03
      WW5 = ((((((((((((((( 7.28996979748849D-19*Y-1.26518146195173D-17)
     +     *Y+1.886145834486D-16)*Y-2.876728287383D-15)*Y+
     +     4.114588668138D-14)*Y-5.44436631413933D-13 )*Y+
     +     6.64976446790959D-12 )*Y-7.44560069974940D-11 )*Y+
     +     7.57553198166848D-10 )*Y-6.92956101109829D-09 )*Y+
     +     5.62222859033624D-08 )*Y-3.97500114084351D-07 )*Y+
     +     2.39039126138140D-06 )*Y-1.18023950002105D-05 )*Y+
     +     4.52254031046244D-05 )*Y-1.21113782150370D-04 )*Y+
     +     1.75013126731224D-04
      RETURN
C
C     X=10.0 TO 15.0                             NROOTS = 5
  160 Y = X-12.5D+00
      RT1 = ((((((((((-4.16387977337393D-17*Y+7.20872997373860D-16)*Y+
     +     1.395993802064D-14)*Y+3.660484641252D-14)*Y-
     +     4.154857548139D-12)*Y+2.301379846544D-11)*Y-
     +     1.033307012866D-09)*Y+3.997777641049D-08)*Y-
     +     9.35118186333939D-07 )*Y+2.38589932752937D-05 )*Y-
     +     5.35185183652937D-04 )*Y+8.85218988709735D-03
      RT2 = ((((((((((-4.56279214732217D-16*Y+6.24941647247927D-15)*Y+
     +     1.737896339191D-13)*Y+8.964205979517D-14)*Y-
     +     3.538906780633D-11)*Y+9.561341254948D-11)*Y-
     +     9.772831891310D-09)*Y+4.240340194620D-07)*Y-
     +     1.02384302866534D-05 )*Y+2.57987709704822D-04 )*Y-
     +     5.54735977651677D-03 )*Y+8.68245143991948D-02
      RT3 = ((((((((((-2.52879337929239D-15*Y+2.13925810087833D-14)*Y+
     +     7.884307667104D-13)*Y-9.023398159510D-13)*Y-
     +     5.814101544957D-11)*Y-1.333480437968D-09)*Y-
     +     2.217064940373D-08)*Y+1.643290788086D-06)*Y-
     +     4.39602147345028D-05 )*Y+1.08648982748911D-03 )*Y-
     +     2.13014521653498D-02 )*Y+2.94150684465425D-01
      RT4 = ((((((((((-6.42391438038888D-15*Y+5.37848223438815D-15)*Y+
     +     8.960828117859D-13)*Y+5.214153461337D-11)*Y-
     +     1.106601744067D-10)*Y-2.007890743962D-08)*Y+
     +     1.543764346501D-07)*Y+4.520749076914D-06)*Y-
     +     1.88893338587047D-04 )*Y+4.73264487389288D-03 )*Y-
     +     7.91197893350253D-02 )*Y+8.60057928514554D-01
      RT5 = (((((((((((-2.24366166957225D-14*Y+4.87224967526081D-14)*Y+
     +     5.587369053655D-12)*Y-3.045253104617D-12)*Y-
     +     1.223983883080D-09)*Y-2.05603889396319D-09 )*Y+
     +     2.58604071603561D-07 )*Y+1.34240904266268D-06 )*Y-
     +     5.72877569731162D-05 )*Y-9.56275105032191D-04 )*Y+
     +     4.23367010370921D-02 )*Y-5.76800927133412D-01 )*Y+
     +     3.87328263873381D+00
      WW1 = ((((((((( 8.98007931950169D-15*Y+7.25673623859497D-14)*Y+
     +     5.851494250405D-14)*Y-4.234204823846D-11)*Y+
     +     3.911507312679D-10)*Y-9.65094802088511D-09 )*Y+
     +     3.42197444235714D-07 )*Y-7.51821178144509D-06 )*Y+
     +     1.94218051498662D-04 )*Y-5.38533819142287D-03 )*Y+
     +     1.68122596736809D-01
      WW2 = ((((((((((-1.05490525395105D-15*Y+1.96855386549388D-14)*Y-
     +     5.500330153548D-13)*Y+1.003849567976D-11)*Y-
     +     1.720997242621D-10)*Y+3.533277061402D-09)*Y-
     +     6.389171736029D-08)*Y+1.046236652393D-06)*Y-
     +     1.73148206795827D-05 )*Y+2.57820531617185D-04 )*Y-
     +     3.46188265338350D-03 )*Y+7.03302497508176D-02
      WW3 = ((((((((((( 3.60020423754545D-16*Y-6.24245825017148D-15)*Y+
     +     9.945311467434D-14)*Y-1.749051512721D-12)*Y+
     +     2.768503957853D-11)*Y-4.08688551136506D-10 )*Y+
     +     6.04189063303610D-09 )*Y-8.23540111024147D-08 )*Y+
     +     1.01503783870262D-06 )*Y-1.20490761741576D-05 )*Y+
     +     1.26928442448148D-04 )*Y-1.05539461930597D-03 )*Y+
     +     1.15543698537013D-02
      WW4 = ((((((((((((( 2.51163533058925D-18*Y-4.31723745510697D-17)*
     +     Y+6.557620865832D-16)*Y-1.016528519495D-14)*Y+
     +     1.491302084832D-13)*Y-2.06638666222265D-12 )*Y+
     +     2.67958697789258D-11 )*Y-3.23322654638336D-10 )*Y+
     +     3.63722952167779D-09 )*Y-3.75484943783021D-08 )*Y+
     +     3.49164261987184D-07 )*Y-2.92658670674908D-06 )*Y+
     +     2.12937256719543D-05 )*Y-1.19434130620929D-04 )*Y+
     +     6.45524336158384D-04
      WW5 = ((((((((((((((-1.29043630202811D-19*Y+2.16234952241296D-18)*
     +     Y-3.107631557965D-17)*Y+4.570804313173D-16)*Y-
     +     6.301348858104D-15)*Y+8.031304476153D-14)*Y-
     +     9.446196472547D-13)*Y+1.018245804339D-11)*Y-
     +     9.96995451348129D-11 )*Y+8.77489010276305D-10 )*Y-
     +     6.84655877575364D-09 )*Y+4.64460857084983D-08 )*Y-
     +     2.66924538268397D-07 )*Y+1.24621276265907D-06 )*Y-
     +     4.30868944351523D-06 )*Y+9.94307982432868D-06
      RETURN
C
  180 IF (X .GT. 25.0D+00) GO TO 220
      IF (X .GT. 20.0D+00) GO TO 200
C     X=15.0 TO 20.0                             NROOTS = 5
      Y = X-17.5D+00
      RT1 = (((((((((( 1.91875764545740D-16*Y+7.8357401095707D-16)*Y-
     +     3.260875931644D-14)*Y-1.186752035569D-13)*Y+
     +     4.275180095653D-12)*Y+3.357056136731D-11)*Y-
     +     1.123776903884D-09)*Y+1.231203269887D-08)*Y-
     +     3.99851421361031D-07 )*Y+1.45418822817771D-05 )*Y-
     +     3.49912254976317D-04 )*Y+6.67768703938812D-03
      RT2 = (((((((((( 2.02778478673555D-15*Y+1.01640716785099D-14)*Y-
     +     3.385363492036D-13)*Y-1.615655871159D-12)*Y+
     +     4.527419140333D-11)*Y+3.853670706486D-10)*Y-
     +     1.184607130107D-08)*Y+1.347873288827D-07)*Y-
     +     4.47788241748377D-06 )*Y+1.54942754358273D-04 )*Y-
     +     3.55524254280266D-03 )*Y+6.44912219301603D-02
      RT3 = (((((((((( 7.79850771456444D-15*Y+6.00464406395001D-14)*Y-
     +     1.249779730869D-12)*Y-1.020720636353D-11)*Y+
     +     1.814709816693D-10)*Y+1.766397336977D-09)*Y-
     +     4.603559449010D-08)*Y+5.863956443581D-07)*Y-
     +     2.03797212506691D-05 )*Y+6.31405161185185D-04 )*Y-
     +     1.30102750145071D-02 )*Y+2.10244289044705D-01
      RT4 = (((((((((((-2.92397030777912D-15*Y+1.94152129078465D-14)*Y+
     +     4.859447665850D-13)*Y-3.217227223463D-12)*Y-
     +     7.484522135512D-11)*Y+7.19101516047753D-10 )*Y+
     +     6.88409355245582D-09 )*Y-1.44374545515769D-07 )*Y+
     +     2.74941013315834D-06 )*Y-1.02790452049013D-04 )*Y+
     +     2.59924221372643D-03 )*Y-4.35712368303551D-02 )*Y+
     +     5.62170709585029D-01
      RT5 = ((((((((((( 1.17976126840060D-14*Y+1.24156229350669D-13)*Y-
     +     3.892741622280D-12)*Y-7.755793199043D-12)*Y+
     +     9.492190032313D-10)*Y-4.98680128123353D-09 )*Y-
     +     1.81502268782664D-07 )*Y+2.69463269394888D-06 )*Y+
     +     2.50032154421640D-05 )*Y-1.33684303917681D-03 )*Y+
     +     2.29121951862538D-02 )*Y-2.45653725061323D-01 )*Y+
     +     1.89999883453047D+00
      WW1 = (((((((((( 1.74841995087592D-15*Y-6.95671892641256D-16)*Y-
     +     3.000659497257D-13)*Y+2.021279817961D-13)*Y+
     +     3.853596935400D-11)*Y+1.461418533652D-10)*Y-
     +     1.014517563435D-08)*Y+1.132736008979D-07)*Y-
     +     2.86605475073259D-06 )*Y+1.21958354908768D-04 )*Y-
     +     3.86293751153466D-03 )*Y+1.45298342081522D-01
      WW2 = ((((((((((-1.11199320525573D-15*Y+1.85007587796671D-15)*Y+
     +     1.220613939709D-13)*Y+1.275068098526D-12)*Y-
     +     5.341838883262D-11)*Y+6.161037256669D-10)*Y-
     +     1.009147879750D-08)*Y+2.907862965346D-07)*Y-
     +     6.12300038720919D-06 )*Y+1.00104454489518D-04 )*Y-
     +     1.80677298502757D-03 )*Y+5.78009914536630D-02
      WW3 = ((((((((((-9.49816486853687D-16*Y+6.67922080354234D-15)*Y+
     +     2.606163540537D-15)*Y+1.983799950150D-12)*Y-
     +     5.400548574357D-11)*Y+6.638043374114D-10)*Y-
     +     8.799518866802D-09)*Y+1.791418482685D-07)*Y-
     +     2.96075397351101D-06 )*Y+3.38028206156144D-05 )*Y-
     +     3.58426847857878D-04 )*Y+8.39213709428516D-03
      WW4 = ((((((((((( 1.33829971060180D-17*Y-3.44841877844140D-16)*Y+
     +     4.745009557656D-15)*Y-6.033814209875D-14)*Y+
     +     1.049256040808D-12)*Y-1.70859789556117D-11 )*Y+
     +     2.15219425727959D-10 )*Y-2.52746574206884D-09 )*Y+
     +     3.27761714422960D-08 )*Y-3.90387662925193D-07 )*Y+
     +     3.46340204593870D-06 )*Y-2.43236345136782D-05 )*Y+
     +     3.54846978585226D-04
      WW5 = ((((((((((((( 2.69412277020887D-20*Y-4.24837886165685D-19)*
     +     Y+6.030500065438D-18)*Y-9.069722758289D-17)*Y+
     +     1.246599177672D-15)*Y-1.56872999797549D-14 )*Y+
     +     1.87305099552692D-13 )*Y-2.09498886675861D-12 )*Y+
     +     2.11630022068394D-11 )*Y-1.92566242323525D-10 )*Y+
     +     1.62012436344069D-09 )*Y-1.23621614171556D-08 )*Y+
     +     7.72165684563049D-08 )*Y-3.59858901591047D-07 )*Y+
     +     2.43682618601000D-06
      RETURN
C
C     X=20.0 TO 25.0                             NROOTS = 5
  200 Y = X-22.5D+00
      RT1 = (((((((((-1.13927848238726D-15*Y+7.39404133595713D-15)*Y+
     +     1.445982921243D-13)*Y-2.676703245252D-12)*Y+
     +     5.823521627177D-12)*Y+2.17264723874381D-10 )*Y+
     +     3.56242145897468D-09 )*Y-3.03763737404491D-07 )*Y+
     +     9.46859114120901D-06 )*Y-2.30896753853196D-04 )*Y+
     +     5.24663913001114D-03
      RT2 = (((((((((( 2.89872355524581D-16*Y-1.22296292045864D-14)*Y+
     +     6.184065097200D-14)*Y+1.649846591230D-12)*Y-
     +     2.729713905266D-11)*Y+3.709913790650D-11)*Y+
     +     2.216486288382D-09)*Y+4.616160236414D-08)*Y-
     +     3.32380270861364D-06 )*Y+9.84635072633776D-05 )*Y-
     +     2.30092118015697D-03 )*Y+5.00845183695073D-02
      RT3 = (((((((((( 1.97068646590923D-15*Y-4.89419270626800D-14)*Y+
     +     1.136466605916D-13)*Y+7.546203883874D-12)*Y-
     +     9.635646767455D-11)*Y-8.295965491209D-11)*Y+
     +     7.534109114453D-09)*Y+2.699970652707D-07)*Y-
     +     1.42982334217081D-05 )*Y+3.78290946669264D-04 )*Y-
     +     8.03133015084373D-03 )*Y+1.58689469640791D-01
      RT4 = (((((((((( 1.33642069941389D-14*Y-1.55850612605745D-13)*Y-
     +     7.522712577474D-13)*Y+3.209520801187D-11)*Y-
     +     2.075594313618D-10)*Y-2.070575894402D-09)*Y+
     +     7.323046997451D-09)*Y+1.851491550417D-06)*Y-
     +     6.37524802411383D-05 )*Y+1.36795464918785D-03 )*Y-
     +     2.42051126993146D-02 )*Y+3.97847167557815D-01
      RT5 = ((((((((((-6.07053986130526D-14*Y+1.04447493138843D-12)*Y-
     +     4.286617818951D-13)*Y-2.632066100073D-10)*Y+
     +     4.804518986559D-09)*Y-1.835675889421D-08)*Y-
     +     1.068175391334D-06)*Y+3.292234974141D-05)*Y-
     +     5.94805357558251D-04 )*Y+8.29382168612791D-03 )*Y-
     +     9.93122509049447D-02 )*Y+1.09857804755042D+00
      WW1 = (((((((((-9.10338640266542D-15*Y+1.00438927627833D-13)*Y+
     +     7.817349237071D-13)*Y-2.547619474232D-11)*Y+
     +     1.479321506529D-10)*Y+1.52314028857627D-09 )*Y+
     +     9.20072040917242D-09 )*Y-2.19427111221848D-06 )*Y+
     +     8.65797782880311D-05 )*Y-2.82718629312875D-03 )*Y+
     +     1.28718310443295D-01
      WW2 = ((((((((( 5.52380927618760D-15*Y-6.43424400204124D-14)*Y-
     +     2.358734508092D-13)*Y+8.261326648131D-12)*Y+
     +     9.229645304956D-11)*Y-5.68108973828949D-09 )*Y+
     +     1.22477891136278D-07 )*Y-2.11919643127927D-06 )*Y+
     +     4.23605032368922D-05 )*Y-1.14423444576221D-03 )*Y+
     +     5.06607252890186D-02
      WW3 = ((((((((( 3.99457454087556D-15*Y-5.11826702824182D-14)*Y-
     +     4.157593182747D-14)*Y+4.214670817758D-12)*Y+
     +     6.705582751532D-11)*Y-3.36086411698418D-09 )*Y+
     +     6.07453633298986D-08 )*Y-7.40736211041247D-07 )*Y+
     +     8.84176371665149D-06 )*Y-1.72559275066834D-04 )*Y+
     +     7.16639814253567D-03
      WW4 = (((((((((((-2.14649508112234D-18*Y-2.45525846412281D-18)*Y+
     +     6.126212599772D-16)*Y-8.526651626939D-15)*Y+
     +     4.826636065733D-14)*Y-3.39554163649740D-13 )*Y+
     +     1.67070784862985D-11 )*Y-4.42671979311163D-10 )*Y+
     +     6.77368055908400D-09 )*Y-7.03520999708859D-08 )*Y+
     +     6.04993294708874D-07 )*Y-7.80555094280483D-06 )*Y+
     +     2.85954806605017D-04
      WW5 = ((((((((((((-5.63938733073804D-21*Y+6.92182516324628D-20)*Y-
     +     1.586937691507D-18)*Y+3.357639744582D-17)*Y-
     +     4.810285046442D-16)*Y+5.386312669975D-15)*Y-
     +     6.117895297439D-14)*Y+8.441808227634D-13)*Y-
     +     1.18527596836592D-11 )*Y+1.36296870441445D-10 )*Y-
     +     1.17842611094141D-09 )*Y+7.80430641995926D-09 )*Y-
     +     5.97767417400540D-08 )*Y+1.65186146094969D-06
      RETURN
C
  220 WW1 = SQRT(PIE4/X)
      IF (X .GT. 40.0D+00) GO TO 240
C     X=25.0 TO 40.0                             NROOTS = 5
      E = EXP(-X)
      RT1 = ((((((((-1.73363958895356D-06*X+1.19921331441483D-04)*X -
     +     1.59437614121125D-02)*X+1.13467897349442D+00)*X -
     +     4.47216460864586D+01)*X+1.06251216612604D+03)*X -
     +     1.52073917378512D+04)*X+1.20662887111273D+05)*X -
     +     4.07186366852475D+05)*E + R15/(X-R15)
      RT2 = ((((((((-1.60102542621710D-05*X+1.10331262112395D-03)*X -
     +     1.50043662589017D-01)*X+1.05563640866077D+01)*X -
     +     4.10468817024806D+02)*X+9.62604416506819D+03)*X -
     +     1.35888069838270D+05)*X+1.06107577038340D+06)*X -
     +     3.51190792816119D+06)*E + R25/(X-R25)
      RT3 = ((((((((-4.48880032128422D-05*X+2.69025112122177D-03)*X -
     +     4.01048115525954D-01)*X+2.78360021977405D+01)*X -
     +     1.04891729356965D+03)*X+2.36985942687423D+04)*X -
     +     3.19504627257548D+05)*X+2.34879693563358D+06)*X -
     +     7.16341568174085D+06)*E + R35/(X-R35)
      RT4 = ((((((((-6.38526371092582D-05*X-2.29263585792626D-03)*X -
     +     7.65735935499627D-02)*X+9.12692349152792D+00)*X -
     +     2.32077034386717D+02)*X+2.81839578728845D+02)*X +
     +     9.59529683876419D+04)*X-1.77638956809518D+06)*X +
     +     1.02489759645410D+07)*E + R45/(X-R45)
      RT5 = ((((((((-3.59049364231569D-05*X-2.25963977930044D-02)*X +
     +     1.12594870794668D+00)*X-4.56752462103909D+01)*X +
     +     1.05804526830637D+03)*X-1.16003199605875D+04)*X -
     +     4.07297627297272D+04)*X+2.22215528319857D+06)*X -
     +     1.61196455032613D+07)*E + R55/(X-R55)
      WW5 = (((((((((-4.61100906133970D-10*X+1.43069932644286D-07)*X -
     +     1.63960915431080D-05)*X+1.15791154612838D-03)*X -
     +     5.30573476742071D-02)*X+1.61156533367153D+00)*X -
     +     3.23248143316007D+01)*X+4.12007318109157D+02)*X -
     +     3.02260070158372D+03)*X+9.71575094154768D+03)*E + W55*WW1
      WW4 = (((((((((-2.40799435809950D-08*X+8.12621667601546D-06)*X -
     +     9.04491430884113D-04)*X+6.37686375770059D-02)*X -
     +     2.96135703135647D+00)*X+9.15142356996330D+01)*X -
     +     1.86971865249111D+03)*X+2.42945528916947D+04)*X -
     +     1.81852473229081D+05)*X+5.96854758661427D+05)*E + W45*WW1
      WW3 = (((((((( 1.83574464457207D-05*X-1.54837969489927D-03)*X +
     +     1.18520453711586D-01)*X-6.69649981309161D+00)*X +
     +     2.44789386487321D+02)*X-5.68832664556359D+03)*X +
     +     8.14507604229357D+04)*X-6.55181056671474D+05)*X +
     +     2.26410896607237D+06)*E + W35*WW1
      WW2 = (((((((( 2.77778345870650D-05*X-2.22835017655890D-03)*X +
     +     1.61077633475573D-01)*X-8.96743743396132D+00)*X +
     +     3.28062687293374D+02)*X-7.65722701219557D+03)*X +
     +     1.10255055017664D+05)*X-8.92528122219324D+05)*X +
     +     3.10638627744347D+06)*E + W25*WW1
      WW1 = WW1-0.01962D+00*E-WW2-WW3-WW4-WW5
      RETURN
C
  240 IF (X .GT. 59.0D+00) GO TO 260
C     X=40.0 TO 59.0                             NROOTS = 5
      XXX = X**3
      E = XXX*EXP(-X)
      RT1 = (((-2.43758528330205D-02*X+2.07301567989771D+00)*X -
     +     6.45964225381113D+01)*X+7.14160088655470D+02)*E + R15/(X-R15)
      RT2 = (((-2.28861955413636D-01*X+1.93190784733691D+01)*X -
     +     5.99774730340912D+02)*X+6.61844165304871D+03)*E + R25/(X-R25)
      RT3 = (((-6.95053039285586D-01*X+5.76874090316016D+01)*X -
     +     1.77704143225520D+03)*X+1.95366082947811D+04)*E + R35/(X-R35)
      RT4 = (((-1.58072809087018D+00*X+1.27050801091948D+02)*X -
     +     3.86687350914280D+03)*X+4.23024828121420D+04)*E + R45/(X-R45)
      RT5 = (((-3.33963830405396D+00*X+2.51830424600204D+02)*X -
     +     7.57728527654961D+03)*X+8.21966816595690D+04)*E + R55/(X-R55)
      E = XXX*E
      WW5 = (( 1.35482430510942D-08*X-3.27722199212781D-07)*X +
     +     2.41522703684296D-06)*E + W55*WW1
      WW4 = (( 1.23464092261605D-06*X-3.55224564275590D-05)*X +
     +     3.03274662192286D-04)*E + W45*WW1
      WW3 = (( 1.34547929260279D-05*X-4.19389884772726D-04)*X +
     +     3.87706687610809D-03)*E + W35*WW1
      WW2 = (( 2.09539509123135D-05*X-6.87646614786982D-04)*X +
     +     6.68743788585688D-03)*E + W25*WW1
      WW1 = WW1-WW2-WW3-WW4-WW5
      RETURN
C
C     X=59.0 TO INFINITY                         NROOTS = 5
  260 RT1 = R15/(X-R15)
      RT2 = R25/(X-R25)
      RT3 = R35/(X-R35)
      RT4 = R45/(X-R45)
      RT5 = R55/(X-R55)
      WW2 = W25*WW1
      WW3 = W35*WW1
      WW4 = W45*WW1
      WW5 = W55*WW1
      WW1 = WW1-WW2-WW3-WW4-WW5
      RETURN
      END
c-------------------------------------------------------------------
C*MODULE RYSPOL  *DECK ROOT6
      SUBROUTINE xROOT6
CDIR$ INLINEALWAYS XROOT6
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION C(10,10),S(10,10),A(10),RT(10),R(9,9),W(9,9),FF(19)
C
      COMMON /xROOT  / XX,UF(13),WF(13),NROOTS
C
      PARAMETER (PT5=0.5D+00, ZERO=0.0D+00, ONE=1.0D+00, FOUR=4.0D+00)
C
C     ITH ROOT OF THE JTH RYS POLYNOMIAL IS RETURNED IN R(I,J) WITH
C     THE CORRESPONDING WEIGHT FACTOR IN W(I,J).   J=1,2,...,N
C     THIS VERSION USES CHRISTOFFEL FORMULA FOR WEIGHTS, AND IT IS
C     COMPLETELY GENERAL UP TO NROOTS=9.  SEE THE REFERENCE
C        "NUMERICAL INTEGRATION USING RYS POLYNOMIALS"
C                   H.F.KING AND M.DUPUIS
C              J.COMPUT.PHYS. 21, 144-165(1976)
C
      N=NROOTS
      X=XX
      IF(N.LT.2) N=2
      N1=N+1
      NN=N+N
      CALL xRYSFUN(X,NN,FF)
C
      DO 11 I=1,N1
         DO 10 J=1,N1
            S(I,J)=FF(I+J-1)
   10    CONTINUE
   11 CONTINUE
      CALL xRYSSMT(C,S,N1)
C
      DO 20 I=1,N
         DO 20 J=1,I
            W(I,J)= ZERO
            R(I,J)= ZERO
   20    CONTINUE
   21 CONTINUE
      WSUM=FF(1)
      W(1,1)=WSUM
      R(1,1)=FF(2)/WSUM
      DUM= SQRT(C(2,3)**2-FOUR *C(1,3)*C(3,3))
      R(1,2)= PT5*(-C(2,3)-DUM)/C(3,3)
      R(2,2)= PT5*(-C(2,3)+DUM)/C(3,3)
      IF(N.EQ.2) GO TO 70
      DO 25 I=3,N1
   25 RT(I)= ONE
      RT(1)=R(1,2)
      RT(2)=R(2,2)
C
      DO 60 K=3,N
         K1=K+1
         DO 30 I=1,K1
            A(I)=C(I,K1)
   30    CONTINUE
         CALL xRYSNOD(A,RT,K)
         DO 50 I=1,K
            R(I,K)=RT(I)
   50    CONTINUE
   60 CONTINUE
C
   70 CONTINUE
      DO 150 K=2,N
         JMAX=K-1
         DO 140 I=1,K
            ROOT=R(I,K)
            DUM = ONE/FF(1)
            DO 110 J=1,JMAX
               J1=J+1
               POLY=C(J1,J1)
               DO 100 M=1,J
                  POLY=POLY*ROOT+C(J1-M,J1)
  100          CONTINUE
               DUM=DUM+POLY*POLY
  110       CONTINUE
            W(I,K) = ONE/DUM
  140    CONTINUE
  150 CONTINUE
C
      DO 160 K=1,NROOTS
         DUM=R(K,NROOTS)
         UF(K)=DUM/(ONE-DUM)
         WF(K)=W(K,NROOTS)
  160 CONTINUE
      RETURN
      END
c------------------------------------------------------------------
C*MODULE RYSPOL  *DECK RT123
      SUBROUTINE xRT123
CDIR$ INLINEALWAYS XRT123
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C             *****   VERSION FEBRUARY 13,1975   *****
      COMMON /xROOT  / X,U(13),W(13),NROOTS
      EQUIVALENCE (U(1),RT1),(U(2),RT2),(U(3),RT3),(U(4),RT4),(U(5),RT5)
      EQUIVALENCE (W(1),WW1),(W(2),WW2),(W(3),WW3),(W(4),WW4),(W(5),WW5)
      DATA R12,PIE4/2.75255128608411D-01, 7.85398163397448D-01/
      DATA R22,W22/ 2.72474487139158D+00, 9.17517095361369D-02/
      DATA R13/     1.90163509193487D-01/
      DATA R23,W23/ 1.78449274854325D+00, 1.77231492083829D-01/
      DATA R33,W33/ 5.52534374226326D+00, 5.11156880411248D-03/
      IF (X .GT. 5.0D+00) GO TO 400
      IF (X .GT. 1.0D+00) GO TO 280
      IF (X .GT. 3.0D-07) GO TO 180
C     X IS APPROXIMATELY ZERO.         NROOTS=1,2, OR 3
      IF (NROOTS-2) 120,140,160
  120 RT1 = 0.5D+00 -X/5.0D+00
      WW1 = 1.0D+00 -X/3.0D+00
      RETURN
  140 RT1 = 1.30693606237085D-01 -2.90430236082028D-02 *X
      RT2 = 2.86930639376291D+00 -6.37623643058102D-01 *X
      WW1 = 6.52145154862545D-01 -1.22713621927067D-01 *X
      WW2 = 3.47854845137453D-01 -2.10619711404725D-01 *X
      RETURN
  160 RT1 = 6.03769246832797D-02 -9.28875764357368D-03 *X
      RT2 = 7.76823355931043D-01 -1.19511285527878D-01 *X
      RT3 = 6.66279971938567D+00 -1.02504611068957D+00 *X
      WW1 = 4.67913934572691D-01 -5.64876917232519D-02 *X
      WW2 = 3.60761573048137D-01 -1.49077186455208D-01 *X
      WW3 = 1.71324492379169D-01 -1.27768455150979D-01 *X
      RETURN
C     X = 0.0 TO 1.0                   NROOTS=1,2, OR 3
  180 IF (NROOTS .EQ. 3) GO TO 220
      F1 = ((((((((-8.36313918003957D-08*X+1.21222603512827D-06 )*X-
     +     1.15662609053481D-05 )*X+9.25197374512647D-05 )*X-
     +     6.40994113129432D-04 )*X+3.78787044215009D-03 )*X-
     +     1.85185172458485D-02 )*X+7.14285713298222D-02 )*X-
     +     1.99999999997023D-01 )*X+3.33333333333318D-01
      WW1 = (X+X)*F1+EXP(-X)
      IF (NROOTS .EQ. 2) GO TO 200
      RT1 = F1/(WW1-F1)
      RETURN
  200 RT1 = (((((((-2.35234358048491D-09*X+2.49173650389842D-08)*X-
     +     4.558315364581D-08)*X-2.447252174587D-06)*X+
     +     4.743292959463D-05)*X-5.33184749432408D-04 )*X+
     +     4.44654947116579D-03 )*X-2.90430236084697D-02 )*X+
     +     1.30693606237085D-01
      RT2 = (((((((-2.47404902329170D-08*X+2.36809910635906D-07)*X+
     +     1.835367736310D-06)*X-2.066168802076D-05)*X-
     +     1.345693393936D-04)*X-5.88154362858038D-05 )*X+
     +     5.32735082098139D-02 )*X-6.37623643056745D-01 )*X+
     +     2.86930639376289D+00
      WW2 = ((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1 = WW1-WW2
      RETURN
  220 RT1 = ((((((-5.10186691538870D-10*X+2.40134415703450D-08)*X-
     +     5.01081057744427D-07 )*X+7.58291285499256D-06 )*X-
     +     9.55085533670919D-05 )*X+1.02893039315878D-03 )*X-
     +     9.28875764374337D-03 )*X+6.03769246832810D-02
      RT2 = ((((((-1.29646524960555D-08*X+7.74602292865683D-08)*X+
     +     1.56022811158727D-06 )*X-1.58051990661661D-05 )*X-
     +     3.30447806384059D-04 )*X+9.74266885190267D-03 )*X-
     +     1.19511285526388D-01 )*X+7.76823355931033D-01
      RT3 = ((((((-9.28536484109606D-09*X-3.02786290067014D-07)*X-
     +     2.50734477064200D-06 )*X-7.32728109752881D-06 )*X+
     +     2.44217481700129D-04 )*X+4.94758452357327D-02 )*X-
     +     1.02504611065774D+00 )*X+6.66279971938553D+00
      F2 = ((((((((-7.60911486098850D-08*X+1.09552870123182D-06 )*X-
     +     1.03463270693454D-05 )*X+8.16324851790106D-05 )*X-
     +     5.55526624875562D-04 )*X+3.20512054753924D-03 )*X-
     +     1.51515139838540D-02 )*X+5.55555554649585D-02 )*X-
     +     1.42857142854412D-01 )*X+1.99999999999986D-01
  240 E = EXP(-X)
      F1 = ((X+X)*F2+E)/3.0D+00
      WW1 = (X+X)*F1+E
  260 T1 = RT1/(RT1+1.0D+00)
      T2 = RT2/(RT2+1.0D+00)
      T3 = RT3/(RT3+1.0D+00)
      A2 = F2-T1*F1
      A1 = F1-T1*WW1
      WW3 = (A2-T2*A1)/((T3-T2)*(T3-T1))
      WW2 = (T3*A1-A2)/((T3-T2)*(T2-T1))
      WW1 = WW1-WW2-WW3
      RETURN
  280 IF (X .GT. 3.0D+00) GO TO 340
C     X = 1.0 TO 3.0                   NROOTS=1,2, OR 3
      Y = X-2.0D+00
      IF (NROOTS .EQ. 3) GO TO 320
      F1 = ((((((((((-1.61702782425558D-10*Y+1.96215250865776D-09 )*Y-
     +     2.14234468198419D-08 )*Y+2.17216556336318D-07 )*Y-
     +     1.98850171329371D-06 )*Y+1.62429321438911D-05 )*Y-
     +     1.16740298039895D-04 )*Y+7.24888732052332D-04 )*Y-
     +     3.79490003707156D-03 )*Y+1.61723488664661D-02 )*Y-
     +     5.29428148329736D-02 )*Y+1.15702180856167D-01
      WW1 = (X+X)*F1+EXP(-X)
      IF (NROOTS .EQ. 2) GO TO 300
      RT1 = F1/(WW1-F1)
      RETURN
  300 RT1 = (((((((((-6.36859636616415D-12*Y+8.47417064776270D-11)*Y-
     +     5.152207846962D-10)*Y-3.846389873308D-10)*Y+
     +     8.472253388380D-08)*Y-1.85306035634293D-06 )*Y+
     +     2.47191693238413D-05 )*Y-2.49018321709815D-04 )*Y+
     +     2.19173220020161D-03 )*Y-1.63329339286794D-02 )*Y+
     +     8.68085688285261D-02
      RT2 = ((((((((( 1.45331350488343D-10*Y+2.07111465297976D-09)*Y-
     +     1.878920917404D-08)*Y-1.725838516261D-07)*Y+
     +     2.247389642339D-06)*Y+9.76783813082564D-06 )*Y-
     +     1.93160765581969D-04 )*Y-1.58064140671893D-03 )*Y+
     +     4.85928174507904D-02 )*Y-4.30761584997596D-01 )*Y+
     +     1.80400974537950D+00
      WW2 = ((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1 = WW1-WW2
      RETURN
  320 RT1 = (((((((( 1.44687969563318D-12*Y+4.85300143926755D-12)*Y-
     +     6.55098264095516D-10 )*Y+1.56592951656828D-08 )*Y-
     +     2.60122498274734D-07 )*Y+3.86118485517386D-06 )*Y-
     +     5.13430986707889D-05 )*Y+6.03194524398109D-04 )*Y-
     +     6.11219349825090D-03 )*Y+4.52578254679079D-02
      RT2 = ((((((( 6.95964248788138D-10*Y-5.35281831445517D-09)*Y-
     +     6.745205954533D-08)*Y+1.502366784525D-06)*Y+
     +     9.923326947376D-07)*Y-3.89147469249594D-04 )*Y+
     +     7.51549330892401D-03 )*Y-8.48778120363400D-02 )*Y+
     +     5.73928229597613D-01
      RT3 = ((((((((-2.81496588401439D-10*Y+3.61058041895031D-09)*Y+
     +     4.53631789436255D-08 )*Y-1.40971837780847D-07 )*Y-
     +     6.05865557561067D-06 )*Y-5.15964042227127D-05 )*Y+
     +     3.34761560498171D-05 )*Y+5.04871005319119D-02 )*Y-
     +     8.24708946991557D-01 )*Y+4.81234667357205D+00
      F2 = ((((((((((-1.48044231072140D-10*Y+1.78157031325097D-09 )*Y-
     +     1.92514145088973D-08 )*Y+1.92804632038796D-07 )*Y-
     +     1.73806555021045D-06 )*Y+1.39195169625425D-05 )*Y-
     +     9.74574633246452D-05 )*Y+5.83701488646511D-04 )*Y-
     +     2.89955494844975D-03 )*Y+1.13847001113810D-02 )*Y-
     +     3.23446977320647D-02 )*Y+5.29428148329709D-02
      GO TO 240
C     X = 3.0 TO 5.0                   NROOTS =1,2, OR 3
  340 Y = X-4.0D+00
      IF (NROOTS .EQ. 3) GO TO 380
      F1 = ((((((((((-2.62453564772299D-11*Y+3.24031041623823D-10 )*Y-
     +     3.614965656163D-09)*Y+3.760256799971D-08)*Y-
     +     3.553558319675D-07)*Y+3.022556449731D-06)*Y-
     +     2.290098979647D-05)*Y+1.526537461148D-04)*Y-
     +     8.81947375894379D-04 )*Y+4.33207949514611D-03 )*Y-
     +     1.75257821619926D-02 )*Y+5.28406320615584D-02
      WW1 = (X+X)*F1+EXP(-X)
      IF (NROOTS .EQ. 2) GO TO 360
      RT1 = F1/(WW1-F1)
      RETURN
  360 RT1 = ((((((((-4.11560117487296D-12*Y+7.10910223886747D-11)*Y-
     +     1.73508862390291D-09 )*Y+5.93066856324744D-08 )*Y-
     +     9.76085576741771D-07 )*Y+1.08484384385679D-05 )*Y-
     +     1.12608004981982D-04 )*Y+1.16210907653515D-03 )*Y-
     +     9.89572595720351D-03 )*Y+6.12589701086408D-02
      RT2 = (((((((((-1.80555625241001D-10*Y+5.44072475994123D-10)*Y+
     +     1.603498045240D-08)*Y-1.497986283037D-07)*Y-
     +     7.017002532106D-07)*Y+1.85882653064034D-05 )*Y-
     +     2.04685420150802D-05 )*Y-2.49327728643089D-03 )*Y+
     +     3.56550690684281D-02 )*Y-2.60417417692375D-01 )*Y+
     +     1.12155283108289D+00
      WW2 = ((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1 = WW1-WW2
      RETURN
  380 RT1 = ((((((( 1.44265709189601D-11*Y-4.66622033006074D-10)*Y+
     +     7.649155832025D-09)*Y-1.229940017368D-07)*Y+
     +     2.026002142457D-06)*Y-2.87048671521677D-05 )*Y+
     +     3.70326938096287D-04 )*Y-4.21006346373634D-03 )*Y+
     +     3.50898470729044D-02
      RT2 = ((((((((-2.65526039155651D-11*Y+1.97549041402552D-10)*Y+
     +     2.15971131403034D-09 )*Y-7.95045680685193D-08 )*Y+
     +     5.15021914287057D-07 )*Y+1.11788717230514D-05 )*Y-
     +     3.33739312603632D-04 )*Y+5.30601428208358D-03 )*Y-
     +     5.93483267268959D-02 )*Y+4.31180523260239D-01
      RT3 = ((((((((-3.92833750584041D-10*Y-4.16423229782280D-09)*Y+
     +     4.42413039572867D-08 )*Y+6.40574545989551D-07 )*Y-
     +     3.05512456576552D-06 )*Y-1.05296443527943D-04 )*Y-
     +     6.14120969315617D-04 )*Y+4.89665802767005D-02 )*Y-
     +     6.24498381002855D-01 )*Y+3.36412312243724D+00
      F2 = ((((((((((-2.36788772599074D-11*Y+2.89147476459092D-10 )*Y-
     +     3.18111322308846D-09 )*Y+3.25336816562485D-08 )*Y-
     +     3.00873821471489D-07 )*Y+2.48749160874431D-06 )*Y-
     +     1.81353179793672D-05 )*Y+1.14504948737066D-04 )*Y-
     +     6.10614987696677D-04 )*Y+2.64584212770942D-03 )*Y-
     +     8.66415899015349D-03 )*Y+1.75257821619922D-02
      GO TO 240
  400 IF (X .GT. 15.0D+00) GO TO 560
      E = EXP(-X)
      IF (X .GT. 10.0D+00) GO TO 480
C     X = 5.0 TO 10.0                  NROOTS =1,2, OR 3
      WW1 = (((((( 4.6897511375022D-01/X-6.9955602298985D-01)/X +
     +     5.3689283271887D-01)/X-3.2883030418398D-01)/X +
     +     2.4645596956002D-01)/X-4.9984072848436D-01)/X -
     +     3.1501078774085D-06)*E + SQRT(PIE4/X)
      F1 = (WW1-E)/(X+X)
      IF (NROOTS-2) 420,440,460
  420 RT1 = F1/(WW1-F1)
      RETURN
  440 Y = X-7.5D+00
      RT1 = (((((((((((((-1.43632730148572D-16*Y+2.38198922570405D-16)*
     +     Y+1.358319618800D-14)*Y-7.064522786879D-14)*Y-
     +     7.719300212748D-13)*Y+7.802544789997D-12)*Y+
     +     6.628721099436D-11)*Y-1.775564159743D-09)*Y+
     +     1.713828823990D-08)*Y-1.497500187053D-07)*Y+
     +     2.283485114279D-06)*Y-3.76953869614706D-05 )*Y+
     +     4.74791204651451D-04 )*Y-4.60448960876139D-03 )*Y+
     +     3.72458587837249D-02
      RT2 = (((((((((((( 2.48791622798900D-14*Y-1.36113510175724D-13)*Y-
     +     2.224334349799D-12)*Y+4.190559455515D-11)*Y-
     +     2.222722579924D-10)*Y-2.624183464275D-09)*Y+
     +     6.128153450169D-08)*Y-4.383376014528D-07)*Y-
     +     2.49952200232910D-06 )*Y+1.03236647888320D-04 )*Y-
     +     1.44614664924989D-03 )*Y+1.35094294917224D-02 )*Y-
     +     9.53478510453887D-02 )*Y+5.44765245686790D-01
      WW2 = ((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1 = WW1-WW2
      RETURN
  460 F2 = (F1+F1+F1-E)/(X+X)
      Y = X-7.5D+00
      RT1 = ((((((((((( 5.74429401360115D-16*Y+7.11884203790984D-16)*Y-
     +     6.736701449826D-14)*Y-6.264613873998D-13)*Y+
     +     1.315418927040D-11)*Y-4.23879635610964D-11 )*Y+
     +     1.39032379769474D-09 )*Y-4.65449552856856D-08 )*Y+
     +     7.34609900170759D-07 )*Y-1.08656008854077D-05 )*Y+
     +     1.77930381549953D-04 )*Y-2.39864911618015D-03 )*Y+
     +     2.39112249488821D-02
      RT2 = ((((((((((( 1.13464096209120D-14*Y+6.99375313934242D-15)*Y-
     +     8.595618132088D-13)*Y-5.293620408757D-12)*Y-
     +     2.492175211635D-11)*Y+2.73681574882729D-09 )*Y-
     +     1.06656985608482D-08 )*Y-4.40252529648056D-07 )*Y+
     +     9.68100917793911D-06 )*Y-1.68211091755327D-04 )*Y+
     +     2.69443611274173D-03 )*Y-3.23845035189063D-02 )*Y+
     +     2.75969447451882D-01
      RT3 = (((((((((((( 6.66339416996191D-15*Y+1.84955640200794D-13)*Y-
     +     1.985141104444D-12)*Y-2.309293727603D-11)*Y+
     +     3.917984522103D-10)*Y+1.663165279876D-09)*Y-
     +     6.205591993923D-08)*Y+8.769581622041D-09)*Y+
     +     8.97224398620038D-06 )*Y-3.14232666170796D-05 )*Y-
     +     1.83917335649633D-03 )*Y+3.51246831672571D-02 )*Y-
     +     3.22335051270860D-01 )*Y+1.73582831755430D+00
      GO TO 260
C     X = 10.0 TO 15.0                 NROOTS=1,2, OR 3
  480 WW1 = (((-1.8784686463512D-01/X+2.2991849164985D-01)/X -
     +     4.9893752514047D-01)/X-2.1916512131607D-05)*E + SQRT(PIE4/X)
      F1 = (WW1-E)/(X+X)
      IF (NROOTS-2) 500,520,540
  500 RT1 = F1/(WW1-F1)
      RETURN
  520 RT1 = ((((-1.01041157064226D-05*X+1.19483054115173D-03)*X -
     +     6.73760231824074D-02)*X+1.25705571069895D+00)*X + (((-
     +     8.57609422987199D+03/X+5.91005939591842D+03)/X -
     +     1.70807677109425D+03)/X+2.64536689959503D+02)/X -
     +     2.38570496490846D+01)*E + R12/(X-R12)
      RT2 = ((( 3.39024225137123D-04*X-9.34976436343509D-02)*X -
     +     4.22216483306320D+00)*X + (((-2.08457050986847D+03/X -
     +     1.04999071905664D+03)/X+3.39891508992661D+02)/X -
     +     1.56184800325063D+02)/X+8.00839033297501D+00)*E + R22/(X-R22)
      WW2 = ((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1 = WW1-WW2
      RETURN
  540 F2 = (F1+F1+F1-E)/(X+X)
      Y = X-12.5D+00
      RT1 = ((((((((((( 4.42133001283090D-16*Y-2.77189767070441D-15)*Y-
     +     4.084026087887D-14)*Y+5.379885121517D-13)*Y+
     +     1.882093066702D-12)*Y-8.67286219861085D-11 )*Y+
     +     7.11372337079797D-10 )*Y-3.55578027040563D-09 )*Y+
     +     1.29454702851936D-07 )*Y-4.14222202791434D-06 )*Y+
     +     8.04427643593792D-05 )*Y-1.18587782909876D-03 )*Y+
     +     1.53435577063174D-02
      RT2 = ((((((((((( 6.85146742119357D-15*Y-1.08257654410279D-14)*Y-
     +     8.579165965128D-13)*Y+6.642452485783D-12)*Y+
     +     4.798806828724D-11)*Y-1.13413908163831D-09 )*Y+
     +     7.08558457182751D-09 )*Y-5.59678576054633D-08 )*Y+
     +     2.51020389884249D-06 )*Y-6.63678914608681D-05 )*Y+
     +     1.11888323089714D-03 )*Y-1.45361636398178D-02 )*Y+
     +     1.65077877454402D-01
      RT3 = (((((((((((( 3.20622388697743D-15*Y-2.73458804864628D-14)*Y-
     +     3.157134329361D-13)*Y+8.654129268056D-12)*Y-
     +     5.625235879301D-11)*Y-7.718080513708D-10)*Y+
     +     2.064664199164D-08)*Y-1.567725007761D-07)*Y-
     +     1.57938204115055D-06 )*Y+6.27436306915967D-05 )*Y-
     +     1.01308723606946D-03 )*Y+1.13901881430697D-02 )*Y-
     +     1.01449652899450D-01 )*Y+7.77203937334739D-01
      GO TO 260
  560 IF (X .GT. 33.0D+00) GO TO 660
C     X = 15.0 TO 33.0                 NROOTS=1,2, OR 3
      E = EXP(-X)
      WW1 = (( 1.9623264149430D-01/X-4.9695241464490D-01)/X -
     +     6.0156581186481D-05)*E + SQRT(PIE4/X)
      F1 = (WW1-E)/(X+X)
      IF (NROOTS-2) 580,600,620
  580 RT1 = F1/(WW1-F1)
      RETURN
  600 RT1 = ((((-1.14906395546354D-06*X+1.76003409708332D-04)*X -
     +     1.71984023644904D-02)*X-1.37292644149838D-01)*X + (-
     +     4.75742064274859D+01/X+9.21005186542857D+00)/X -
     +     2.31080873898939D-02)*E + R12/(X-R12)
      RT2 = ((( 3.64921633404158D-04*X-9.71850973831558D-02)*X -
     +     4.02886174850252D+00)*X + (-1.35831002139173D+02/X -
     +     8.66891724287962D+01)/X+2.98011277766958D+00)*E + R22/(X-R22)
      WW2 = ((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1 = WW1-WW2
      RETURN
  620 F2 = (F1+F1+F1-E)/(X+X)
      IF (X .GT. 20.0D+00) GO TO 640
      RT1 = ((((((-2.43270989903742D-06*X+3.57901398988359D-04)*X -
     +     2.34112415981143D-02)*X+7.81425144913975D-01)*X -
     +     1.73209218219175D+01)*X+2.43517435690398D+02)*X + (-
     +     1.97611541576986D+04/X+9.82441363463929D+03)/X -
     +     2.07970687843258D+03)*E + R13/(X-R13)
      RT2 = (((((-2.62627010965435D-04*X+3.49187925428138D-02)*X -
     +     3.09337618731880D+00)*X+1.07037141010778D+02)*X -
     +     2.36659637247087D+03)*X + ((-2.91669113681020D+06/X +
     +     1.41129505262758D+06)/X-2.91532335433779D+05)/X +
     +     3.35202872835409D+04)*E + R23/(X-R23)
      RT3 = ((((( 9.31856404738601D-05*X-2.87029400759565D-02)*X -
     +     7.83503697918455D-01)*X-1.84338896480695D+01)*X +
     +     4.04996712650414D+02)*X + (-1.89829509315154D+05/X +
     +     5.11498390849158D+04)/X-6.88145821789955D+03)*E + R33/(X-R33)
      GO TO 260
  640 RT1 = ((((-4.97561537069643D-04*X-5.00929599665316D-02)*X +
     +     1.31099142238996D+00)*X-1.88336409225481D+01)*X -
     +     6.60344754467191D+02 /X+1.64931462413877D+02)*E + R13/(X-R13)
      RT2 = ((((-4.48218898474906D-03*X-5.17373211334924D-01)*X +
     +     1.13691058739678D+01)*X-1.65426392885291D+02)*X -
     +     6.30909125686731D+03 /X+1.52231757709236D+03)*E + R23/(X-R23)
      RT3 = ((((-1.38368602394293D-02*X-1.77293428863008D+00)*X +
     +     1.73639054044562D+01)*X-3.57615122086961D+02)*X -
     +     1.45734701095912D+04 /X+2.69831813951849D+03)*E + R33/(X-R33)
      GO TO 260
C     X = 33.0 TO INFINITY             NROOTS=1,2, OR 3
  660 WW1 = SQRT(PIE4/X)
      IF (NROOTS-2) 680,700,720
  680 RT1 = 0.5D+00/(X-0.5D+00)
      RETURN
  700 IF (X .GT. 40.0D+00) GO TO 740
      E = EXP(-X)
      RT1 = (-8.78947307498880D-01*X+1.09243702330261D+01)*E + R12/(X-
     +     R12)
      RT2 = (-9.28903924275977D+00*X+8.10642367843811D+01)*E + R22/(X-
     +     R22)
      WW2 = ( 4.46857389308400D+00*X-7.79250653461045D+01)*E + W22*WW1
      WW1 = WW1-WW2
      RETURN
  720 IF (X .GT. 47.0D+00) GO TO 760
      E = EXP(-X)
      RT1 = ((-7.39058467995275D+00*X+3.21318352526305D+02)*X -
     +     3.99433696473658D+03)*E + R13/(X-R13)
      RT2 = ((-7.38726243906513D+01*X+3.13569966333873D+03)*X -
     +     3.86862867311321D+04)*E + R23/(X-R23)
      RT3 = ((-2.63750565461336D+02*X+1.04412168692352D+04)*X -
     +     1.28094577915394D+05)*E + R33/(X-R33)
      WW3 = ((( 1.52258947224714D-01*X-8.30661900042651D+00)*X +
     +     1.92977367967984D+02)*X-1.67787926005344D+03)*E + W33*WW1
      WW2 = (( 6.15072615497811D+01*X-2.91980647450269D+03)*X +
     +     3.80794303087338D+04)*E + W23*WW1
      WW1 = WW1-WW2-WW3
      RETURN
  740 RT1 = R12/(X-R12)
      RT2 = R22/(X-R22)
      WW2 = W22*WW1
      WW1 = WW1-WW2
      RETURN
  760 RT1 = R13/(X-R13)
      RT2 = R23/(X-R23)
      RT3 = R33/(X-R33)
      WW2 = W23*WW1
      WW3 = W33*WW1
      WW1 = WW1-WW2-WW3
      RETURN
      END
c-----------------------------------------------------------------
C*MODULE INT2B   *DECK GENR70
      SUBROUTINE xGENR70(IEXCH,SCHSKP)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL SCHSKP
C
C     GAUSSIAN TWO ELECTRON INTEGRAL PACKAGE
C
C     SCHSKP (SCHWARZ SKIP) MEANS ONLY IB INDEXING, AND ZEROING IS DONE.
C
C     GENR70 EVALUATES THOSE INTEGRALS INVOLVING ONLY S AND P FUNCTIONS.
C     MAIN LOOP OVER SHELLS ... SEE EXPLANATION IN ROUTINE GINPUT ...
C     LINK 301
C     ACCEPTS NUMBERS OF FOUR SHELLS: ISHELL  JSHELL  KSHELL  LSHELL
C     FINDS THEIR ANGULAR QUANTUM NUMBERS AND, BASED ON THIS
C     ORDERS SHELLS IN A STANDARD MANNER: INEW  JNEW  KNEW  LNEW
C
C     ONLY POSSIBILITIES ALLOWED FOR ANGULAR QUANTUM NUMBERS ARE THEN
C     0000  0001  0011  0101  0111  1111
C     DETERMINES TYPE OF INTEGRAL SET BASED ON THE ABOVE NUMBERS
C     CALLS THE FOLLOWING ROUTINES IN THE ORDER GIVEN (THESE ROUTINES
C     WERE LIFTED FROM GAUSSIAN76):
C
C     SINFO, TO PRESET INTEGRAL ACCURACY LIMITS
C
C     SGEOM, OBTAINS GEOMETRICAL INFORMATION ABOUT THE FOUR CENTERS
C     AND FINDS TWO SETS OF LOCAL AXES FOR CENTERS:
C     (A AND B) P SET AND (C AND D) Q SET
C
C     PINF, OBTAINS INFORMATION ABOUT GAUSSIAN FUNCTIONS CONNECTED
C     WITH THE P SET OF AXES
C
C     AT THIS POINT GENR70 OBTAINS INFORMATION ABOUT THE GAUSSIAN
C     FUNCTIONS CONNECTED WITH THE Q SET OF AXES
C
C     SP0000 TO SP1111, OBTAIN UP TO 88 INTEGRALS REFERRED
C     TO AXES A B AND Q
C
C     ROT2, ROTATES THESE INTEGRALS TO UP TO 160 INTEGRALS ON A B AND Q
C
C     TQ0011 TO TQ1111, TRANSLATE THESE INTEGRALS ON A B AND Q
C     TO UP TO 256 INTEGRALS ON A B C AND D
C
C     R30001 TO R31111, ROTATE UP TO 256 INTEGRALS ON A B C AND D
C     TO THE SAME NUMBER REFERRED TO THE FIXED SPACE AXES
C
C          THE POPLE/HEHRE ROTATED AXIS CODE WAS EXTENSIVELY REWRITTEN
C          IN JULY 2001 BY JOSE M. SIERRA OF SYNSTAR, SA, MADRID SPAIN.
C     THE CHANGES INCLUDE USING DO LOOPS INSTEAD OF UNROLLED CODE, AND
C     ROUTINES SINFO, SGEOM AND PINF HAVE BEEN IN-LINED INTO GENR70
C     ROUTINE  R30S1S ENCOMPASSES SUBS FROM R30001 TO R31111
C     ROUTINE  SP0S1S ENCOMPASSES SUBS FROM SP0000 TO SP1111
C     ROUTINE  TQ0S1S ENCOMPASSES SUBS FROM TQ0011 TO TQ1111
C
#include "gamess.par"
C
      COMMON /xB     / CO(MXSH,3)
      COMMON /xFLIPS / IB(4,3)
      COMMON /xGEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),
     4                CONP(MXG2),NGANGB
      COMMON /xGOUT  / GPOPLE(256*3),NORG
      COMMON /xMAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     2                CF(MXGTOT),CG(MXGTOT),
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /xSHLG70/ ISHELL,JSHELL,KSHELL,LSHELL,INEW,JNEW,KNEW,LNEW
      COMMON /xSHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH),
     2                NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH),
     3                NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH),
     4                NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
      COMMON /xSHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     2                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     3                NIJ,IJ,KL,IJKL
C
C     ==================================================================
C     NOTE THAT GAMGEN IS THE PLACE THAT SETS COMMON BLOCK /INTACC/
C     THIS GROUP OF STATEMENTS IS INTERRELATED FOR ROUTINES
C     GAMGEN, GENR70 AND SP0S1S IN INT2B, AND JKDG80, GRDG80 AND GRDPER
C     IN GRD2B. IF YOU CHANGE THESE, CHANGE ALL SIX PLACES EQUIVALENTLY.
C     GAMGEN HAS COMMENTS EXPLAINING THIS IN THE MOST DETAIL.
C
      PARAMETER (M3R=401)
      COMMON /xINTACC/ EI1,EI2,CUF,RFI,CUG,RGI,CUX,R(3,0:5,0:M3R)
C     ==================================================================
C
CJMS  LABELLED COMMON JMSGYH DEFINED FOR COMPUTATIONAL EFFICIENCY.
CJMS  IT IS ONLY USED IN THIS MODULE INT2B.
C
      COMMON /xJMSGYH/ G(256),H(256),RQ(0:1,0:1),SQ(0:1,0:1)
C
      DIMENSION P12(3,3),P34(3,3),P(3,3),T(3)
C
      PARAMETER (ZER=0.0D+00, ONE=1.0D+00)
      PARAMETER (PT7=0.7D+00, PT9=0.9D+00)
      PARAMETER (TENM04=1.0D-04, TENM12=1.0D-12)
C
CJMS  PARAMETER (PIDIV4=0.785398163397448279D+00)
CC    PARAMETER (PIDIV4=0.785398163397448309615660845819D+00)
CC    PARAMETER (PIDIV4=0.785398163397448D+00)
      PARAMETER (PITO52=0.349868366552497D+02)
C
CC    PIDIV4= ATAN(ONE)
CC    PI4= PIDIV4
CC    PI = PI4+PI4+PI4+PI4
CC    PITO52=(PI +PI )*PI * SQRT(PI )
      PIF= PITO52
C
      LAT= KTYPE(ISHELL)-1
      LBT= KTYPE(JSHELL)-1
      LCT= KTYPE(KSHELL)-1
      LDT= KTYPE(LSHELL)-1
      ITYPE= 1+LDT+2*(LCT+2*(LBT+2*LAT))
C
      IF(ITYPE.EQ. 1 .OR. ITYPE.EQ. 2 .OR. ITYPE.EQ. 4 .OR.
     2   ITYPE.EQ. 6 .OR. ITYPE.EQ. 8 .OR. ITYPE.EQ.16) THEN
C
C TYPES 0000,0001,0011,0101,0111,1111 ARE UNALTERED
C
         INEW= ISHELL
         JNEW= JSHELL
         KNEW= KSHELL
         LNEW= LSHELL
         LA= LAT
         LB= LBT
         LC= LCT
         LD= LDT
         IB(1,IEXCH)= 1
         IB(2,IEXCH)= 2
         IB(3,IEXCH)= 3
         IB(4,IEXCH)= 4
      ELSEIF(ITYPE.EQ.10 .OR. ITYPE.EQ.12) THEN
C
C TYPES 1001,1011 HAVE IJ SWITCHED
C
         INEW= JSHELL
         JNEW= ISHELL
         KNEW= KSHELL
         LNEW= LSHELL
         LA= LBT
         LB= LAT
         LC= LCT
         LD= LDT
         IB(1,IEXCH)= 2
         IB(2,IEXCH)= 1
         IB(3,IEXCH)= 3
         IB(4,IEXCH)= 4
      ELSEIF(ITYPE.EQ. 3 .OR. ITYPE.EQ. 7) THEN
C
C TYPES 0010,0110 HAVE KL SWITCHED
C
         INEW= ISHELL
         JNEW= JSHELL
         KNEW= LSHELL
         LNEW= KSHELL
         LA= LAT
         LB= LBT
         LC= LDT
         LD= LCT
         IB(1,IEXCH)= 1
         IB(2,IEXCH)= 2
         IB(3,IEXCH)= 4
         IB(4,IEXCH)= 3
      ELSEIF(ITYPE.EQ. 5 .OR. ITYPE.EQ.13 .OR. ITYPE.EQ.14) THEN
C
C TYPES 0100,1100,1101 HAVE PAIRS IJ AND KL SWITCHED
C
         INEW= KSHELL
         JNEW= LSHELL
         KNEW= ISHELL
         LNEW= JSHELL
         LA= LCT
         LB= LDT
         LC= LAT
         LD= LBT
         IB(1,IEXCH)= 3
         IB(2,IEXCH)= 4
         IB(3,IEXCH)= 1
         IB(4,IEXCH)= 2
      ELSEIF(ITYPE.EQ.11) THEN
C
C TYPE 1010 HAS IJ SWITCHED AND KL SWITCHED
C
         INEW= JSHELL
         JNEW= ISHELL
         KNEW= LSHELL
         LNEW= KSHELL
         LA= LBT
         LB= LAT
         LC= LDT
         LD= LCT
         IB(1,IEXCH)= 2
         IB(2,IEXCH)= 1
         IB(3,IEXCH)= 4
         IB(4,IEXCH)= 3
      ELSEIF(ITYPE.EQ. 9) THEN
C
C TYPE 1000  HAS PAIRS IJ AND KL SWITCHED FOLLOWED BY KL SWITCH
C
         INEW= KSHELL
         JNEW= LSHELL
         KNEW= JSHELL
         LNEW= ISHELL
         LA= LCT
         LB= LDT
         LC= LBT
         LD= LAT
         IB(1,IEXCH)= 4
         IB(2,IEXCH)= 3
         IB(3,IEXCH)= 1
         IB(4,IEXCH)= 2
      ELSEIF(ITYPE.EQ.15) THEN
C
C TYPE 1110 HAS PAIRS IJ AND KL SWITCHED FOLLOWED BY IJ SWITCH
C
         INEW= LSHELL
         JNEW= KSHELL
         KNEW= ISHELL
         LNEW= JSHELL
         LA= LDT
         LB= LCT
         LC= LAT
         LD= LBT
         IB(1,IEXCH)= 3
         IB(2,IEXCH)= 4
         IB(3,IEXCH)= 2
         IB(4,IEXCH)= 1
      ENDIF
C
C ONLY 6 STANDARD TYPES REMAIN. 0000,0001,0011,0101,0111,1111
C SPECIFY THESE BY JTYPE
C
      IF(ITYPE.EQ. 1) THEN
         JTYPE= 1
         NGOUT= 1
      ELSEIF(ITYPE.EQ. 2 .OR. ITYPE.EQ. 3 .OR.
     2       ITYPE.EQ. 5 .OR. ITYPE.EQ. 9) THEN
         JTYPE= 2
         NGOUT= 4
      ELSEIF(ITYPE.EQ. 4 .OR. ITYPE.EQ.13) THEN
         JTYPE= 3
         NGOUT= 16
      ELSEIF(ITYPE.EQ. 6 .OR. ITYPE.EQ. 7 .OR.
     2       ITYPE.EQ.10 .OR. ITYPE.EQ.11) THEN
         JTYPE= 4
         NGOUT= 64
      ELSEIF(ITYPE.EQ. 8 .OR. ITYPE.EQ.12 .OR.
     2       ITYPE.EQ.14 .OR. ITYPE.EQ.15) THEN
         JTYPE= 5
         NGOUT= 64
      ELSEIF(ITYPE.EQ.16) THEN
         JTYPE= 6
         NGOUT= 256
      ENDIF
C
C EMPTY COMMON GOUT
C
      CALL xVCLR(GPOPLE(1+NORG),1,NGOUT)
C
      IF(SCHSKP) GO TO 999
C
C OBTAIN INFORMATION ABOUT SHELLS: INEW, KNEW, JNEW, LNEW
C NUMBER OF GAUSSIANS GO INTO NGA,... IN COMMON SHLLFO
C SHELL ANGULAR QUANTUM NUMBERS LA,... GO INTO COMMON SHLLFO
C GAUSSIAN EXPONENTS GO INTO ARRAYS EXA,EXB,EXC,EXD IN COMMON SHLLFO
C GAUSSIAN COEFFICIENTS GO INTO ARRAYS CSA,CPA,... IN COMMON SHLLFO
C
C NUMBERS OF GAUSSIAN FUNCTIONS IN SHELLS INEW JNEW KNEW AND LNEW
C
      NGA= KNG(INEW)
      NGB= KNG(JNEW)
      NGC= KNG(KNEW)
      NGD= KNG(LNEW)
C
C STARTING LOCATIONS OF SHELLS INEW JNEW KNEW AND LNEW IN LIST
C OF GAUSSIAN FUNCTIONS
C
      I= KSTART(INEW)-1
      J= KSTART(JNEW)-1
      K= KSTART(KNEW)-1
      L= KSTART(LNEW)-1
C
C LOOP OVER GAUSSIANS IN EACH SHELL
C FIRST SHELL INEW
C
      DO 110 NI=1,NGA
         N=I+NI
C
C THE MAXIMUM COEFFICIENT ASSOCIATED WITH SHELL IS
C USED TO DETERMINE IF ANY OF THE INTEGRALS ASSOCIATED WITH A SET
C OF SHELLS IS LARGE ENOUGH TO WARRANT EVALUATION OF THE ENTIRE SET
C
         CMAXA(NI)= CMAX(N)
C
C GAUSSIAN EXPONENTS
C
         EXA(NI)= EX(N)
C
C S COEFFICIENTS
C
         CSA(NI)= CS(N)
C
C P COEFFICIENTS
C
         CPA(NI)= CP(N)
  110 CONTINUE
C
C REPEAT PROCEDURE FOR SHELLS JNEW KNEW AND LNEW
C
      DO 120 NJ=1,NGB
         N=J+NJ
         CMAXB(NJ)= CMAX(N)
         EXB(NJ)= EX(N)
         CSB(NJ)= CS(N)
         CPB(NJ)= CP(N)
  120 CONTINUE
      DO 130 NK=1,NGC
         N=K+NK
         CMAXC(NK)=CMAX(N)
         EXC(NK)=EX(N)
         CSC(NK)=CS(N)*QQ4
         CPC(NK)=CP(N)*QQ4
  130 CONTINUE
      DO 140 NL=1,NGD
         N=L+NL
         CMAXD(NL)=CMAX(N)
         EXD(NL)=EX(N)
         CSD(NL)=CS(N)
         CPD(NL)=CP(N)
  140 CONTINUE
C
      MAB= LA+LB-1
      MCD= LC+LD-1
      NGANGB=NGA*NGB
C
C COORDINATES OF ATOMS ASSOCIATED WITH SHELLS INEW JNEW KNEW AND LNEW
C
      R12= ZER
      R34= ZER
      DO 150 N=1,3
         P12(N,1)= CO(INEW,N)
         P12(N,2)= CO(JNEW,N)
         P12(N,3)= P12(N,2)-P12(N,1)
      R12= R12+P12(N,3)*P12(N,3)
         P34(N,1)= CO(KNEW,N)
         P34(N,2)= CO(LNEW,N)
         P34(N,3)= P34(N,2)-P34(N,1)
  150 R34= R34+P34(N,3)*P34(N,3)
C
C FIND DIRECTION COSINES OF PENULTIMATE AXES FROM COORDINATES OF AB
C P(1,1),P(1,2),... ARE DIRECTION COSINES OF AXES AT P.  Z-AXIS ALONG AB
C T(1),T(2),T(3)... ARE DIRECTION COSINES OF AXES AT Q.  Z-AXIS ALONG CD
C
C FIND DIRECTION COSINES OF AB AND CD. THESE ARE LOCAL Z-AXES.
C IF INDETERMINATE TAKE ALONG SPACE Z-AXIS
C
      P(1,3)= ZER
      P(2,3)= ZER
      P(3,3)= ONE
      RAB= SQRT(R12)
      IF(RAB.NE.ZER) THEN
         TMP= ONE/RAB
         P(1,3)= P12(1,3)*TMP
         P(2,3)= P12(2,3)*TMP
         P(3,3)= P12(3,3)*TMP
      ENDIF
C
      T(1)= ZER
      T(2)= ZER
      T(3)= ONE
      RCD= SQRT(R34)
      IF(RCD.NE.ZER) THEN
         TMP= ONE/RCD
         T(1)= P34(1,3)*TMP
         T(2)= P34(2,3)*TMP
         T(3)= P34(3,3)*TMP
      ENDIF
C
C FIND LOCAL Y-AXIS AS COMMON PERPENDICULAR TO AB AND CD
C IF INDETERMINATE TAKE PERPENDICULAR TO AB AND SPACE Z-AXIS
C IF STILL INDETERMINATE TAKE PERPENDICULAR TO AB AND SPACE X-AXIS
C
      COSG= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
      COSG= MIN( ONE,COSG)
      COSG= MAX(-ONE,COSG)
      SING= SQRT(ONE-COSG*COSG)
C
C MODIFIED ROTATION TESTING.
C THIS FIX CURES THE SMALL ANGLE PROBLEM.
C
      P(1,2)= T(3)*P(2,3)-T(2)*P(3,3)
      P(2,2)= T(1)*P(3,3)-T(3)*P(1,3)
      P(3,2)= T(2)*P(1,3)-T(1)*P(2,3)
      IF( ABS(COSG).GT.PT9)
     1   SING= SQRT(P(1,2)*P(1,2)+P(2,2)*P(2,2)+P(3,2)*P(3,2))
      IF( ABS(COSG).LE.PT9 .OR. SING.GE.TENM12) THEN
         TMP= ONE/SING
         P(1,2)= P(1,2)*TMP
         P(2,2)= P(2,2)*TMP
         P(3,2)= P(3,2)*TMP
      ELSE
         I=3
         IF( ABS(P(1,3)).LE.PT7) I=1
         TMP = P(I,3)*P(I,3)
         TMP = MIN( ONE,TMP)
         TMP = SQRT(ONE-TMP)
         IF(TMP.NE.ZER) TMP= ONE/TMP
         IF( ABS(P(1,3)).LE.PT7) THEN
            P(1,2)= ZER
            P(2,2)= P(3,3)*TMP
            P(3,2)=-P(2,3)*TMP
         ELSE
            P(1,2)= P(2,3)*TMP
            P(2,2)=-P(1,3)*TMP
            P(3,2)= ZER
         ENDIF
      ENDIF
C
C FIND DIRECTION COSINES OF LOCAL X-AXES
C
      P(1,1)= P(2,2)*P(3,3)-P(3,2)*P(2,3)
      P(2,1)= P(3,2)*P(1,3)-P(1,2)*P(3,3)
      P(3,1)= P(1,2)*P(2,3)-P(2,2)*P(1,3)
C
C FIND COORDINATES OF C RELATIVE TO LOCAL AXES AT A
C
      T(1)= P34(1,1)-P12(1,1)
      T(2)= P34(2,1)-P12(2,1)
      T(3)= P34(3,1)-P12(3,1)
      ACX = T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
      ACY = T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
      ACZ = T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
C
C SET ACY= 0  IF CLOSE
C
      IF( ABS(ACY).LE.TENM04) THEN
         ACY =ZER
         ACY2=ZER
      ELSE
         ACY2=ACY*ACY
      ENDIF
C
C DIRECTION COSINES OF CD LOCAL AXES WITH RESPECT TO AB LOCAL AXES
C (COSG,0,-SING)  (0,1,0)  (SING,0,COSG)
C
C PRELIMINARY P LOOP
C
C FILL GEOMPQ WITH INFORMATION ABOUT P IN PRELIMINARY P-LOOP
C
      JI= 1
      DO 170 I=1,NGA
         X01= EXA(I)
         DO 170 J=1,NGB
            X02= EXB(J)
            X12= X01+X02
            X21= ONE/X12
            Y01= X01*X21
            Y02= ONE-Y01
            Y12= Y01*X02
            TX12(JI)= X12
            TX21(JI)= X21
            TY02(JI)= Y02*RAB
CJMS        TY01(JI)=-Y01*RAB
            TY01(JI)= TY02(JI)-RAB
            IF(JTYPE.EQ.4 .OR. JTYPE.EQ.5) TY01(JI)= TY01(JI)*X12
            ISMLP(JI)=2
            R12Y12= R12*Y12
            IF(R12Y12.GT.CUX) GO TO 160
            E12= X21* EXP(-R12Y12)
            TST= E12*CMAXA(I)*CMAXB(J)
            ISMLP(JI)=0
            IF(TST.LE.EI1) ISMLP(JI)=1
            IF(TST.LE.EI2) ISMLP(JI)=2
            E12= PIF*E12
C
C FOR TYPES 0000,0001,0011 ONLY D00P NEEDED
C
            D00P(JI)= E12*CSA(I)*CSB(J)
            IF(JTYPE.GT.3) THEN
               D01P(JI)= E12*CSA(I)*CPB(J)
               IF(JTYPE.LE.5) THEN
                  CONP(JI)= D01P(JI)*X21
                  IF(CONP(JI).NE.ZER) THEN
                     D00P(JI)= D00P(JI)/CONP(JI)
                  ENDIF
               ELSE
                  D10P(JI)= E12*CPA(I)*CSB(J)
                  D11P(JI)= E12*CPA(I)*CPB(J)
                  CONP(JI)= D11P(JI)
                  IF(CONP(JI).NE.ZER) THEN
                     TMP = ONE/CONP(JI)
                     D00P(JI)= D00P(JI)*TMP
                     D01P(JI)= D01P(JI)*TMP
                     D10P(JI)= D10P(JI)*TMP
                  ENDIF
               ENDIF
            ENDIF
  160       CONTINUE
  170 JI=JI+1
C
C IF 0000 INITIALIZE
C
      IF(JTYPE.EQ.1) G(1)= ZER
C
C BEGIN Q LOOP
C
      DO 190 K=1,NGC
         X03= EXC(K)
         DO 190 L=1,NGD
            X04= EXD(L)
            X34= X03+X04
            X43= ONE/X34
            Y03= X03*X43
            Y04= ONE-Y03
            Y34= Y03*X04
            R34Y34= R34*Y34
            IF(R34Y34.GT.CUX) GO TO 180
            E34= X43* EXP(-R34Y34)
            TST= E34*CMAXC(K)*CMAXD(L)
            ISMLQ= 0
            IF(TST.LE.EI1) ISMLQ= 1
            IF(TST.LE.EI2) GO TO 180
C
C RQ(0,0) = COMPONENT OF CQ ALONG PENULTIMATE X-AXIS
C RQ(1,0) = COMPONENT OF CQ ALONG PENULTIMATE Z-AXIS
C RQ(0,1) = COMPONENT OF DQ ALONG PENULTIMATE X-AXIS
C RQ(1,1) = COMPONENT OF DQ ALONG PENULTIMATE Z-AXIS
C
            CQ = RCD*Y04
            DQ =-RCD+CQ
            RQ(0,0)= CQ *SING
            RQ(1,0)= CQ *COSG
            RQ(0,1)= DQ *SING
            RQ(1,1)= DQ *COSG
C
C FIND COORDINATES OF Q RELATIVE TO AXES AT A
C QPR IS PERPENDICULAR FROM Q TO AB
C
            AQX= ACX+RQ(0,0)
            AQZ= ACZ+RQ(1,0)
            QPS= AQX*AQX+ACY2
            QPR= SQRT(QPS)
C
            IF(JTYPE.NE.1) THEN
C
C PHI IS 180 - AZIMUTHAL ANGLE FOR AQ IN AB LOCAL AXIS SYSTEM
C
               COSP= ONE
               SINP= ZER
               IF(QPR.GT.TENM12) THEN
                  TMP = ONE/QPR
                  COSP=-AQX*TMP
                  SINP=-ACY*TMP
               ENDIF
               SQ(0,0)= E34*CSC(K)*CSD(L)
               SQ(1,0)= E34*CSC(K)*CPD(L)
               SQ(0,1)= E34*CPC(K)*CSD(L)
               SQ(1,1)= E34*CPC(K)*CPD(L)
            ENDIF
C
C USE SPECIAL FAST ROUTINE FOR INNER LOOPS FOR 0000 ... 1111
C
            CALL xSP0S1S(JTYPE)
C
            IF(JTYPE.EQ.1) THEN
               SQ(0,0)= E34*CSC(K)*CSD(L)
               G(1)= G(1)+H(1)*SQ(0,0)
            ELSEIF(JTYPE.EQ.2) THEN
               G(2)= H(1)*RQ(0,1)+H(2)*COSP
               G(3)=              H(2)*SINP
               G(4)= H(1)*RQ(1,1)+H(4)
               GPOPLE(1+NORG)= GPOPLE(1+NORG)+H(1)*SQ(0,0)
               GPOPLE(2+NORG)= GPOPLE(2+NORG)+G(2)*SQ(1,0)
               GPOPLE(3+NORG)= GPOPLE(3+NORG)+G(3)*SQ(1,0)
               GPOPLE(4+NORG)= GPOPLE(4+NORG)+G(4)*SQ(1,0)
            ELSE
               CALL xROT2  (JTYPE,MAB,MCD,COSP,SINP)
               CALL xTQ0S1S(JTYPE,R34,GPOPLE(1+NORG))
            ENDIF
  180       CONTINUE
  190 CONTINUE
C
      IF(JTYPE.EQ.1) THEN
         GPOPLE(1+NORG)= G(1)
      ELSE
C
CJMS  NOW, THE TRANSPOSE OF P TO BE USED FOR COMPUTATIONAL EFFICIENCY
C
         DO 195 J=1,2
            DO 195 I=J+1,3
               TMP= P(I,J)
               P(I,J)= P(J,I)
               P(J,I)= TMP
  195    CONTINUE
C
         CALL xR30S1S(JTYPE,GPOPLE(1+NORG),P)
      ENDIF
  999 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------
C*MODULE INT2B   *DECK R30S1S
      SUBROUTINE xR30S1S(JTYPE,F,P)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C ROTATE UP TO 256 INTEGRALS TO SPACE FIXED AXES
C INCOMING AND OUTGOING INTEGRALS IN F
C P(1,1),...ARE DIRECTION COSINES OF SPACE FIXED AXES WRT AXES AT P
C
      DIMENSION F(0:3,0:3,0:3,0:3),P(3,3)
C
      DIMENSION T(3)
C
      IF(JTYPE.EQ.2) THEN
         T(1)= F(1,0,0,0)
         T(2)= F(2,0,0,0)
         T(3)= F(3,0,0,0)
         F(1,0,0,0)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
         F(2,0,0,0)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
         F(3,0,0,0)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
      ELSEIF(JTYPE.EQ.3) THEN
         DO 310 I=0,3
            T(1)= F(I,1,0,0)
            T(2)= F(I,2,0,0)
            T(3)= F(I,3,0,0)
            F(I,1,0,0)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
            F(I,2,0,0)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
            F(I,3,0,0)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
  310    CONTINUE
         DO 320 J=0,3
            T(1)= F(1,J,0,0)
            T(2)= F(2,J,0,0)
            T(3)= F(3,J,0,0)
            F(1,J,0,0)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
            F(2,J,0,0)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
            F(3,J,0,0)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
  320    CONTINUE
      ELSEIF(JTYPE.EQ.4) THEN
         DO 410 I=0,3
            T(1)= F(I,0,1,0)
            T(2)= F(I,0,2,0)
            T(3)= F(I,0,3,0)
            F(I,0,1,0)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
            F(I,0,2,0)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
            F(I,0,3,0)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
  410    CONTINUE
         DO 420 K=0,3
            T(1)= F(1,0,K,0)
            T(2)= F(2,0,K,0)
            T(3)= F(3,0,K,0)
            F(1,0,K,0)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
            F(2,0,K,0)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
            F(3,0,K,0)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
  420    CONTINUE
      ELSEIF(JTYPE.EQ.5) THEN
         DO 510 J=0,3
            DO 510 I=0,3
               T(1)= F(I,J,1,0)
               T(2)= F(I,J,2,0)
               T(3)= F(I,J,3,0)
               F(I,J,1,0)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
               F(I,J,2,0)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
               F(I,J,3,0)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
  510    CONTINUE
         DO 520 K=0,3
            DO 520 I=0,3
               T(1)= F(I,1,K,0)
               T(2)= F(I,2,K,0)
               T(3)= F(I,3,K,0)
               F(I,1,K,0)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
               F(I,2,K,0)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
               F(I,3,K,0)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
  520    CONTINUE
         DO 530 K=0,3
            DO 530 J=0,3
               T(1)= F(1,J,K,0)
               T(2)= F(2,J,K,0)
               T(3)= F(3,J,K,0)
               F(1,J,K,0)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
               F(2,J,K,0)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
               F(3,J,K,0)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
  530    CONTINUE
      ELSEIF(JTYPE.EQ.6) THEN
         DO 610 K=0,3
            DO 610 J=0,3
               DO 610 I=0,3
                  T(1)= F(I,J,K,1)
                  T(2)= F(I,J,K,2)
                  T(3)= F(I,J,K,3)
                  F(I,J,K,1)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
                  F(I,J,K,2)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
                  F(I,J,K,3)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
  610    CONTINUE
         DO 620 L=0,3
            DO 620 J=0,3
               DO 620 I=0,3
                  T(1)= F(I,J,1,L)
                  T(2)= F(I,J,2,L)
                  T(3)= F(I,J,3,L)
                  F(I,J,1,L)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
                  F(I,J,2,L)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
                  F(I,J,3,L)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
  620    CONTINUE
         DO 630 L=0,3
            DO 630 K=0,3
               DO 630 I=0,3
                  T(1)= F(I,1,K,L)
                  T(2)= F(I,2,K,L)
                  T(3)= F(I,3,K,L)
                  F(I,1,K,L)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
                  F(I,2,K,L)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
                  F(I,3,K,L)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
  630    CONTINUE
         DO 640 L=0,3
            DO 640 K=0,3
               DO 640 J=0,3
                  T(1)= F(1,J,K,L)
                  T(2)= F(2,J,K,L)
                  T(3)= F(3,J,K,L)
                  F(1,J,K,L)= T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
                  F(2,J,K,L)= T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
                  F(3,J,K,L)= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
  640    CONTINUE
      ENDIF
      RETURN
      END
c-----------------------------------------------------------------
C*MODULE INT2B   *DECK ROT2
      SUBROUTINE xROT2  (JTYPE,MAB,MCD,COSP,SINP)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /xJMSGYH/ G(0:3,0:3,0:3,0:3),H(0:3,0:3,0:3,0:3),
     2                RQ(2,2),SQ(2,2)
C
      DIMENSION RP(2,2)
C
      PARAMETER (ZER=0.0D+00, ONE=1.0D+00)
C
C NO ROTATION IF ABCD COPLANAR
C
      K2=0
      IF(MAB.GE.0) K2=3
      L2=0
      IF(MAB.GT.0) L2=3
      G(0,0,0,0)= H(0,0,0,0)
      IF(COSP.NE.ZER .AND. SINP.EQ.ZER) THEN
         G(1,0,0,0)= H(1,0,0,0)
         G(2,0,0,0)= ZER
         G(3,0,0,0)= H(3,0,0,0)
         IF(JTYPE.LT.4) THEN
            G(1,1,0,0)= H(1,1,0,0)
            G(2,1,0,0)= ZER
            G(3,1,0,0)= H(3,1,0,0)
            G(2,2,0,0)= H(2,2,0,0)
            G(3,2,0,0)= ZER
            G(3,3,0,0)= H(3,3,0,0)
            IF(COSP.LT.ZER) THEN
               G(1,0,0,0)=-G(1,0,0,0)
               G(3,1,0,0)=-G(3,1,0,0)
            ENDIF
            GO TO 999
         ELSEIF(JTYPE.LT.5) THEN
            G(0,0,1,0)= H(0,0,1,0)
            G(1,0,1,0)= H(1,0,1,0)
            G(2,0,1,0)= ZER
            G(3,0,1,0)= H(3,0,1,0)
            G(0,0,2,0)= ZER
            G(1,0,2,0)= ZER
            G(2,0,2,0)= H(2,0,2,0)
            G(3,0,2,0)= ZER
            G(0,0,3,0)= H(0,0,3,0)
            G(1,0,3,0)= H(1,0,3,0)
            G(2,0,3,0)= ZER
            G(3,0,3,0)= H(3,0,3,0)
            IF(COSP.LT.ZER) THEN
               G(1,0,0,0)=-G(1,0,0,0)
               G(0,0,1,0)=-G(0,0,1,0)
               G(3,0,1,0)=-G(3,0,1,0)
               G(1,0,3,0)=-G(1,0,3,0)
            ENDIF
            GO TO 999
         ENDIF
         DO 110 L=0,L2
            DO 110 K=0,K2
               IF((K.NE.2.AND.L.NE.2) .OR.
     1            (K.EQ.2.AND.L.EQ.2)) THEN
                  G(0,0,K,L)= H(0,0,K,L)
                  G(1,0,K,L)= H(1,0,K,L)
                  G(2,0,K,L)= ZER
                  G(3,0,K,L)= H(3,0,K,L)
                  G(1,1,K,L)= H(1,1,K,L)
                  G(2,1,K,L)= ZER
                  G(3,1,K,L)= H(3,1,K,L)
                  G(2,2,K,L)= H(2,2,K,L)
                  G(3,2,K,L)= ZER
                  G(3,3,K,L)= H(3,3,K,L)
               ELSE
                  G(0,0,K,L)= ZER
                  G(1,0,K,L)= ZER
                  G(2,0,K,L)= H(2,0,K,L)
                  G(3,0,K,L)= ZER
                  G(1,1,K,L)= ZER
                  G(2,1,K,L)= H(2,1,K,L)
                  G(3,1,K,L)= ZER
                  G(2,2,K,L)= ZER
                  G(3,2,K,L)= H(3,2,K,L)
                  G(3,3,K,L)= ZER
               ENDIF
  110    CONTINUE
         IF(COSP.LT.ZER) THEN
            G(1,0,0,0)=-G(1,0,0,0)
            G(3,1,0,0)=-G(3,1,0,0)
            G(0,0,1,0)=-G(0,0,1,0)
            G(3,0,1,0)=-G(3,0,1,0)
            G(1,0,3,0)=-G(1,0,3,0)
            G(1,1,1,0)=-G(1,1,1,0)
            G(2,2,1,0)=-G(2,2,1,0)
            G(3,3,1,0)=-G(3,3,1,0)
            G(2,1,2,0)=-G(2,1,2,0)
            G(3,1,3,0)=-G(3,1,3,0)
            IF(JTYPE.LT.6) GO TO 999
            G(0,0,0,1)=-G(0,0,0,1)
            G(3,0,0,1)=-G(3,0,0,1)
            G(1,1,0,1)=-G(1,1,0,1)
            G(2,2,0,1)=-G(2,2,0,1)
            G(3,3,0,1)=-G(3,3,0,1)
            G(1,0,1,1)=-G(1,0,1,1)
            G(3,1,1,1)=-G(3,1,1,1)
            G(2,0,2,1)=-G(2,0,2,1)
            G(3,2,2,1)=-G(3,2,2,1)
            G(0,0,3,1)=-G(0,0,3,1)
            G(3,0,3,1)=-G(3,0,3,1)
            G(1,1,3,1)=-G(1,1,3,1)
            G(2,2,3,1)=-G(2,2,3,1)
            G(3,3,3,1)=-G(3,3,3,1)
            G(2,1,0,2)=-G(2,1,0,2)
            G(2,0,1,2)=-G(2,0,1,2)
            G(3,2,1,2)=-G(3,2,1,2)
            G(1,0,2,2)=-G(1,0,2,2)
            G(3,1,2,2)=-G(3,1,2,2)
            G(2,1,3,2)=-G(2,1,3,2)
            G(1,0,0,3)=-G(1,0,0,3)
            G(3,1,0,3)=-G(3,1,0,3)
            G(0,0,1,3)=-G(0,0,1,3)
            G(3,0,1,3)=-G(3,0,1,3)
            G(1,1,1,3)=-G(1,1,1,3)
            G(2,2,1,3)=-G(2,2,1,3)
            G(3,3,1,3)=-G(3,3,1,3)
            G(2,1,2,3)=-G(2,1,2,3)
            G(1,0,3,3)=-G(1,0,3,3)
            G(3,1,3,3)=-G(3,1,3,3)
         ENDIF
      ELSE
         I2=0
         IF(MCD.GE.0) I2=3
         J2=0
         IF(MCD.GT.0) J2=3
         RP(1,1)= COSP
         RP(2,1)= SINP
         RP(1,2)=-SINP
         RP(2,2)= COSP
         IF(I2.EQ.3) THEN
            IF(J2.EQ.3) THEN
               V44= RP(1,1)*RP(1,1)
               V47= ONE-V44
               V54= RP(1,1)*RP(2,1)
               V45=-V54-V54
               V55= V44-V47
               IF(JTYPE.LT.4) THEN
                  G(1,0,0,0)= H(1,0,0,0)*RP(1,1)
                  G(2,0,0,0)= H(1,0,0,0)*RP(2,1)
                  G(3,0,0,0)= H(3,0,0,0)
                  TMP1= H(1,1,0,0)-H(2,2,0,0)
                  TMP2= TMP1*V47
                  G(1,1,0,0)= H(1,1,0,0)-TMP2
                  G(2,1,0,0)= TMP1*V54
                  G(3,1,0,0)= H(3,1,0,0)*RP(1,1)
                  G(2,2,0,0)= H(2,2,0,0)+TMP2
                  G(3,2,0,0)= H(3,1,0,0)*RP(2,1)
                  G(3,3,0,0)= H(3,3,0,0)
                  GO TO 999
               ENDIF
               DO 120 L=0,L2
                  DO 120 K=0,K2
                     IF((K.NE.2.AND.L.NE.2) .OR.
     1                  (K.EQ.2.AND.L.EQ.2)) THEN
                        TMP1= H(1,1,K,L)-H(2,2,K,L)
                        TMP2= TMP1*V47
                        G(1,1,K,L)= H(1,1,K,L)-TMP2
                        G(2,1,K,L)= TMP1*V54
                        G(3,1,K,L)= H(3,1,K,L)*RP(1,1)
                        G(2,2,K,L)= H(2,2,K,L)+TMP2
                        G(3,2,K,L)= H(3,1,K,L)*RP(2,1)
                        G(3,3,K,L)= H(3,3,K,L)
                     ELSE
                        TMP = H(2,1,K,L)*V45
                        G(1,1,K,L)= TMP
                        G(2,1,K,L)= H(2,1,K,L)*V55
                        G(3,1,K,L)= H(3,2,K,L)*RP(1,2)
                        G(2,2,K,L)=-TMP
                        G(3,2,K,L)= H(3,2,K,L)*RP(2,2)
                        G(3,3,K,L)= ZER
                     ENDIF
  120          CONTINUE
            ENDIF
            IF(JTYPE.LT.6) THEN
               DO 125 K=0,K2
                  N=1
                  IF(K.EQ.2) N=3-N
                  G(1,0,K,0)= H(N,0,K,0)*RP(1,N)
                  G(2,0,K,0)= H(N,0,K,0)*RP(2,N)
                  G(3,0,K,0)= H(3,0,K,0)
  125          CONTINUE
            ELSE
               DO 130 L=0,L2
                  M=1
                  IF(L.EQ.2) M=3-M
                  DO 130 K=0,K2
                     N=M
                     IF(K.EQ.2) N=3-N
                     G(1,0,K,L)= H(N,0,K,L)*RP(1,N)
                     G(2,0,K,L)= H(N,0,K,L)*RP(2,N)
                     G(3,0,K,L)= H(3,0,K,L)
  130          CONTINUE
               G(3,0,2,1)= ZER
               G(3,0,0,2)= ZER
               G(3,0,1,2)= ZER
               G(3,0,3,2)= ZER
               G(3,0,2,3)= ZER
            ENDIF
         ENDIF
         IF(K2.EQ.3) THEN
            G(3,0,2,0)= ZER
            G(0,0,1,0)= H(0,0,1,0)*RP(1,1)
            G(0,0,2,0)= H(0,0,1,0)*RP(2,1)
            G(0,0,3,0)= H(0,0,3,0)
            DO 140 I=1,I2
               T01= G(I,0,1,0)
               T02= G(I,0,2,0)
               G(I,0,1,0)= T01*RP(1,1)+T02*RP(1,2)
               G(I,0,2,0)= T01*RP(2,1)+T02*RP(2,2)
  140       CONTINUE
            IF(JTYPE.LT.5) GO TO 999
            DO 150 J=1,J2
               DO 150 I=J,I2
                  T01= G(I,J,1,0)
                  T02= G(I,J,2,0)
                  G(I,J,1,0)= T01*RP(1,1)+T02*RP(1,2)
                  G(I,J,2,0)= T01*RP(2,1)+T02*RP(2,2)
  150       CONTINUE
            IF(JTYPE.LT.6) GO TO 999
            IF(L2.EQ.3) THEN
               G(0,0,0,1)= H(0,0,0,1)*RP(1,1)
               G(0,0,0,2)= H(0,0,0,1)*RP(2,1)
               G(0,0,0,3)= H(0,0,0,3)
               DO 160 I=1,I2
                  T01= G(I,0,0,1)
                  T02= G(I,0,0,2)
                  G(I,0,0,1)= T01*RP(1,1)+T02*RP(1,2)
                  G(I,0,0,2)= T01*RP(2,1)+T02*RP(2,2)
  160          CONTINUE
               G(0,0,1,1)= H(0,0,1,1)
               G(0,0,2,1)= ZER
               G(0,0,3,1)= H(0,0,3,1)
               G(0,0,1,2)= ZER
               G(0,0,2,2)= H(0,0,2,2)
               G(0,0,3,2)= ZER
               G(0,0,1,3)= H(0,0,1,3)
               G(0,0,2,3)= ZER
               G(0,0,3,3)= H(0,0,3,3)
               DO 170 J=1,J2
                  DO 170 I=J,I2
                     T01= G(I,J,0,1)
                     T02= G(I,J,0,2)
                     G(I,J,0,1)= T01*RP(1,1)+T02*RP(1,2)
                     G(I,J,0,2)= T01*RP(2,1)+T02*RP(2,2)
  170          CONTINUE
               DO 180 K=1,K2
                  DO 180 J=0,J2
                     DO 180 I=J,I2
                        T01= G(I,J,K,1)
                        T02= G(I,J,K,2)
                        H(I,J,K,1)= T01*RP(1,1)+T02*RP(1,2)
                        H(I,J,K,2)= T01*RP(2,1)+T02*RP(2,2)
                        H(I,J,K,3)= G(I,J,K,3)
  180          CONTINUE
               DO 190 L=1,L2
                  DO 190 J=0,J2
                     DO 190 I=J,I2
                        T01= H(I,J,1,L)
                        T02= H(I,J,2,L)
                        G(I,J,1,L)= T01*RP(1,1)+T02*RP(1,2)
                        G(I,J,2,L)= T01*RP(2,1)+T02*RP(2,2)
                        G(I,J,3,L)= H(I,J,3,L)
  190          CONTINUE
            ENDIF
         ENDIF
      ENDIF
  999 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------
C*MODULE INT2B   *DECK TQ0S1S
      SUBROUTINE xTQ0S1S(JTYPE,R34,F)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION F(0:3,0:3,0:3,0:3)
C
C TRANSLATE UP TO 160 INTEGRALS ON A B AND Q
C        TO UP TO 256 INTEGRALS ON A B C AND D
C
      COMMON /xJMSGYH/ G(0:3,0:3,0:3,0:3),H(256),
     2                RQ(0:1,0:1),SQ(0:1,0:1)
C
      PARAMETER (ZER=0.0D+00)
C
      IF(JTYPE.EQ.3) THEN
         G(0,1,0,0)= G(1,0,0,0)
         G(0,2,0,0)= G(2,0,0,0)
         G(0,3,0,0)= G(3,0,0,0)
         G(1,2,0,0)= G(2,1,0,0)
         G(1,3,0,0)= G(3,1,0,0)
         G(2,3,0,0)= G(3,2,0,0)
         IF(R34.GT.ZER) THEN
            G(0,1,0,0)= G(0,1,0,0)+G(0,0,0,0)*RQ(0,0)
            G(0,3,0,0)= G(0,3,0,0)+G(0,0,0,0)*RQ(1,0)
            DO 310 I=1,3
               G(I,1,0,0)= G(I,1,0,0)+G(I,0,0,0)*RQ(0,0)
               G(I,3,0,0)= G(I,3,0,0)+G(I,0,0,0)*RQ(1,0)
               G(1,I,0,0)= G(1,I,0,0)+G(0,I,0,0)*RQ(0,1)
               G(3,I,0,0)= G(3,I,0,0)+G(0,I,0,0)*RQ(1,1)
  310       CONTINUE
            G(1,0,0,0)= G(1,0,0,0)+G(0,0,0,0)*RQ(0,1)
            G(3,0,0,0)= G(3,0,0,0)+G(0,0,0,0)*RQ(1,1)
         ENDIF
         N= 0
         DO 320 J=0,3
            F(0,J,0,0)= F(0,J,0,0)+G(0,J,0,0)*SQ(0,N)
            F(1,J,0,0)= F(1,J,0,0)+G(1,J,0,0)*SQ(1,N)
            F(2,J,0,0)= F(2,J,0,0)+G(2,J,0,0)*SQ(1,N)
            F(3,J,0,0)= F(3,J,0,0)+G(3,J,0,0)*SQ(1,N)
  320    N= 1
      ELSEIF(JTYPE.EQ.4) THEN
         DO 410 K=0,3
            IF(R34.GT.ZER) THEN
               G(1,0,K,0)= G(1,0,K,0)+G(0,0,K,0)*RQ(0,1)
               G(3,0,K,0)= G(3,0,K,0)+G(0,0,K,0)*RQ(1,1)
            ENDIF
            F(0,0,K,0)= F(0,0,K,0)+G(0,0,K,0)*SQ(0,0)
            F(1,0,K,0)= F(1,0,K,0)+G(1,0,K,0)*SQ(1,0)
            F(2,0,K,0)= F(2,0,K,0)+G(2,0,K,0)*SQ(1,0)
            F(3,0,K,0)= F(3,0,K,0)+G(3,0,K,0)*SQ(1,0)
  410    CONTINUE
      ELSEIF(JTYPE.EQ.5) THEN
         DO 530 K=0,3
            G(0,1,K,0)= G(1,0,K,0)
            G(0,2,K,0)= G(2,0,K,0)
            G(0,3,K,0)= G(3,0,K,0)
            G(1,2,K,0)= G(2,1,K,0)
            G(1,3,K,0)= G(3,1,K,0)
            G(2,3,K,0)= G(3,2,K,0)
            IF(R34.GT.ZER) THEN
               G(0,1,K,0)= G(0,1,K,0)+G(0,0,K,0)*RQ(0,0)
               G(0,3,K,0)= G(0,3,K,0)+G(0,0,K,0)*RQ(1,0)
               DO 510 I=1,3
                  G(I,1,K,0)= G(I,1,K,0)+G(I,0,K,0)*RQ(0,0)
                  G(I,3,K,0)= G(I,3,K,0)+G(I,0,K,0)*RQ(1,0)
                  G(1,I,K,0)= G(1,I,K,0)+G(0,I,K,0)*RQ(0,1)
                  G(3,I,K,0)= G(3,I,K,0)+G(0,I,K,0)*RQ(1,1)
  510          CONTINUE
               G(1,0,K,0)= G(1,0,K,0)+G(0,0,K,0)*RQ(0,1)
               G(3,0,K,0)= G(3,0,K,0)+G(0,0,K,0)*RQ(1,1)
            ENDIF
            N= 0
            DO 520 J=0,3
               F(0,J,K,0)= F(0,J,K,0)+G(0,J,K,0)*SQ(0,N)
               F(1,J,K,0)= F(1,J,K,0)+G(1,J,K,0)*SQ(1,N)
               F(2,J,K,0)= F(2,J,K,0)+G(2,J,K,0)*SQ(1,N)
               F(3,J,K,0)= F(3,J,K,0)+G(3,J,K,0)*SQ(1,N)
  520       N= 1
  530    CONTINUE
      ELSEIF(JTYPE.EQ.6) THEN
         DO 630 L=0,3
            DO 630 K=0,3
               G(0,1,K,L)= G(1,0,K,L)
               G(0,2,K,L)= G(2,0,K,L)
               G(0,3,K,L)= G(3,0,K,L)
               G(1,2,K,L)= G(2,1,K,L)
               G(1,3,K,L)= G(3,1,K,L)
               G(2,3,K,L)= G(3,2,K,L)
               IF(R34.GT.ZER) THEN
                  G(0,1,K,L)= G(0,1,K,L)+G(0,0,K,L)*RQ(0,0)
                  G(0,3,K,L)= G(0,3,K,L)+G(0,0,K,L)*RQ(1,0)
                  DO 610 I=1,3
                     G(I,1,K,L)= G(I,1,K,L)+G(I,0,K,L)*RQ(0,0)
                     G(I,3,K,L)= G(I,3,K,L)+G(I,0,K,L)*RQ(1,0)
                     G(1,I,K,L)= G(1,I,K,L)+G(0,I,K,L)*RQ(0,1)
                     G(3,I,K,L)= G(3,I,K,L)+G(0,I,K,L)*RQ(1,1)
  610             CONTINUE
                  G(1,0,K,L)= G(1,0,K,L)+G(0,0,K,L)*RQ(0,1)
                  G(3,0,K,L)= G(3,0,K,L)+G(0,0,K,L)*RQ(1,1)
               ENDIF
               N= 0
               DO 620 J=0,3
                  F(0,J,K,L)= F(0,J,K,L)+G(0,J,K,L)*SQ(0,N)
                  F(1,J,K,L)= F(1,J,K,L)+G(1,J,K,L)*SQ(1,N)
                  F(2,J,K,L)= F(2,J,K,L)+G(2,J,K,L)*SQ(1,N)
                  F(3,J,K,L)= F(3,J,K,L)+G(3,J,K,L)*SQ(1,N)
  620          N= 1
  630    CONTINUE
      ENDIF
      RETURN
      END
c-----------------------------------------------------------------
C*MODULE INT2B   *DECK SP0S1S
      SUBROUTINE xSP0S1S(JTYPE)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "gamess.par"
C
      COMMON /xGEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),
     4                CONP(MXG2),NGANGB
      COMMON /xMAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
C
C     ==================================================================
C     NOTE THAT GAMGEN IS THE PLACE THAT SETS COMMON BLOCK /INTACC/
C     THIS GROUP OF STATEMENTS IS INTERRELATED FOR ROUTINES
C     GAMGEN, GENR70 AND SP0S1S IN INT2B, AND JKDG80, GRDG80 AND GRDPER
C     IN GRD2B. IF YOU CHANGE THESE, CHANGE ALL SIX PLACES EQUIVALENTLY.
C     GAMGEN HAS COMMENTS EXPLAINING THIS IN THE MOST DETAIL.
C
      PARAMETER (M3R=401)
      COMMON /xINTACC/ EI1,EI2,CUF,RFI,CUG,RGI,CUX,R(3,0:5,0:M3R)
C     ==================================================================
C
      COMMON /xJMSGYH/ G(256),H(0:3,0:3,0:3,0:3),RQ(2,2),SQ(2,2)
C
      DIMENSION U(0:5),XOV(0:5,2),YOC(0:5,2),ZOW(0:8,2)
      DIMENSION AOB(0:8,2),S(0:13),T(0:13)
      DIMENSION FQ(0:5),FQP(0:4,0:4),FQX(0:2,0:2),HQX(0:2,0:2)
      DIMENSION ARRY1(4,2),ARRY2(4,2),ARRY3(4,2)
C
      PARAMETER (ZER=0.0D+00, PT5=0.5D+00, ONE=1.0D+00, TWO=2.0D+00)
CJMS  PARAMETER (PIDIV4=0.7853981625D+00)
CC    PARAMETER (PIDIV4=0.785398163397448309615660845819D+00)
      PARAMETER (PIDIV4=0.785398163397448D+00)
C
CC    PIDIV4= ATAN(ONE)
      PI4= PIDIV4
C
      MFQ= JTYPE-2
      IF(JTYPE.LE.3) MFQ= JTYPE-1
C
      H(0,0,0,0)= ZER
      IF(JTYPE.EQ.2) THEN
         H(1,0,0,0)= ZER
         H(3,0,0,0)= ZER
      ELSEIF(JTYPE.EQ.3) THEN
         H(1,0,0,0)= ZER
         H(3,0,0,0)= ZER
         H(1,1,0,0)= ZER
         H(3,1,0,0)= ZER
         H(3,3,0,0)= ZER
      ELSEIF(JTYPE.EQ.4) THEN
         H(1,0,0,0)= ZER
         H(3,0,0,0)= ZER
         H(0,0,1,0)= ZER
         H(1,0,1,0)= ZER
         H(3,0,1,0)= ZER
         H(0,0,3,0)= ZER
         H(1,0,3,0)= ZER
         H(3,0,3,0)= ZER
      ELSEIF(JTYPE.GT.4) THEN
         DO 170 J=1,JTYPE-4
            DO 150 I=0,5
               XOV(I,J)= ZER
               YOC(I,J)= ZER
  150       CONTINUE
            DO 160 I=0,8
  160       ZOW(I,J)= ZER
  170    CONTINUE
         IF(JTYPE.EQ.6) THEN
            DO 180 I=0,13
               S(I)= ZER
               T(I)= ZER
  180       CONTINUE
         ENDIF
      ENDIF
C
C BEGIN P LOOP
C
      DO 300 I=1,NGANGB
         ISML= ISMLQ+ISMLP(I)
         IF(ISML.GT.1) GO TO 300
         X12= TX12(I)
         X21= TX21(I)
         Y01= TY01(I)
         Y02= TY02(I)
         X14= X12+X34
         X41= ONE/X14
         PQR= AQZ-Y02
         PQS= PQR*PQR
         RHO= X12*X34*X41
         XVA=(PQS+QPS)*RHO
         IF(XVA.LE.CUF) THEN
            TMP= XVA*RFI
            N  = INT(TMP)
            C02= TMP-N
            TMP= C02*(ONE-C02)
            C03=-TMP*(TWO-C02)
            C04=-TMP*(ONE+C02)
            DO 210 M=0,MFQ
               FQ(M)= R(1,M,N)+C02*R(2,M,N)+C03*R(3,M,N)+
     1                C04*R(3,M,N+1)
  210       CONTINUE
         ELSE
            XIN= ONE/XVA
            FQ(0)= SQRT(PI4*XIN)
            TMP= PT5*XIN
            DO 220 M=1,MFQ
               FQ(M)= FQ(M-1)*TMP
  220       TMP= TMP+XIN
         ENDIF
C
         IF(JTYPE.EQ.3) RHO= RHO*X43
         FQ0= D00P(I)
         IF(JTYPE.GT.3) FQ0= CONP(I)
         TMP= FQ0* SQRT(X41)
         DO 225 M=0,MFQ
            FQP(M,0)=FQ(M)*TMP
  225    TMP= TMP*RHO
C
         IF(JTYPE.LE.3) THEN
            H(0,0,0,0)= H(0,0,0,0)+FQP(0,0)
            IF(JTYPE.GT.1) THEN
               H(1,0,0,0)= H(1,0,0,0)+FQP(1,0)
               H(3,0,0,0)= H(3,0,0,0)-FQP(1,0)*PQR
               IF(JTYPE.GT.2) THEN
                  H(1,1,0,0)= H(1,1,0,0)+FQP(2,0)
                  H(3,1,0,0)= H(3,1,0,0)-FQP(2,0)*PQR
                  H(3,3,0,0)= H(3,3,0,0)+FQP(2,0)*PQS
               ENDIF
            ENDIF
         ELSE
            DP00= D00P(I)
            M=1
            DO 235 L=1,MFQ
               DO 230 K=M,MFQ
  230          FQP(K,L)= FQP(K,L-1)*PQR
  235       M=L
            IF(JTYPE.EQ.4) THEN
               H(0,0,0,0)= H(0,0,0,0)+FQP(0,0)*DP00
               H(1,0,0,0)= H(1,0,0,0)+FQP(1,0)*DP00
               H(3,0,0,0)= H(3,0,0,0)-FQP(1,1)*DP00
               H(0,0,1,0)= H(0,0,1,0)-FQP(1,0)
               H(1,0,1,0)= H(1,0,1,0)-FQP(2,0)
               H(3,0,1,0)= H(3,0,1,0)+FQP(2,1)
               H(0,0,3,0)= H(0,0,3,0)+FQP(0,0)*Y01+FQP(1,1)
               H(1,0,3,0)= H(1,0,3,0)+FQP(1,0)*Y01
               H(3,0,3,0)= H(3,0,3,0)-FQP(1,1)*Y01-FQP(2,2)
            ELSE
               U(0)= FQP(0,0)
               U(1)= FQP(1,0)
               U(2)= FQP(2,0)
               U(3)= FQP(1,1)
               U(4)= FQP(2,1)
               U(5)= FQP(2,2)
               N=1
               M=0
               DO 245 L=1,4
                  DO 240 K=N,3
                     AOB(M,1)= FQP(K,L-1)
  240             M=M+1
  245          N=L
               IF(JTYPE.EQ.5) THEN
                  DO 250 N=0,5
                     XOV(N,1)= XOV(N,1)+U(N)*DP00
                     YOC(N,1)= YOC(N,1)+U(N)*Y01
  250             CONTINUE
                  DO 255 N=0,8
  255             ZOW(N,1)= ZOW(N,1)+AOB(N,1)
               ELSEIF(JTYPE.EQ.6) THEN
                  DO 260 N=0,5
  260             XOV(N,1)= XOV(N,1)+U(N)*DP00
                  DP01= D01P(I)
                  DP10= D10P(I)
                  IF(R12.NE.ZER) THEN
                     TVN= Y02*DP10
                     TYN= Y01*DP01
                     TCN= Y01*Y02
                     DO 265 N=0,5
                        XOV(N,2)= XOV(N,2)+U(N)*TVN
                        YOC(N,1)= YOC(N,1)+U(N)*TYN
                        YOC(N,2)= YOC(N,2)+U(N)*TCN
  265                CONTINUE
                  ENDIF
                  TZN= X21*DP01
                  TWN= X21*DP10
                  DO 270 N=0,8
                     ZOW(N,1)= ZOW(N,1)+AOB(N,1)*TZN
                     ZOW(N,2)= ZOW(N,2)+AOB(N,1)*TWN
  270             CONTINUE
                  X2S= X21*X21
                  N=1
                  M=0
                  DO 285 L=1,4
                     DO 280 K=N,3
                        AOB(M,1)= FQP(K  ,L)*X21
                        AOB(M,2)= FQP(K+1,L)*X2S
  280                M=M+1
  285             N=L
                  DO 290 N=0,3
                     S(N)= S(N)+FQP(N,0)*X21
                     T(N)= T(N)+FQP(N,0)*X2S
  290             CONTINUE
                  T(4)= T(4)+FQP(4,0)*X2S
                  DO 295 N=5,13
                     S(N)= S(N)+AOB(N-5,1)
                     T(N)= T(N)+AOB(N-5,2)
  295             CONTINUE
               ENDIF
            ENDIF
         ENDIF
  300 CONTINUE
C
      IF(JTYPE.EQ.1) GO TO 999
C
      FQX(1,0)= QPR
      FQX(0,1)= X43
      IF(JTYPE.GT.2) H43= PT5*X43
      IF(JTYPE.EQ.2) THEN
         H(1,0,0,0)= H(1,0,0,0)*X43*QPR
         H(3,0,0,0)= H(3,0,0,0)*X43
      ELSEIF(JTYPE.EQ.3) THEN
         TMP   =(H(0,0,0,0)-H(1,0,0,0))*H43
         H(1,0,0,0)= H(1,0,0,0)*QPR
         H(1,1,0,0)= H(1,1,0,0)*QPS+TMP
         H(3,1,0,0)= H(3,1,0,0)*QPR
         H(2,2,0,0)= TMP
         H(3,3,0,0)= H(3,3,0,0)+TMP
      ELSEIF(JTYPE.EQ.4) THEN
         FQX(1,1)= FQX(1,0)*FQX(0,1)
         FQX(2,1)= FQX(1,0)*FQX(1,1)
         TEMP0 =-H(0,0,1,0)*H43
         TEMP1 = H(3,0,1,0)*FQX(1,1)
         H(1,0,0,0)= H(1,0,0,0)*FQX(1,1)
         H(3,0,0,0)= H(3,0,0,0)*FQX(0,1)
         H(0,0,1,0)= H(0,0,1,0)*FQX(1,0)
         H(1,0,1,0)= H(1,0,1,0)*FQX(2,1)+TEMP0
         H(3,0,1,0)= TEMP1
         H(2,0,2,0)= TEMP0
         H(1,0,3,0)= H(1,0,3,0)*FQX(1,1)+TEMP1
         H(3,0,3,0)= H(3,0,3,0)*FQX(0,1)+TEMP0
      ELSE
         FQX(2,0)= QPS
         HQX(0,1)= H43
         FQX(1,1)= FQX(1,0)*FQX(0,1)
         HQX(1,1)= PT5*FQX(1,1)
         FQX(2,1)= FQX(2,0)*FQX(0,1)
         HQX(2,1)= PT5*FQX(2,1)
         FQX(0,2)= FQX(0,1)*FQX(0,1)
         HQX(0,2)= PT5*FQX(0,2)
         FQX(1,2)= FQX(1,0)*FQX(0,2)
         HQX(1,2)= PT5*FQX(1,2)
         FQX(2,2)= FQX(2,0)*FQX(0,2)
         HQX(2,2)= PT5*FQX(2,2)
         DO 310 K=0,5
  310    YOC(K,1)= YOC(K,1)+ZOW(3+K,1)
         YOC(5,1)= YOC(5,1)-ZOW(4  ,1)
         IF(JTYPE.EQ.6) THEN
            DO 320 K=0,5
  320       XOV(K,2)= XOV(K,2)+ZOW(3+K,2)
            XOV(5,2)= XOV(5,2)-ZOW(4  ,2)
            BQZ= AQZ-RAB
            N=0
            M=1
            DO 340 L=1,3
               DO 330 K=N,N+2
                  AOB(K,1)= AQZ*S(M+K)-(S(5+K)-T(5+K))
                  AOB(K,2)= BQZ*S(M+K)-(S(5+K)-T(5+K))
  330          CONTINUE
            M=M+1
  340       N=N+3
            AOB(8,1)= AQZ*S(12)-(S(13)-T(13))
            AOB(8,2)= BQZ*S(12)-(S(13)-T(13))
            DO 350 K=0,5
  350       YOC(K,2)= YOC(K,2)+AOB(3+K,1)+AOB(3+K,2)-T(8+K)
            U(0)  = PT5*(S(0)-T(1))
            U(1)  = PT5*(S(1)-T(2))
            U(2)  = PT5*(S(2)-T(3))
            TEMP0 = H43*(U(0)-X43*U(1))
            TEMP1 = T(2)-QPS*T(3)
            TEMP2 = FQX(2,2)*U(2)+TEMP0
            T2211 = HQX(2,1)*(T(2)-X43*T(3))
            ARRY1(2,2)=-HQX(0,1)*(S(5)-T(5))
            ARRY1(3,2)=-HQX(1,2)*(S(6)-T(6))
            ARRY1(4,2)= HQX(0,2)*(S(9)-T(9))
            H(0,0,1,1)= FQX(2,0)*T(2)+U(0)
            H(1,0,1,1)=-FQX(1,1)*(TEMP1-U(1))
            H(3,0,1,1)=-FQX(2,1)*T(6)+ARRY1(2,2)
            H(1,1,1,1)= HQX(0,2)*TEMP1 +TEMP2 +HQX(2,1)*T(2)-
     2                  FQX(2,2)*(T(3)-QPS*T(4)+T(3))
            H(3,1,1,1)= FQX(1,2)*(T(6)-QPS*T(7))+ARRY1(3,2)
            H(2,2,1,1)= TEMP0 +T2211
            H(3,3,1,1)= TEMP0 +T2211 +ARRY1(4,2)+FQX(2,2)*T(10)
            T2021 =-HQX(1,1)*T(2)
            T2121 = HQX(0,2)*TEMP1-H43*H43*T(2)
            T3221 = HQX(1,2)*T(6)
            H(2,0,2,1)= T2021
            H(2,1,2,1)= T2121
            H(3,2,2,1)= T3221
            H(2,0,1,2)= T2021
            H(2,1,1,2)= T2121
            H(3,2,1,2)= T3221
            H(0,0,2,2)= U(0)
            H(1,0,2,2)= FQX(1,1)*U(1)
            H(3,0,2,2)= ARRY1(2,2)
            H(1,1,2,2)= FQX(2,2)*U(2)+TEMP0
            H(3,1,2,2)= ARRY1(3,2)
            H(2,2,2,2)= HQX(0,2)*T(2)+TEMP0
            H(3,3,2,2)= ARRY1(4,2)+TEMP0
         ENDIF
         MM=JTYPE-4
C
         L=0
         DO 360 M=1,MM
            K=2-L
            ZOW(2,M)  = ZOW(1,M)-ZOW(2,M)*QPS
            ARRY2(2,M)= HQX(0,1)*ZOW(0,M)
            ARRY2(3,M)= HQX(1,2)*ZOW(1,M)
            ARRY2(4,M)=-HQX(0,2)*ZOW(4,M)
            H(2,0,K,L)= ARRY2(2,M)
            H(2,1,K,L)= ARRY2(3,M)
            H(3,2,K,L)= ARRY2(4,M)
            IF(M.EQ.1) THEN
               ARRY1(2,1)= ARRY2(2,1)
               ARRY1(3,1)= ARRY2(3,1)
               ARRY1(4,1)= ARRY2(4,1)
            ELSE
               ARRY3(2,2)= ARRY2(2,2)
               ARRY3(3,2)= ARRY2(3,2)
               ARRY3(4,2)= ARRY2(4,2)
            ENDIF
  360    L=K
         I=0
         L=0
         DO 380 M=1,MM
            ARRY3(1,M)= H43*(XOV(0,M)-X43*XOV(1,M))
            H(0,0,0,L)= XOV(0,M)
            H(1,0,0,L)= FQX(1,1)*XOV(1,M)
            H(3,0,0,L)=-FQX(0,1)*XOV(3,M)
            H(1,1,0,L)= FQX(2,2)*XOV(2,M)+ARRY3(1,M)
            H(3,1,0,L)=-FQX(1,2)*XOV(4,M)
            H(2,2,0,L)= ARRY3(1,M)
            H(3,3,0,L)= FQX(0,2)*XOV(5,M)+ARRY3(1,M)
            IF(M.EQ.2) THEN
               H(3,0,0,L)= H(3,0,0,L)+ARRY3(2,M)
               H(3,1,0,L)= H(3,1,0,L)+ARRY3(3,M)
            ENDIF
            ARRY1(1,M)= H43*(YOC(0,M)-X43*YOC(1,M))
            IF(M.EQ.2) THEN
               ARRY1(1,M)= ARRY1(1,M) +TEMP0
               DO 370 K=0,2
                  YOC(  K,M)= YOC(  K,M)+U(  K)
                  U(3+K)=PT5*(AOB(K,1)+AOB(K,2))
                  YOC(3+K,M)= YOC(3+K,M)-U(3+K)
  370          CONTINUE
            ENDIF
            H(0,0,3,L)= YOC(0,M)
            H(1,0,3,L)= FQX(1,1)*YOC(1,M)
            H(3,0,3,L)=-FQX(0,1)*YOC(3,M)+ARRY1(2,M)
            H(1,1,3,L)= FQX(2,2)*YOC(2,M)+ARRY1(1,M)
            H(3,1,3,L)=-FQX(1,2)*YOC(4,M)+ARRY1(3,M)
            H(2,2,3,L)= ARRY1(1,M)
            H(3,3,3,L)= FQX(0,2)*YOC(5,M)+ARRY1(1,M)
            IF(M.EQ.2)
     1      H(3,3,3,L)= H(3,3,3,L)+ARRY1(4,2)+
     2                  FQX(0,2)*(PT5*T(2)+U(5)-(AOB(4,1)+AOB(4,2)))
            ARRY2(1,M)=-HQX(1,1)*(ZOW(0,M)-X43*ZOW(1,M))
            J=1-I
            H(0,0,J,I)=-FQX(1,0)*ZOW(0,M)
            H(1,0,J,I)=-FQX(2,1)*ZOW(1,M)+ARRY2(2,M)
            H(3,0,J,I)= FQX(1,1)*ZOW(4,M)
            H(1,1,J,I)= FQX(1,2)*ZOW(2,M)+ARRY2(1,M)
            H(3,1,J,I)= FQX(2,2)*ZOW(5,M)+ARRY2(4,M)
            H(2,2,J,I)= ARRY2(1,M)
            H(3,3,J,I)=-FQX(1,2)*ZOW(7,M)+ARRY2(1,M)
         I=J
  380    L=3
         IF(JTYPE.EQ.6) THEN
            I=3
            L=3
            DO 390 M=1,MM
               ARRY3(1,M)= HQX(1,1)*(AOB(0,M)-X43*AOB(1,M))
               ARRY3(2,M)= PT5*T(2)-AOB(4,M)
               ARRY3(3,M)= PT5*T(3)-AOB(5,M)
               K=4-L
               H(0,0,K,L)=-FQX(1,0)*AOB(0,M)
               H(1,0,K,L)= FQX(0,1)*(AOB(0,M)-QPS*AOB(1,M))-H43*AOB(0,M)
               H(3,0,K,L)=-FQX(1,1)*ARRY3(2,M)
               H(1,1,K,L)= FQX(1,2)*(AOB(1,M)-QPS*AOB(2,M))-ARRY3(1,M)
               H(3,1,K,L)= HQX(0,2)*ARRY3(2,M)-FQX(2,2)*ARRY3(3,M)
               H(2,2,K,L)=-ARRY3(1,M)
               H(3,3,K,L)= FQX(1,2)*(T(6)-AOB(7,M))-ARRY3(1,M)
               J=5-I
               H(2,0,J,I)= H43*AOB(0,M)
               H(2,1,J,I)= HQX(1,2)*AOB(1,M)
               H(3,2,J,I)= HQX(0,2)*ARRY3(2,M)
            I=J
  390       L=K
         ENDIF
      ENDIF
  999 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------
C*MODULE RYSPOL  *DECK RYSSMT
      SUBROUTINE xRYSSMT(C,S,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(10,10),S(10,10),V(10),Y(10)
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00)
C
C     ROUTINE RETURNS AN N BY N TRIANGULAR MATRIX C SUCH THAT
C     C(TRANSPOSE)SC=I,  WHERE I IS AN N BY N IDENTITY MATRIX.
C
      DO 10 I=1,N
      DO 10 J=1,I
   10 C(I,J)= ZERO
C
      DO 100 J=1,N
         KMAX=J-1
         FAC=S(J,J)
         IF(KMAX.EQ.0) GO TO 60
C
         DO 20 K=1,KMAX
            V(K)= ZERO
            Y(K)=S(K,J)
   20    CONTINUE
         DO 50 K=1,KMAX
            DOT= ZERO
            DO 30 I=1,K
               DOT=C(I,K)*Y(I)+DOT
   30       CONTINUE
            DO 40 I=1,K
               V(I)=V(I)-DOT*C(I,K)
   40       CONTINUE
            FAC=FAC-DOT*DOT
   50    CONTINUE
C
   60    CONTINUE
C
C               THE ORIGINAL CODE LEADS TO DIVIDE BY
C               ZERO ON VAX COMPUTERS FOR X>3100 OR SO.
C--      FAC=ONE/SQRT(FAC)
C
         IF(FAC.GT.ZERO) THEN
            FAC=  ONE/SQRT(FAC)
         ELSE
            FAC = 1.0D+25
         END IF
C
         C(J,J)=FAC
         IF(KMAX.EQ.0) GO TO 100
         DO 70 K=1,KMAX
            C(K,J)=FAC*V(K)
   70    CONTINUE
  100 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------
C*MODULE RYSPOL  *DECK RYSFUN
      SUBROUTINE xRYSFUN(X,N,FF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FF(19)
      DIMENSION XMAX(9)
      PARAMETER (ZERO=0.0D+00, PT5=0.5D+00, ONE=1.0D+00, TWO=2.0D+00)
C          THE FIRST FIVE ROOTS ARE DONE IN SPECIALIZED ROUTINES.
      DATA XMAX/5*0.0D+00,
     *          9000000.0D+00, 750000.0D+00, 150000.0D+00, 50000.0D+00/
C
      TOL   = 1.0D-14
      XX    = X+X
      FACMIN= XX
C
C        NOTE ADDED BY MWS, SEPTEMBER 1990:
C        VALUES X.GT.3200 AND NROOTS.EQ.6 GENERATE NEGATIVE ARGUMENTS
C        TO THE SQRT FUNCTION IN RYSSMT WHEN ONE USES E=10**(-35).
C        THE PROBLEM GOES AWAY, AND THE ROOTS AND WEIGHTS STABILIZE
C        FOR X OUT TO 5000 WHEN 10**(-40) TO 10**(-60) IS USED.
C        SINCE VAX UNDERFLOWS AT AROUND 10**(-38) WE HAVE TO BE MORE
C        CAREFUL THAN JUST SETTING E A FEW POWERS OF TEN DOWNWARD.
C        MICHEL'S ORIGINAL CODE IS SHOWN FIRST.   THE VALUE OF 360
C        IS WHERE IBMS GET INTO UNDERFLOWS (ABOUT 10**(-78)).  THE
C        1990 CHANGES INCLUDE THE COMPUTATION OF E IN THE ASYMPTOTIC
C        CODE A FEW LINES DOWN.
C
C--   E     = 1.0D-35
C--   IF(FACMIN.LT.360.0D+00) E=EXP(-X)
C--   IF(FACMIN.GT.80.0D+00) GO TO 100
C
      E = ZERO
      IF(FACMIN.LT.360.0D+00) E=EXP(-X)
      IF(FACMIN.GT.80.0D+00) GO TO 100
C
      TERM= ONE
      SUM = ONE
      FAC = N
      FAC = FAC+PT5
C
   10 FAC = FAC+ONE
      TERM= TERM*X/FAC
      SUM = SUM+TERM
      IF(FAC.LE.FACMIN) GO TO 10
      T=TERM
      S=SUM
      IF(T.GT.S*TOL) GO TO 10
C
      FAC=N+N+1
      FF(N+1)=SUM*E/FAC
      M=N-1
      FAC=M+M+1
C
   20 IF(M.LT.0) RETURN
      FF(M+1)=(E+XX*FF(M+2))/FAC
      M=M-1
      FAC=FAC-TWO
      GO TO 20
C
C     NOTE BY MWS, NOVEMBER 1997:  THIS CODE HAD BEEN SETTING -E- TO
C     1.0D-35 WITH DOWNWARDS ADJUSTMENTS TO 1.0D-45 ON MOST MACHINES.
C     MICHEL USUALLY USES A VALUE OF 1.0D-78 HERE, WHERE THE LOWER
C     END OF THE IBM MAINFRAME'S DYNAMIC RANGE LIES.  MIKE'S HACK AT
C     TRYING TO ACCOMODATE BOTH THE VAX, AND THE NEED FOR A LOWER GUESS
C     AT E IS SHOWN COMMENTED OUT.  INSTEAD WE WILL NOW TEST FOR THE
C     MAXIMUM X FOR WHICH 1D-78 WILL RETURN 9 DIGIT ACCURACY ROOTS AND
C     WEIGHTS FOR.  ANYONE WHO FINDS THEMSELVES HERE MIGHT BE ABLE TO
C     INTEGRATE ACCURATELY BY CHOOSING A STILL SMALLER VALUE OF -E-,
C     BUT SHOULD BE CAREFUL TO TEST THE NUMERICS OF THE SITUATION
C     INSTEAD OF JUST BLINDLY LOWERING THE THRESHHOLD FURTHER.
C
C  TO CONVERGE 9 SIGNIFICANT DIGITS (IS THIS ENOUGH???) FOR NROOTS=6
C          X    ROOT CONV.  WEIGHT CONV
C        5,000      -43        -45       (NUMBER IS POWER OF TEN FOR E)
C       10,000      -47        -48
C       20,000      -49        -50
C       50,000      -54        -54
C      100,000      -57        -58
C      200,000      -60        -61
C    1,000,000      -67        -69
C    9,000,000      -76        -78
C   10,000,000      -79        -79
C     SIMILAR TESTS SHOW THAT HIGHER POLYNOMIALS GO TO SMALLER X VALS
C
C     NO. OF SIGNIFICANT DIGITS ACCURACY USING E=1.0D-45 FOR NROOTS=6
C          X       ROOTS    WEIGHTS
C        5,000        9         9
C       10,000        8         8
C       15,000        7         5
C       20,000        4         4
C       30,000        3         2
C       40,000        2         2
C       50,000        1         0
C    THAT IS, FOR X=40000, ROOT 1 IS 2.442D-6 BUT SHOULD BE 2.469D-6,
C    SO ONLY TWO DIGITS OF ACCURACY IS ATTAINED.  AT X=58000, THE
C    GUESS OF 1D-45 FOR -E- FAILS TO FIND THE ROOTS AT ALL.
C
C     --- USE ASYMPTOTIC EXPANSION FOR LARGE ARGUMENTS ---
C
  100 CONTINUE
C---      IF(E.EQ.ZERO) THEN
C---         ESAVE = 1.0D-35
C---         DO 5 I=1,10
C---            E = ESAVE*0.1D+00
C---            IF(E.EQ.ZERO) GO TO 6
C---            ESAVE = E
C---    5    CONTINUE
C---    6    CONTINUE
C---         E = ESAVE
C---      END IF
C
      IF(E.EQ.ZERO) THEN
         NROOTS = N/2
         IF(X.GT.XMAX(NROOTS)) THEN
            WRITE(6,9010) NROOTS,X,XMAX(NROOTS)
cwc            CALL ABRT
            call aces_exit(1)
         END IF
         E = 1.0D-78
      END IF
C
      A= SQRT(0.7853981633974483096156608D+00/X)
      TMAX= A*TOL/E
      TERM= ONE/XX
      SUM = TERM
      FAC = ONE
C
  110 FAC = FAC-TWO
      TERM= FAC*TERM/XX
      SUM = TERM+SUM
      T   = TERM
      IF(ABS(T).GT.TMAX) GO TO 110
C
      FF(1)= A-E*SUM
      FAC  = -ONE
      M    = 0
C
  120 IF(M.EQ.N) RETURN
      M      = M+1
      FAC    = FAC+TWO
      FF(M+1)= (FAC*FF(M)-E)/XX
      GO TO 120
C
 9010 FORMAT(1X,'RYSFUN: UNABLE TO FIND ROOTS OF RYS POLYNOMIAL.'/
     *       1X,'NROOTS=',I4,' X=',F20.5,' XMAX=',F20.5)
      END
c-----------------------------------------------------------------
C*MODULE RYSPOL  *DECK RYSNOD
      SUBROUTINE xRYSNOD(A,RT,K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL GOPARR,DSKWRK,MASWRK
      DIMENSION A(10),RT(10)
      COMMON /xPAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      PARAMETER (ZERO=0.0D+00)
C
C     RETURNS IN RT(I) THE ITH ROOT OF A POLYNOMIAL OF ORDER K
C     WHOSE MTH COEFFICIENT IS STORED IN A(M+1).  IT IS ASSUMED
C     THAT THE INITIAL VALUES IN RT BRACKET THE FINAL VALUES.
C
      TOL=1.0D-11
      K1=K+1
      R2= ZERO
      P2=A(1)
      DO 100 M=1,K
         R1=R2
         P1=P2
         R2=RT(M)
         P2=A(K1)
         DO 10 I=1,K
            P2=P2*R2+A(K1-I)
   10    CONTINUE
         PROD=P1*P2
         IF(PROD.GE.ZERO) THEN
            WRITE(6,15) M,K,ME
cwc            CALL ABRT
            call aces_exit(1)
         END IF
         R5=R1
         P5=P1
         R6=R2
         P6=P2
   30    CONTINUE
         R3=R5
         P3=P5
         R4=R6
         P4=P6
         R =(R3*P4-R4*P3)/(P4-P3)
         DR=R4-R3
         DELTA=DR
         IF(ABS(DELTA).LT.TOL) GO TO 90
         DR=0.0625D+00*DR
         R5=R-DR
         IF(R5.LT.R3) R5=R3
         R6=R+DR
         IF(R6.GT.R4) R6=R4
         P5=A(K1)
         P6=P5
         DO 40 I=1,K
            P5=P5*R5+A(K1-I)
            P6=P6*R6+A(K1-I)
   40    CONTINUE
   45    CONTINUE
         PROD=P5*P6
         IF(PROD.LT. ZERO) GO TO 30
         PROD=P3*P5
         IF(PROD.GT. ZERO) GO TO 60
         R5=0.25D+00*R3+0.75D+00*R5
         P5=A(K1)
         DO 50 I=1,K
            P5=P5*R5 + A(K1-I)
   50    CONTINUE
         GO TO 45
C
   60    CONTINUE
         R6=0.25D+00*R4 + 0.75D+00*R6
         P6=A(K1)
         DO 70 I=1,K
            P6=P6*R6 + A(K1-I)
   70    CONTINUE
         GO TO 45
C
   90    CONTINUE
         RT(M)=R
  100 CONTINUE
      RETURN
C
   15 FORMAT(//1X,'RYSNOD: ROOT NUMBER ',I3,
     *            ' WAS NOT FOUND FOR POLYNOMIAL OF ORDER ',I3,
     *            ' ON NODE',I4//)
      END
c----------------------------------------------------------------------
C*MODULE INT2A   *DECK FORMS
      SUBROUTINE xFORMS(GHONDO)
CDIR$ INLINEALWAYS XFORMS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION GHONDO(*)
C
      COMMON /xDENS  / DKL(784),DIJ(784)
      COMMON /xINTDEX/ IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784),
     *                KLGT(784),KLX(784),KLY(784),KLZ(784)
      COMMON /xSHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     *                NIJ,IJ,KL,IJKL
      COMMON /xROOT  / XX,U(13),W(13),NROOTS
      COMMON /xXYZ   / XIN(31213),YIN(31213),ZIN(31213)
C
C     ----- FORM INTEGRALS OVER FUNCTIONS -----
C     DIMENSIONING XIN(81,5), AND ROLLING UP THE COMPUTATION
C     OF GHONDO IN A LOOP OF LENGTH NROOTS ADDS 33 SECONDS TO
C     A 240 SECOND INTEGRAL COMPUTATION JOB.  LEAVE IT UNROLLED.
C
      GO TO (100,140,180,220,260,300,340,380,420),NROOTS
C
  100 DO 120 I = 1,IJ
      D1 = DIJ(I)
      NX = IJX(I)
      NY = IJY(I)
      NZ = IJZ(I)
      N1 = IJGT(I)
      MAX = IK(I)
      DO 120 K = 1,MAX
      MX = NX+KLX(K)
      MY = NY+KLY(K)
      MZ = NZ+KLZ(K)
      N = N1+KLGT(K)
  120 GHONDO(N) = ( XIN(MX )*YIN(MY )*ZIN(MZ ) )*D1*DKL(K)+GHONDO(N)
      RETURN
C
  140 DO 160 I = 1,IJ
      D1 = DIJ(I)
      NX = IJX(I)
      NY = IJY(I)
      NZ = IJZ(I)
      N1 = IJGT(I)
      MAX = IK(I)
      DO 160 K = 1,MAX
      MX = NX+KLX(K)
      MY = NY+KLY(K)
      MZ = NZ+KLZ(K)
      N = N1+KLGT(K)
  160 GHONDO(N) = ( XIN(MX )*YIN(MY )*ZIN(MZ ) +XIN(MX+625)*YIN(MY+625)*
     +     ZIN(MZ+625) )*D1*DKL(K)+GHONDO(N)
      RETURN
C
  180 DO 200 I = 1,IJ
      D1 = DIJ(I)
      NX = IJX(I)
      NY = IJY(I)
      NZ = IJZ(I)
      N1 = IJGT(I)
      MAX = IK(I)
      DO 200 K = 1,MAX
      MX = NX+KLX(K)
      MY = NY+KLY(K)
      MZ = NZ+KLZ(K)
      N = N1+KLGT(K)
  200 GHONDO(N) = ( XIN(MX )*YIN(MY )*ZIN(MZ ) +XIN(MX+625)*YIN(MY+625)*
     +     ZIN(MZ+625) +XIN(MX+1250)*YIN(MY+1250)*ZIN(MZ+1250) )*D1*
     +     DKL(K)+GHONDO(N)
      RETURN
C
  220 DO 240 I = 1,IJ
      D1 = DIJ(I)
      NX = IJX(I)
      NY = IJY(I)
      NZ = IJZ(I)
      N1 = IJGT(I)
      MAX = IK(I)
      DO 240 K = 1,MAX
      MX = NX+KLX(K)
      MY = NY+KLY(K)
      MZ = NZ+KLZ(K)
      N = N1+KLGT(K)
  240 GHONDO(N) = ( XIN(MX )*YIN(MY )*ZIN(MZ ) +XIN(MX+625)*YIN(MY+625)*
     +     ZIN(MZ+625) +XIN(MX+1250)*YIN(MY+1250)*ZIN(MZ+1250) +XIN(MX+
     +     1875)*YIN(MY+1875)*ZIN(MZ+1875) )*D1*DKL(K)+GHONDO(N)
      RETURN
C
  260 DO 280 I = 1,IJ
      D1 = DIJ(I)
      NX = IJX(I)
      NY = IJY(I)
      NZ = IJZ(I)
      N1 = IJGT(I)
      MAX = IK(I)
      DO 280 K = 1,MAX
      MX = NX+KLX(K)
      MY = NY+KLY(K)
      MZ = NZ+KLZ(K)
      N = N1+KLGT(K)
  280 GHONDO(N) = ( XIN(MX )*YIN(MY )*ZIN(MZ ) +XIN(MX+625)*YIN(MY+625)*
     +     ZIN(MZ+625) +XIN(MX+1250)*YIN(MY+1250)*ZIN(MZ+1250) +XIN(MX+
     +     1875)*YIN(MY+1875)*ZIN(MZ+1875) +XIN(MX+2500)*YIN(MY+2500)*
     +     ZIN(MZ+2500) )*D1*DKL(K)+GHONDO(N)
      RETURN
C
  300 DO 320 I = 1,IJ
      D1 = DIJ(I)
      NX = IJX(I)
      NY = IJY(I)
      NZ = IJZ(I)
      N1 = IJGT(I)
      MAX = IK(I)
      DO 320 K = 1,MAX
      MX = NX+KLX(K)
      MY = NY+KLY(K)
      MZ = NZ+KLZ(K)
      N = N1+KLGT(K)
  320 GHONDO(N) = ( XIN(MX )*YIN(MY )*ZIN(MZ ) +XIN(MX+625)*YIN(MY+625)*
     +     ZIN(MZ+625) +XIN(MX+1250)*YIN(MY+1250)*ZIN(MZ+1250) +XIN(MX+
     +     1875)*YIN(MY+1875)*ZIN(MZ+1875) +XIN(MX+2500)*YIN(MY+2500)*
     +     ZIN(MZ+2500) +XIN(MX+3125)*YIN(MY+3125)*ZIN(MZ+3125))
     +     *D1*DKL(K)+GHONDO(N)
      RETURN
C
  340 DO 360 I = 1,IJ
      D1 = DIJ(I)
      NX = IJX(I)
      NY = IJY(I)
      NZ = IJZ(I)
      N1 = IJGT(I)
      MAX = IK(I)
      DO 360 K = 1,MAX
      MX = NX+KLX(K)
      MY = NY+KLY(K)
      MZ = NZ+KLZ(K)
      N = N1+KLGT(K)
  360 GHONDO(N) = ( XIN(MX )*YIN(MY )*ZIN(MZ ) +XIN(MX+625)*YIN(MY+625)*
     +     ZIN(MZ+625) +XIN(MX+1250)*YIN(MY+1250)*ZIN(MZ+1250) +XIN(MX+
     +     1875)*YIN(MY+1875)*ZIN(MZ+1875) +XIN(MX+2500)*YIN(MY+2500)*
     +     ZIN(MZ+2500) +XIN(MX+3125)*YIN(MY+3125)*ZIN(MZ+3125) +XIN(MX+
     +     3750)*YIN(MY+3750)*ZIN(MZ+3750))
     +     *D1*DKL(K)+GHONDO(N)
      RETURN
C
  380 DO 400 I = 1,IJ
      D1 = DIJ(I)
      NX = IJX(I)
      NY = IJY(I)
      NZ = IJZ(I)
      N1 = IJGT(I)
      MAX = IK(I)
      DO 400 K = 1,MAX
      MX = NX+KLX(K)
      MY = NY+KLY(K)
      MZ = NZ+KLZ(K)
      N = N1+KLGT(K)
  400 GHONDO(N) = ( XIN(MX )*YIN(MY )*ZIN(MZ ) +XIN(MX+625)*YIN(MY+625)*
     +     ZIN(MZ+625) +XIN(MX+1250)*YIN(MY+1250)*ZIN(MZ+1250) +XIN(MX+
     +     1875)*YIN(MY+1875)*ZIN(MZ+1875) +XIN(MX+2500)*YIN(MY+2500)*
     +     ZIN(MZ+2500) +XIN(MX+3125)*YIN(MY+3125)*ZIN(MZ+3125) +XIN(MX+
     +     3750)*YIN(MY+3750)*ZIN(MZ+3750) +XIN(MX+4375)*YIN(MY+4375)*
     +     ZIN(MZ+4375))*D1*DKL(K)+GHONDO(N)
      RETURN
C
  420 DO 440 I = 1,IJ
      D1 = DIJ(I)
      NX = IJX(I)
      NY = IJY(I)
      NZ = IJZ(I)
      N1 = IJGT(I)
      MAX = IK(I)
      DO 440 K = 1,MAX
      MX = NX+KLX(K)
      MY = NY+KLY(K)
      MZ = NZ+KLZ(K)
      N = N1+KLGT(K)
  440 GHONDO(N) = ( XIN(MX )*YIN(MY )*ZIN(MZ ) +XIN(MX+625)*YIN(MY+625)*
     +     ZIN(MZ+625) +XIN(MX+1250)*YIN(MY+1250)*ZIN(MZ+1250) +XIN(MX+
     +     1875)*YIN(MY+1875)*ZIN(MZ+1875) +XIN(MX+2500)*YIN(MY+2500)*
     +     ZIN(MZ+2500) +XIN(MX+3125)*YIN(MY+3125)*ZIN(MZ+3125) +XIN(MX+
     +     3750)*YIN(MY+3750)*ZIN(MZ+3750) +XIN(MX+4375)*YIN(MY+4375)*
     +     ZIN(MZ+4375) +XIN(MX+5000)*YIN(MY+5000)*ZIN(MZ+5000))
     +     *D1*DKL(K)+GHONDO(N)
      RETURN
      END
c-----------------------------------------------------------------------
C*MODULE INT2A   *DECK XYZINT
      SUBROUTINE xXYZINT
CDIR$ INLINEALWAYS XXYZINT
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL N0,N1,M0,M1
C
      COMMON /xSETINT/ I(13),K(13),NIMAX,NJMAX,NKMAX,NLMAX,NMAX,MMAX
     +               ,BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00
     +               ,DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL
      COMMON /xXYZ   / XINT(31213),YINT(31213),ZINT(31213)
C
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00)
C
      N0 = NMAX .EQ. 0
      N1 = NMAX .LE. 1
      M0 = MMAX .EQ. 0
      M1 = MMAX .LE. 1
C
C     ----- I(0,0) -----
C
      I1 = I(1)
      XINT(I1) = ONE
      YINT(I1) = ONE
      ZINT(I1) = F00
      IF (N0 .AND. M0) RETURN
      I2 = I(2)
      K2 = K(2)
      CP10 = B00
      IF (N0) GO TO 100
C
C     ----- I(1,0) -----
C
      XINT(I2) = XC00
      YINT(I2) = YC00
      ZINT(I2) = ZC00*F00
      IF (M0) GO TO 120
C
C     ----- I(0,1) -----
C
  100 I3 = I1+K2
      XINT(I3) = XCP00
      YINT(I3) = YCP00
      ZINT(I3) = ZCP00*F00
      IF (N0) GO TO 120
C
C     ----- I(1,1) -----
C
      I3 = I2+K2
      XINT(I3) = XCP00*XINT(I2)+CP10
      YINT(I3) = YCP00*YINT(I2)+CP10
      ZINT(I3) = ZCP00*ZINT(I2)+CP10*F00
  120 IF (N1) GO TO 180
      C10 = ZERO
      I3 = I1
      I4 = I2
      DO 160 N = 2,NMAX
      C10 = C10+B10
C
C     ----- I(N,0) -----
C
      I5 = I(N+1)
      XINT(I5) = C10*XINT(I3)+XC00*XINT(I4)
      YINT(I5) = C10*YINT(I3)+YC00*YINT(I4)
      ZINT(I5) = C10*ZINT(I3)+ZC00*ZINT(I4)
      IF (M0) GO TO 140
      CP10 = CP10+B00
C
C     ----- I(N,1) -----
C
      I3 = I5+K2
      XINT(I3) = XCP00*XINT(I5)+CP10*XINT(I4)
      YINT(I3) = YCP00*YINT(I5)+CP10*YINT(I4)
      ZINT(I3) = ZCP00*ZINT(I5)+CP10*ZINT(I4)
  140 I3 = I4
  160 I4 = I5
  180 IF (M1) GO TO 240
      CP01 = ZERO
      C01 = B00
      I3 = I1
      I4 = I1+K2
      DO 220 M = 2,MMAX
      CP01 = CP01+BP01
C
C     ----- I(0,M) -----
C
      I5 = I1+K(M+1)
      XINT(I5) = CP01*XINT(I3)+XCP00*XINT(I4)
      YINT(I5) = CP01*YINT(I3)+YCP00*YINT(I4)
      ZINT(I5) = CP01*ZINT(I3)+ZCP00*ZINT(I4)
      IF (N0) GO TO 200
      C01 = C01+B00
C
C     ----- I(1,M) -----
C
      I3 = I2+K(M+1)
      XINT(I3) = XC00*XINT(I5)+C01*XINT(I4)
      YINT(I3) = YC00*YINT(I5)+C01*YINT(I4)
      ZINT(I3) = ZC00*ZINT(I5)+C01*ZINT(I4)
  200 I3 = I4
  220 I4 = I5
  240 IF (N1 .OR. M1) GO TO 300
C
C     ----- I(N,M) -----
C
      C01 = B00
      K3 = K2
      DO 280 M = 2,MMAX
      K4 = K(M+1)
      C01 = C01+B00
      I3 = I1
      I4 = I2
      C10 = B10
      DO 260 N = 2,NMAX
      I5 = I(N+1)
      XINT(I5+K4) = C10*XINT(I3+K4)+XC00*XINT(I4+K4)+C01*XINT(I4+K3)
      YINT(I5+K4) = C10*YINT(I3+K4)+YC00*YINT(I4+K4)+C01*YINT(I4+K3)
      ZINT(I5+K4) = C10*ZINT(I3+K4)+ZC00*ZINT(I4+K4)+C01*ZINT(I4+K3)
      C10 = C10+B10
      I3 = I4
  260 I4 = I5
  280 K3 = K4
  300 IF (NJMAX .EQ. 0) GO TO 440
C
C     ----- I(NI,NJ,M) -----
C
      M = 0
      I5 = I(NMAX+1)
  320 MIN = NIMAX
      KM = K(M+1)
  340 N = NMAX
      I3 = I5+KM
  360 I4 = I(N)+KM
      XINT(I3) = XINT(I3)+DXIJ*XINT(I4)
      YINT(I3) = YINT(I3)+DYIJ*YINT(I4)
      ZINT(I3) = ZINT(I3)+DZIJ*ZINT(I4)
      I3 = I4
      N = N-1
      IF (N .GT. MIN) GO TO 360
      MIN = MIN+1
      IF (MIN .LT. NMAX) GO TO 340
      IF (NIMAX .EQ. 0) GO TO 420
      I3 = 25+KM+I1
      DO 400 NJ = 1,NJMAX
      I4 = I3
      DO 380 NI = 1,NIMAX
      XINT(I4) = XINT(I4+100)+DXIJ*XINT(I4-25)
      YINT(I4) = YINT(I4+100)+DYIJ*YINT(I4-25)
      ZINT(I4) = ZINT(I4+100)+DZIJ*ZINT(I4-25)
  380 I4 = I4+125
  400 I3 = I3+25
  420 M = M+1
      IF (M .LE. MMAX) GO TO 320
  440 IF (NLMAX .EQ. 0) GO TO 600
C
C     ----- I(NI,NJ,NK,NL) -----
C
      I5 = K(MMAX+1)
      IA = I1
      NI = 0
  460 NJ = 0
      IB = IA
  480 MIN = NKMAX
  500 M = MMAX
      I3 = IB+I5
  520 I4 = IB+K(M)
      XINT(I3) = XINT(I3)+DXKL*XINT(I4)
      YINT(I3) = YINT(I3)+DYKL*YINT(I4)
      ZINT(I3) = ZINT(I3)+DZKL*ZINT(I4)
      I3 = I4
      M = M-1
      IF (M .GT. MIN) GO TO 520
      MIN = MIN+1
      IF (MIN .LT. MMAX) GO TO 500
      IF (NKMAX .EQ. 0) GO TO 580
      I3 = IB+1
      DO 560 NL = 1,NLMAX
      I4 = I3
      DO 540 NK = 1,NKMAX
      XINT(I4) = XINT(I4+4)+DXKL*XINT(I4-1)
      YINT(I4) = YINT(I4+4)+DYKL*YINT(I4-1)
      ZINT(I4) = ZINT(I4+4)+DZKL*ZINT(I4-1)
  540 I4 = I4+5
  560 I3 = I3+1
  580 NJ = NJ+1
      IB = IB+25
      IF (NJ .LE. NJMAX) GO TO 480
      NI = NI+1
      IA = IA+125
      IF (NI .LE. NIMAX) GO TO 460
  600 CONTINUE
      RETURN
      END
c----------------------------------------------------------------------
C*MODULE MTHLIB  *DECK VCLR
      SUBROUTINE xVCLR(A,INCA,N)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION A(*)
C
      PARAMETER (ZERO=0.0D+00)
C
C     ----- ZERO OUT VECTOR -A-, USING INCREMENT -INCA- -----
C
      IF (INCA .NE. 1) GO TO 200
      DO 110 L=1,N
         A(L) = ZERO
  110 CONTINUE
      RETURN
C
  200 CONTINUE
      LA=1-INCA
      DO 210 L=1,N
         LA=LA+INCA
         A(LA) = ZERO
  210 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------
C*MODULE INT2A   *DECK S0000
      SUBROUTINE xS0000(GHONDO,DDIJ)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL IANDJ,KANDL,SAME,OUT
C
#include "gamess.par"
C
      DIMENSION GHONDO(*),DDIJ(16*MXG2)
C
      COMMON /xIJGNRL/ A(MXG2),R(MXG2),X1(MXG2),Y1(MXG2),Z1(MXG2),
     *                IJD(784)
      COMMON /xGOUT  / GPOPLE(768),NORG
      COMMON /xMISC  / IANDJ,KANDL,SAME
      COMMON /xSHLINF/ AG(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH),
     *                CFA(MXGSH),CGA(MXGSH),
     *                BG(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH),
     *                CFB(MXGSH),CGB(MXGSH),
     *                CG(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH),
     *                CFC(MXGSH),CGC(MXGSH),
     *                DG(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH),
     *                CFD(MXGSH),CGD(MXGSH),
     *                XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK,
     *                NGA,NGB,NGC,NGD
      COMMON /xSHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     *                NIJ,IJ,KL,IJKL
      COMMON /xSHLT  / TOL,CUTOFF,OUT,ICOUNT
C
      PARAMETER (PI252=34.986836655250D+00, PIE4=7.85398163397448D-01,
     *           ZERO=0.0D+00, ONE=1.0D+00)
C
C     SPECIAL SSSS INTEGRAL ROUTINE WHEN USING HONDO INTEGRALS
C
      GGOUT = ZERO
      LGMAX = NGD
      DO 300 KG = 1,NGC
      BK = CG(KG)
      BRRK = BK*RRK
      BXK = BK*XK
      BYK = BK*YK
      BZK = BK*ZK
      CSK = CSC(KG)
      IF (KANDL) LGMAX = KG
      DO 280 LG = 1,LGMAX
      BL = DG(LG)
      BB = BK+BL
      BBINV = ONE/BB
      DUM = BL*BRRK*BBINV
      IF (DUM .GT. TOL) GO TO 280
      BBRRK = DUM
      D2 = CSD(LG)*CSK*BBINV
      IF (KANDL .AND. LG .NE. KG) D2 = D2+D2
      BBX = (BXK+BL*XL)*BBINV
      BBY = (BYK+BL*YL)*BBINV
      BBZ = (BZK+BL*ZL)*BBINV
      SUM = ZERO
      NN = 1
      DO 260 N = 1,NIJ
      DUM = BBRRK+R(N)
      IF (DUM .GT. TOL) GO TO 260
      EXPE = EXP(-DUM)
      AA = A(N)
      AB = AA+BB
      DUM = X1(N)-BBX
      XX = DUM*DUM
      DUM = Y1(N)-BBY
      XX = DUM*DUM+XX
      DUM = Z1(N)-BBZ
      XX = DUM*DUM+XX
      X = XX*AA*BB/AB
C
      IF (X .GT. 5.0D+00) GO TO 160
      IF (X .GT. 1.0D+00) GO TO 120
      IF (X .GT. 3.0D-07) GO TO 100
      WW1 = 1.0D+00-X/3.0D+00
      GO TO 240
C
  100 CONTINUE
      F1 = ((((((((-8.36313918003957D-08*X+1.21222603512827D-06 )*X-
     +     1.15662609053481D-05 )*X+9.25197374512647D-05 )*X-
     +     6.40994113129432D-04 )*X+3.78787044215009D-03 )*X-
     +     1.85185172458485D-02 )*X+7.14285713298222D-02 )*X-
     +     1.99999999997023D-01 )*X+3.33333333333318D-01
      WW1 = (X+X)*F1+EXP(-X)
      GO TO 240
C
  120 CONTINUE
      IF (X .GT. 3.0D+00) GO TO 140
      Y = X-2.0D+00
      F1 = ((((((((((-1.61702782425558D-10*Y+1.96215250865776D-09 )*Y-
     +     2.14234468198419D-08 )*Y+2.17216556336318D-07 )*Y-
     +     1.98850171329371D-06 )*Y+1.62429321438911D-05 )*Y-
     +     1.16740298039895D-04 )*Y+7.24888732052332D-04 )*Y-
     +     3.79490003707156D-03 )*Y+1.61723488664661D-02 )*Y-
     +     5.29428148329736D-02 )*Y+1.15702180856167D-01
      WW1 = (X+X)*F1+EXP(-X)
      GO TO 240
C
  140 CONTINUE
      Y = X-4.0D+00
      F1 = ((((((((((-2.62453564772299D-11*Y+3.24031041623823D-10 )*Y-
     +     3.614965656163D-09)*Y+3.760256799971D-08)*Y-
     +     3.553558319675D-07)*Y+3.022556449731D-06)*Y-
     +     2.290098979647D-05)*Y+1.526537461148D-04)*Y-
     +     8.81947375894379D-04 )*Y+4.33207949514611D-03 )*Y-
     +     1.75257821619926D-02 )*Y+5.28406320615584D-02
      WW1 = (X+X)*F1+EXP(-X)
      GO TO 240
C
  160 CONTINUE
      IF (X .GT. 15.0D+00) GO TO 200
      E = EXP(-X)
      IF (X .GT. 10.0D+00) GO TO 180
      XINV = ONE/X
      WW1 = (((((( 4.6897511375022D-01*XINV-6.9955602298985D-01)*XINV +
     +     5.3689283271887D-01)*XINV-3.2883030418398D-01)*XINV +
     +     2.4645596956002D-01)*XINV-4.9984072848436D-01)*XINV -
     +     3.1501078774085D-06)*E + SQRT(PIE4*XINV)
      GO TO 240
C
  180 CONTINUE
      XINV = ONE/X
      WW1 = (((-1.8784686463512D-01*XINV+2.2991849164985D-01)*XINV
     +         -4.9893752514047D-01)*XINV-2.1916512131607D-05)*E
     +         + SQRT(PIE4*XINV)
      GO TO 240
C
  200 CONTINUE
      IF (X .GT. 33.0D+00) GO TO 220
      XINV = ONE/X
      E = EXP(-X)
      WW1 = (( 1.9623264149430D-01*XINV-4.9695241464490D-01)*XINV -
     +     6.0156581186481D-05)*E + SQRT(PIE4*XINV)
      GO TO 240
C
  220 WW1 = SQRT(PIE4/X)
C
  240 SUM = SUM+DDIJ(NN)*WW1*EXPE/SQRT(AB)
  260 NN = NN+16
      GGOUT = GGOUT+D2*SUM
  280 CONTINUE
  300 CONTINUE
      GHONDO(NORG+1) = GGOUT*PI252*QQ4
      RETURN
      END
c-----------------------------------------------------------------
C*MODULE SCFLIB  *DECK SHLDEN
cwc1
cwc      SUBROUTINE SHLDEN(SCFTYP,DA,DB,DGVB,DSH,IA,L1,L2,NSH2,NFO)
      SUBROUTINE xSHLDEN(iuhf,DA,DB,DSH,imap,nbas)
c imap contains the mapping from 2 (nbas*nbas) to 1 dimension
c (specific for ACES)
cwc0
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL ALPHA,BETA
C
cwc1
cwc      DIMENSION DA(L2),DB(L2),DGVB(L2,NFO),DSH(NSH2),IA(L1)
      DIMENSION DA(*),DB(*),DSH(*)
#include "gamess.par"
      COMMON /xIJPAIR/ IA(MXAO)
c iorder contains mapping between the orbital numbers in GAMESS and those
c in the calling program (such as ACES)
      common /GMSorder/ coeint(MXAO),iorder(MXAO),
     >               nAfrst(MXSH),nAstep(MXSH),sphrcl
      dimension imap(nbas,nbas)
      logical sphrcl
cwc0
C
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
C
      PARAMETER (ZERO=0.0D+00)
C
      DATA RHF,UHF,ROHF,GVB/8HRHF     ,8HUHF     ,8HROHF    ,8HGVB     /
C
C     ----- COMPRESS DENSITY OVER AOS TO DENSITY OVER SHELLS -----
C
cwc1
c      write (99,*) 'shlden begin'
      if (iuhf.eq.0) then
         SCFTYP=RHF
      else
         SCFTYP=UHF
      end if
cwc0
      IJSH=0
      ALPHA = SCFTYP.EQ.RHF  .OR.  SCFTYP.EQ.UHF  .OR.  SCFTYP.EQ.ROHF
      BETA  =                      SCFTYP.EQ.UHF  .OR.  SCFTYP.EQ.ROHF
cwc1
      if (sphrcl) then
         call SHDsph (alpha,beta,DA,DB,DSH,imap,nbas)
         return
      end if   
cwc0
      DO 240 ISH=1,NSHELL
         MINI = KLOC(ISH)
         MAXI = MINI + KMAX(ISH) - KMIN(ISH)
         DO 230 JSH=1,ISH
            MINJ = KLOC(JSH)
            MAXJ = MINJ + KMAX(JSH) - KMIN(JSH)
            IJSH = IJSH+1
            DMAX = ZERO
            DO 130 I=MINI,MAXI
cwc1
c I is the orbital number as known to GAMESS, but in the density matrix
c this orbital is located according to ACES's order
               ibas=iorder(i)
cwc0
               IF(ISH.EQ.JSH) MAXJ=I
               DO 120 J=MINJ,MAXJ
                  IJ = IA(I) + J
cwc1
                  jbas=iorder(j)
                  ij=imap(ibas,jbas)
cwc0
                  IF(ALPHA) THEN
                     IF(ABS(DA(IJ)).GT.DMAX) DMAX = ABS(DA(IJ))
                  END IF
                  IF(BETA) THEN
                     IF(ABS(DB(IJ)).GT.DMAX) DMAX = ABS(DB(IJ))
                  END IF
cwc                  IF(SCFTYP.EQ.GVB) THEN
cwc                     DO 110 IHAM=1,NFO
cwc                        IF(ABS(DGVB(IJ,IHAM)).GT.DMAX)
cwc     *                               DMAX=ABS(DGVB(IJ,IHAM))
cwc 110                                    CONTINUE
cwc                  END IF
 120           CONTINUE
 130        CONTINUE
            DSH(IJSH) = DMAX
c            write (99,'(2i3,g20.12)') ish,jsh,dmax
 230     CONTINUE
 240  CONTINUE
c      write (99,*) 'shlden end'
      RETURN
      END
c-----------------------------------------------------------------
C*MODULE INT2A   *DECK SCHWDN
cwc1
cwc      DOUBLE PRECISION FUNCTION SCHWDN(DSH,ISH,JSH,KSH,LSH,IA)
      DOUBLE PRECISION FUNCTION xSCHWDN(DSH,ISH,JSH,KSH,LSH)
cwc0
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
cwc1
cwc      DIMENSION DSH(*),IA(*)
      DIMENSION DSH(*)
#include "gamess.par"
      COMMON /xIJPAIR/ IA(MXAO)
cwc0
C
      PARAMETER (FOUR=4.0D+00)
C
C     ----- FIND MAXIMUM DENSITY CONTRIBUTION TO THIS SHELL SET -----
C     -DSH- IS THE DENSITY MATRIX ALREADY COMPRESSED TO SHELLS
C
      IJ = IA(ISH)+JSH
      IK = IA(ISH)+KSH
      IL = IA(ISH)+LSH
      KL = IA(KSH)+LSH
      JK = IA(JSH)+KSH
      JL = IA(JSH)+LSH
      IF(JSH.LT.KSH) JK=IA(KSH)+JSH
      IF(JSH.LT.LSH) JL=IA(LSH)+JSH
      xSCHWDN=MAX(FOUR*DSH(IJ),FOUR*DSH(KL),
     *           DSH(JL),DSH(JK),DSH(IL),DSH(IK))
      RETURN
      END
c-----------------------------------------------------------------
C*MODULE INT2A   *DECK ZQOUT
      SUBROUTINE xZQOUT(GHONDO)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION GHONDO(*)
C
      LOGICAL IANDJ,KANDL,SAME
C
      COMMON /xINTDEX/ IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784),
     *                KLGT(784),KLX(784),KLY(784),KLZ(784)
      COMMON /xMISC  / IANDJ,KANDL,SAME
      COMMON /xSHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     *                NIJX,IJ,KL,IJKL
C
      PARAMETER (ZERO=0.0D+00)
C
C     ----- ZERO HONDO CONVENTIONAL INTEGRAL OUTPUT REGION -----
C
      IJN = 0
      JMAX = MAXJ
      DO 260 I = MINI,MAXI
         IF (IANDJ) JMAX = I
         DO 240 J = MINJ,JMAX
            IJN = IJN+1
            N1 = IJGT(IJN)
            LMAX = MAXL
            KLN = 0
            DO 220 K =  MINK,MAXK
               IF (KANDL) LMAX = K
               DO 200 L = MINL,LMAX
                  KLN = KLN+1
                  IF (SAME .AND. KLN .GT. IJN) GO TO 240
                  NN = N1+KLGT(KLN)
                  GHONDO(NN) = ZERO
  200          CONTINUE
  220       CONTINUE
  240    CONTINUE
  260 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------
C*MODULE INT2A   *DECK ZPKOUT
      SUBROUTINE xZPKOUT(II,JJ,KK,LL,GHONDO,SKIPA,SKIPB,SKIPC,NPSYM)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION GHONDO(*)
      DIMENSION IB(15),JB(15),KB(15),LB(15)
C
      LOGICAL SKIPA,SKIPB,SKIPC,NPSYM,IANDJ,KANDL,FIRST
C
#include "gamess.par"
C
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),MIN(MXSH),MAX(MXSH),NSHELL
C
      PARAMETER (ZERO=0.0D+00)
C
      SAVE FIRST,NORG1H,NORG2H,NORG3H,IB,JB,KB,LB
C
      DATA FIRST/.TRUE./
C
C     ----- ZERO HONDO PK INTEGRAL FORMATION REGION -----
C
      IF(FIRST) THEN
         FIRST=.FALSE.
         CALL xBASCHK(LMAX)
                       NANGM =  4
         IF(LMAX.EQ.2) NANGM =  6
         IF(LMAX.EQ.3) NANGM = 10
         IF(LMAX.EQ.4) NANGM = 15
         NORG1H = 1
         NORG2H = NORG1H + NANGM**4
         NORG3H = NORG2H + NANGM**4
         DO 100 I=1,NANGM
            LB(I) = I-1
            KB(I) = LB(I) * NANGM
            JB(I) = KB(I) * NANGM
            IB(I) = JB(I) * NANGM
 100     CONTINUE
      END IF
C
      NORG1 = NORG1H
      NORG2 = NORG2H
      NORG3 = NORG3H
C
      LIT = KTYPE(II)
      MINI = MIN(II)
      MAXI = MAX(II)
      MINJ = MIN(JJ)
      MAXJ = MAX(JJ)
      LKT = KTYPE(KK)
      MINK = MIN(KK)
      MAXK = MAX(KK)
      MINL = MIN(LL)
      MAXL = MAX(LL)
      IANDJ = II .EQ. JJ
      KANDL = KK .EQ. LL
C
C     TYPE = 1 FOR (II II II II)
C            2     (II JJ JJ JJ)
C            3     (II II KK KK) AND  LIT.GE.LKT
C            4     (II II KK KK) AND  LIT.LT.LKT
C            5     (II II II LL)
C            6     (II JJ KK KK)
C            7     (II JJ JJ LL)
C            8     (II II KK LL)
C            9     (II JJ KK LL)
C
      NTYP = 0
      IF(II.EQ.JJ .AND. JJ.EQ.KK .AND. KK.EQ.LL) NTYP = 1
      IF(II.GT.JJ .AND. JJ.EQ.KK .AND. KK.EQ.LL) NTYP = 2
      IF(II.EQ.JJ .AND. JJ.GT.KK .AND. KK.EQ.LL
     *                         .AND. LIT.GE.LKT) NTYP = 3
      IF(II.EQ.JJ .AND. JJ.GT.KK .AND. KK.EQ.LL
     *                         .AND. LIT.LT.LKT) NTYP = 4
      IF(II.EQ.JJ .AND. JJ.EQ.KK .AND. KK.GT.LL) NTYP = 5
      IF(II.GT.JJ .AND. JJ.GT.KK .AND. KK.EQ.LL) NTYP = 6
      IF(II.GT.JJ .AND. JJ.EQ.KK .AND. KK.GT.LL) NTYP = 7
      IF(II.EQ.JJ .AND. JJ.GT.KK .AND. KK.GT.LL) NTYP = 8
      IF(II.GT.JJ .AND. JJ.GT.KK .AND. KK.GT.LL) NTYP = 9
C
      IF(SKIPA .AND. .NOT.NPSYM) NORG2 = 1
      IF(SKIPC .AND. .NOT.NPSYM) NORG3 = NORG2H
      IF(SKIPB .AND. .NOT.NPSYM) NORG3 = 1
C
C     ----- BEGIN LOOPS OVER PRIMITIVES IN THIS SHELL -----
C
C     INTEGRAL TYPES N1,G1 FOR (I,J//K,L)
C                    N2,G2 FOR (I,K//J,L)
C                    N3,G3 FOR (I,L//J,K)
C
C                      HONDO INTEGRALS
C        N1 = IB(IA)+JB(JA)+KB(KA)+LB(LA)+NORG1
C        N2 = IB(IA)+JB(KA)+KB(JA)+LB(LA)+NORG2
C
      JMAX = MAXJ
      KMAX = MAXK
      LMAX = MAXL
      DO 860 I = MINI,MAXI
      IA = I-MINI+1
      N1I = NORG1 + IB(IA)
      N2I = NORG2 + IB(IA)
C
      IF (IANDJ) JMAX = I
      DO 840 J = MINJ,JMAX
      IF (JJ .EQ. KK) KMAX = J
      JA = J-MINJ+1
      N1IJ = N1I + JB(JA)
      N2IJ = N2I + KB(JA)
C
      DO 820 K = MINK,KMAX
      KA = K-MINK+1
      N1IJK = N1IJ + KB(KA)
      N2IJK = N2IJ + JB(KA)
C
      IF (KANDL) LMAX = K
      DO 800 L = MINL,LMAX
      LA = L-MINL+1
      N1 = N1IJK + LB(LA)
      N2 = N2IJK + LB(LA)
C
      GO TO (200,220,230,250,270,280,290,300,310),NTYP
  200 IF (IA .EQ. JA) GO TO 210
      N3 = IB(IA)+JB(LA)+KB(JA)+LB(KA)+NORG3
      GO TO 400
  210 N3 = IB(JA)+JB(KA)+KB(IA)+LB(LA)+NORG3
      GO TO 400
  220 N3 = IB(IA)+JB(LA)+KB(JA)+LB(KA)+NORG3
      GO TO 400
  230 IF (IA .EQ. JA) GO TO 240
      N3 = IB(IA)+JB(LA)+KB(JA)+LB(KA)+NORG3
      GO TO 400
  240 N3 = IB(JA)+JB(KA)+KB(IA)+LB(LA)+NORG3
      GO TO 400
  250 IF (KA .EQ. LA) GO TO 260
      N3 = IB(JA)+JB(KA)+KB(IA)+LB(LA)+NORG3
      GO TO 400
  260 N3 = IB(IA)+JB(LA)+KB(JA)+LB(KA)+NORG3
      GO TO 400
  270 N3 = IB(JA)+JB(KA)+KB(IA)+LB(LA)+NORG3
      GO TO 400
  280 N3 = IB(IA)+JB(LA)+KB(JA)+LB(KA)+NORG3
      GO TO 400
  290 N3 = IB(IA)+JB(LA)+KB(JA)+LB(KA)+NORG3
      GO TO 400
  300 N3 = IB(JA)+JB(KA)+KB(IA)+LB(LA)+NORG3
      GO TO 400
  310 N3 = IB(IA)+JB(LA)+KB(JA)+LB(KA)+NORG3
C
  400 CONTINUE
      GHONDO(N1) = ZERO
      GHONDO(N2) = ZERO
      GHONDO(N3) = ZERO
  800 CONTINUE
  820 CONTINUE
  840 CONTINUE
  860 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------
C*MODULE INT2B   *DECK GAMGEN
      SUBROUTINE xGAMGEN(NDER)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         ----- COMPUTES AND TABULATES F0(X) TO F5(X)
C               IN RANGE: 0.0 <= X <= CF1, IN STEPS OF FIN,
C               OR FOR GRADIENTS
C               IN RANGE: 0.0 <= X <= CG1, IN STEPS OF GIN.
C               USED BY THE 1-E ROUTINE SP0S1S IN MODULE INT2B
C               AND  BY THE 2-E ROUTINE JKDG80 IN MODULE GRD2B
C               THE TABLE IS GENERATED ONLY ONCE.
C               FINDS MAXIMUM VALUE OF S AND P COEFFICIENTS.
C               SETS LIMITS DETERMINING HOW ACCURATELY A SET
C               OF INTEGRALS NEED BE EVALUATED TO GUARANTEE
C               OVERALL INTEGRAL ACCURACY -----
C
C     INPUT ARGUMENT:
C           NDER= 0 ... INTEGRAL CALCULATIONS
C           NDER= 1 ... GRADIENT CALCULATIONS
C     (NEW) NDER= 2 ... INTEGRAL CALCULATIONS + FILMAX IN-LINING
C
C===================
C     JMS+MWS, 2001:
C     COMPLETELY REWRITTEN CODE GENERATED IDENTICAL RESULTS (ENERGY
C     AND GRADIENT AND NUMBERS OF INTEGRALS, USING THE 1988 PARAMETERS.
C     THE NEW INTEGRAL CODE MAY BE SLIGHTLY SLOWER THAN THE ORIGINAL,
C     AND THE GRADIENT ABOUT THE SAME SPEED AS THE ORIGINAL, BUT THIS
C     DEPENDS ON EXACTLY WHAT COMPUTER BRAND IS BEING USED.
C
C     INVESTIGATION OF THE ACCURACY PARAMETERS LEAD TO SETTING
C     'EI1'='EI2'=1.E-11, INSTEAD OF 1.E-10, AND IN A FORM WHERE THEY
C     ARE TIED TO THE INPUT VALUE OF CUTOFF. THIS CHANGE IS RESPONSIBLE
C     FOR FINDING SOME OF THE "MISSING INTEGRALS", BRINGING THE ROTATED
C     AXIS CODE'S INTEGRAL COUNT CLOSE TO THE HONDO COUNT.
C     SETTING 'EI1'='EI2'=1.E-12 MAKES THE NUMBER OF INTEGRALS ALMOST
C     IDENTICAL TO THE HONDO COUNT, BUT TAKES TOO MUCH LONGER.
C     FURTHER, SETTING 'EI1'=1.E-11 AND 'EI2'=1.E-13 GETS THE 'CORRECT'
C     NUMBER OF INTEGRALS WITH NOT APPRECIABLE PENALTY IN RUNNING TIME.
C
C     CHANGING 'CUF' FROM 20 TO 25 REMOVES HALF OF THE ENERGY DIFFERENCE
C     WITH HONDO, FURTHER INCREASING 'CUF' DOES NOT SEEM TO REMOVE ANY
C     MORE OF THE ENERGY DIFFERENCE.
C     SETTING 'CUG' NEAR 17, GETS AN RMS GRADIENT CLOSER TO HONDO'S.
C     SETTING 22 <= 'CX1' <= 28 DOES NOT PARTICULARLY AFFECT ACCURACY
C     (IMPLIES AN EXPONENTIAL THRESHOLD OF 50.657 <= 'CUX' <= 64.472).
C
C     CHANGING THE INTERVALS 'FIN' AND 'GIN' FOR THE ERF INTERPOLATION
C     TABLE MAY BE THE MOST CONSTRUCTIVE WAY TO PRODUCE BETTER AGREEMENT
C     WITH THE HONDO INTEGRALS.  BASIS SETS WITH DIFFUSE L SHELLS MAY BE
C     THE PLACE TO EXPLORE THIS, AS EVEN THE NEW VALUES FOR 'EI1/2' AND
C     'CUF' DO NOT PRODUCE ACCURATE VIRTUAL ORBITAL LCAO VALUES.
C     HOWEVER, THE ZTHRSH TEST IN THE MAIN INTEGRAL DRIVER WAS CHANGED
C     TO SEND DIFFUSE L SHELL INTEGRALS TO HONDO ONLY IF THE VIRTUALS
C     ARE OCCUPIED (E.G. MP2 AND CI).
C
C     THE FOLLOWING PARAMETERS ARE THE KEY ONES IN DETERMINING
C     THE ACCURACY OF THE INTEGRALS CALCULATION:
C
C     EI1 SUBSTITUTES ERROR1 DEFINED IN ROUTINE FILMAX IN MODULE INT2B
C     EI2 SUBSTITUTES ERROR2 DEFINED IN ROUTINE FILMAX IN MODULE INT2B
C     FIN SUBSTITUTES P08    INTERVAL FOR INTEGRAL CALCULATION
C     GIN SUBSTITUTES P05    INTERVAL FOR GRADIENT CALCULATION
C     CUF SUBSTITUTES VAR1   DEFINED IN ROUTINE FILMAX IN MODULE INT2B
C     CUF SUBSTITUTES VAR2   DEFINED IN ROUTINE FILMAX IN MODULE INT2B
C     CUG SUBSTITUTES VAR3   THRESHOLD IN GRADIENT CALCULATION
C
C     SENSITIVITY STUDIES SHOW THAT:
C     EI1 AND EI2 ARE IMPORTANT TO GET THE NUMBER OF INTEGRALS CORRECT
C     FIN AND GIN ARE IMPORTANT TO GET A GOOD ESTIMATE OF THE TOTAL
C                 ENERGY AND THE RMS GRADIENT, RESPECTIVELY
C     CUF IS IMPORTANT TO GET A GOOD ESTIMATE OF THE TOTAL ENERGY
C     CUG IS IMPORTANT TO GET A GOOD ESTIMATE OF THE RMS GRADIENT
C     CUX IS NOT PARTICULARLY IMPORTANT TO GET A GOOD ESTIMATE OF THE
C                 TOTAL ENERGY OR THE RMS GRADIENT
C
C     THE FOLLOWING PARAMETERS VALUES SEEM TO PRODUCE A GOOD BALANCE
C     BETWEEN ACCURACY AND SPEED OF THE CALCULATIONS
C
      PARAMETER (EC1= 1.0D-02, EC2= 1.0D-04)
      PARAMETER (FIN=6.25D-02, GIN=0.05D+00)
      PARAMETER (CF1=25.0D+00, CG1=17.0D+00, CX1=25.0D+00)
C
C     AS AN INDICATION OF THE ABOVE, SEE THE FOLLOWING INSTANCES:
C       6-31G(D) BASIS SET    CPU TIMING FROM SUN E450             TOTAL
C     INTTYP       # INTS       TOTAL ENERGY        RMS GRAD.       TIME
C PENTENAL MOLECULE...106 AO-S
C     POPLE      11,833,060   -268.827 020 1639   0.009 403 622    125.5
C     HONDO      11,833,060   -268.827 020 1663   0.009 403 652    162.9
C NICOTINE MOLECULE...208 AO-S
C     POPLE     128,845,447   -495.693 471 9566   0.013 383 774   1524.5
C     HONDO     128,845,447   -495.693 471 9600   0.013 383 734   1809.4
C LUCIFERIN MOLECULE...294 AO-S
C     POPLE     249,632,366  -1549.690 097 1864   0.024 472 402   3528.2
C     HONDO     249,632,365  -1549.690 097 1957   0.024 472 438   3996.6
C
C     HP'S F90 AND LINUX' G77 DO NOT LIKE THE ARITHMETIC BELOW.
C     THE SOLUTION IS TO COMPUTE M3R FROM THESE FORMULAE,
C     AND COPY ITS VALUE TO GAMGEN, JKDG80, AND SP0S1S
C
C     PARAMETER (M3F=INT(CF1/FIN))
C     PARAMETER (M3G=INT(CG1/GIN))
C     PARAMETER (M3R=MAX(M3F,M3G)+1)
C
C     ==================================================================
C     THIS GROUP OF STATEMENTS IS INTERRELATED FOR ROUTINES
C     GAMGEN, GENR70 AND SP0S1S IN INT2B, AND JKDG80, GRDG80 AND GRDPER
C     IN GRD2B. IF YOU CHANGE THESE, CHANGE ALL SIX PLACES EQUIVALENTLY.
C
      PARAMETER (M3R=401)
      COMMON /xINTACC/ EI1,EI2,CUF,RFI,CUG,RGI,CUX,R(3,0:5,0:M3R)
C     ==================================================================
C
C     SETTING 15 <= 'CUT' <= 25 DOES NOT PARTICULARLY AFFECT ACCURACY.
C     GIVEN THAT (GENERALLY) WE WORK IN 64-BITS PRECISION (APPROX.
C     EQUAL TO 14 SIGNIFICANT DIGITS), A RELATIVE TOLERANCE 'RTL'
C     SMALLER THAN 1.E-14 (SEE EXPANSIONS IN LOOPS 110 AND 140)
C     WOULD NOT MAKE ANY SENSE.
C
      PARAMETER (CUT=20.0D+00, RTL=1.0D-14)
C===================
C
#include "gamess.par"
C
      COMMON /xMAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     2                CF(MXGTOT),CG(MXGTOT),
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      LOGICAL OUT
      COMMON /xSHLT  / TOL,CUTOFF,OUT,ICOUNT
C
      DIMENSION F(0:5),X(-2:M3R+2),E(-2:M3R+2),Y(0:5,-2:M3R+2)
C
      PARAMETER (PT5=0.5D+00, ONE=1.0D+00, RLN10=2.30258D+00)
      PARAMETER (ZPT184=0.184D+00, OO6=1.0D+00/6.0D+00)
CC    PARAMETER (PIDIV4=0.785398163397448309615660845819D+00)
      PARAMETER (PIDIV4=0.785398163397448D+00)
C
CC    PIDIV4= ATAN(ONE)
      PI4= PIDIV4
C
C THIS IS A CHECK ON THE DIMENSIONING
C
      M3F=INT(CF1/FIN)
      M3G=INT(CG1/GIN)
      IF(M3R .LT. MAX(M3F,M3G)+1) THEN
         WRITE(6,*)'GAMGEN: DIMENSION PROBLEM, M3F,M3G,M3R=',M3F,M3G,M3R
cwc         CALL ABRT
      END IF
C
      XIN= FIN
      IF(NDER.EQ.1) XIN= GIN
C
C PREPARE THE ARGUMENT AND SAVE A FEW OPERATIONS.
C
      EMI= EXP(-XIN)
      XVA=-XIN-XIN
      EMX= PT5* EXP(-XVA)
      DO 100 K=-2,M3R+2
CC       X(K)= XIN*K
CC       E(K)= PT5* EXP(-X(K))
         X(K)= XVA
         E(K)= EMX
         XVA= XVA+XIN
         EMX= EMX*EMI
  100 CONTINUE
C
      M= 5
      REM= M
      DO 180 K=-2,M3R+2
         XVA= X(K)
         EMX= E(K)
         FAC= REM-PT5
C
C RANGES FOR APPROXIMATE SOLUTION
C
         IF(XVA.LE.CUT) THEN
            DEN= FAC+ONE
            TRM= ONE/DEN
            SUM= TRM
  110       DEN= DEN+ONE
               TRM= TRM*XVA/DEN
               SUM= SUM+TRM
            IF( ABS(TRM).GT. ABS(RTL*SUM)) GO TO 110
C
  120       F(M)= EMX*SUM
         ELSE
            XIN= ONE/XVA
            EMI= EMX*XIN
            FAX= FAC*XIN
            APP= SQRT(PI4*XIN)
            FAC= FAX
            DO 130 I=1,M
               APP= APP*FAC
  130       FAC= FAC-XIN
            TRM= XVA+REM
            NTR= INT(TRM)
            TRM= ONE
            SUM= ONE-APP/EMI
            DO 140 I=2,NTR
               TRM= TRM*FAX
               SUM= SUM+TRM
               IF( ABS(TRM).LE. ABS(RTL*SUM)) GO TO 150
  140       FAX= FAX-XIN
C
  150       F(M)=-EMI*SUM
         ENDIF
C
C RECUR DOWNWARDS TO F(0)
C
         DEN= REM-PT5
         DO 160 J=M,1,-1
            F(J-1)=(F(J)*XVA+EMX)/DEN
  160    DEN= DEN-ONE
C
         DO 170 J=0,M
            Y(J,K)= F(J)
  170    CONTINUE
  180 CONTINUE
C
C NOW, PUT VARIABLES INTO COMMON BLOCK /INTACC/
C
      IF(NDER.EQ.2) THEN
C
C FIND MAXIMUM VALUE OF S AND P COEFFICIENTS
C
         DO 195 I=1,NSHELL
            L = KSTART(I)
            N = L+KNG(I)-1
            DO 190 J=L,N
  190       CMAX(J)= MAX( ABS(CS(J)), ABS(CP(J)))
  195    CONTINUE
         EI1= EC1*CUTOFF
         EI2= EC2*CUTOFF
      ENDIF
      CUF= CF1
      RFI= ONE/FIN
      CUG= CG1
      RGI= ONE/GIN
      CUX= CX1*RLN10
C
      C01 = OO6+ZPT184*PT5
      C02 = OO6*ZPT184
      DO 200 K=0,M3R
         DO 200 J=0,M
            R(1,J,K)= Y(J,K  )
            R(2,J,K)= Y(J,K+1)-Y(J,K  )
            T01 =(Y(J,K+1)-Y(J,K  ))-(Y(J,K  )-Y(J,K-1))
            T02 =(Y(J,K+2)-Y(J,K+1))-(Y(J,K-1)-Y(J,K-2))
            R(3,J,K)= C01*T01- C02*T02
  200 CONTINUE
C
      RETURN
      END
c----------------------- INTER.f end -----------------------------------
