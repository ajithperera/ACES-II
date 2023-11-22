C*MODULE GRD2A   *DECK JKDER
      SUBROUTINE xJKDER
C
C       THE DRIVER FOR THE TWO ELECTRON GRADIENT
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION METHMC
C
cwc1
#undef _PARALLEL_ACES_COM_
#include "parallel_aces.com"
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
cwc0
      LOGICAL SKIPI,SKIPJ,SKIPK,SKIPL,NOCORE,NOOPEN
      LOGICAL HFSCF,UHFTYP,ROGVB,PACK2E,POPLE,HONDO
      LOGICAL GOPARR,DSKWRK,MASWRK,NXT
      LOGICAL SOME,OUT,DBG,MP2,MC,PER1,PER2,PER3,UMP2
      LOGICAL CANONC,FCORE,FORS,NOCI,EKT,LINSER,LCIS
C
#include "gamess.par"
      PARAMETER (MXRT=100, MXNORO=250)
      PARAMETER (MXPT=100, MXFRG=50, MXFGPT=MXPT*MXFRG)
C
      DIMENSION M0(48),M1(48),M2(48),M3(48)
C
      COMMON /xDERSKP/ IIAT,JJAT,KKAT,LLAT,SKIPI,SKIPJ,SKIPK,SKIPL
cwc      COMMON /xDETWFN/ WSTATE(MXRT),SPINS(MXRT),CRIT,PRTTOL,SDET,SZDET,
cwc     *                GRPDET,STSYM,GLIST,
cwc     *                NFLGDM(MXRT),IWTS(MXRT),NCORSV,NCOR,NACT,NORBDT,
cwc     *                NADET,NBDET,KDET,KSTDET,IROOT,IPURES,MAXW1,NITER,
cwc     *                MAXP,NCIDET,IGPDET,KSTSYM,NFTGCI
cwc      COMMON /xDFGRID/DFTTHR,SW0,NDFTFG,NRAD,NTHE,NPHI,NRAD0,NTHE0,NPHI0
cwc      COMMON /xDFTPAR/ DFTTYP(20),EXENA,EXENB,EXENC,IDFT34,NAUXFUN,
cwc     *                                                    NAUXSHL
cwc      COMMON /xDLT   / LAT,LBT,LCT,LDT
      COMMON /xDSHLNO/ LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     *                NIJ,IJD,KLD,IJ,KL
      COMMON /xDSHLT / RTOL,DTOL,VTOL1,VTOL2,VTOLS,OUT,DBG
      COMMON /xFMCOM / X(1)
cwc      COMMON /xFRGINF/ NMPTS(MXFRG),NMTTPT,IEFC,IEFD,IEFQ,IEFO,
cwc     *                NPPTS(MXFRG),NPTTPT,IEFP,
cwc     *                NRPTS(MXFRG),NRTTPT,IREP,ICHGP,NFRG
      COMMON /xGRAD  / DE(3,MXATM)
cwc1
c buffer for mpi_allreduce
      dimension buf(3*MXATM)
cwc0
cwc      COMMON /xGUGWFN/ NFZC,NMCC,NDOC,NAOS,NBOS,NALP,NVAL,NNEXT,NFZV,
cwc     *                IFORS,IEXCIT,ICICI,NOIRR
      COMMON /xIJPAIR/ IA(MXAO)
      COMMON /xINFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM)
      COMMON /xINTFIL/ NINTMX,NHEX,NTUPL,PACK2E,IPOPLE
      COMMON /xINTOPT/ ISCHWZ,IECP,NECP,IEFLD
      COMMON /xIOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
cwc      COMMON /xMCINP / METHMC,CISTEP,ACURCY,ENGTOL,DAMP,MICIT,NWRDMC,
cwc     *                NORBMC,NOROT(2,MXNORO),MOFRZ(15),NPFLG(10),
cwc     *                NOFO,CANONC,FCORE,FORS,NOCI,EKT,LINSER
cwc      COMMON /xMP2DM2/ NOC1,NOC2,ISTEP,NOC3,NOC4
cwc      COMMON /xMP2PAR/ OSPT,TOL,METHOD,NWDMP2,MPPROP,NACORE,NBCORE,
cwc     *                NOA,NOB,NORB,NBF,NOMIT,MOCPHF,MAXITC
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /xOUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /xPAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
cwc      COMMON /xPAULMO/ CENTNM(MXFGPT),CENTCD(3,MXFGPT),NORBEF(MXFRG),
cwc     *                NPBF(MXFRG),NTMO
      COMMON /xRESTAR/ TIMLIM,IREST,NREC,INTLOC,IST,JST,KST,LST
cwc      COMMON /xRUNOPT/ RUNTYP,EXETYP,NEVALS
      COMMON /xSHLBAS/ MAXTYP,MAXNUM
cwc      COMMON /xSIMDAT/ NACC,NREJ,IGOMIN,NRPA,IBWM,NACCT,NREJT,NRPAT,
cwc     *                NPRTGO,IDPUNC,IGOFLG
      COMMON /xSYMTRY/ MAPSHL(MXSH,48),MAPCTR(MXATM,48),
     *                T(432),INVT(48),NT
cwc      COMMON /xTMVALS/ TI,TX,TIM
      COMMON /xWFNOPT/ SCFTYP,CITYP,DFTYPE,CCTYP,MPLEVL,MPCTYP
C
      COMMON/xDERMEM/IWFN,IXCH,INIJG,IGINT,IFINT,ISINT,IIJKLG,
     1 IDAB,ICHRG,IXY,IXZ,IYZ,IX,IY,IZ,ISJ,ISK,ISL,IGIJKL,IGNKL,IGNM,
     2 IDIJ,IDKL,IB00,IB01,IB10,IC00,ID00,IF00,
     3 IAAI,IAAJ,IBBK,IBBL,IFI,IFJ,IFK,IFL,
     4 ISII,ISJJ,ISKK,ISLL,ISIJ,ISIK,ISIL,ISJK,ISJL,ISKL,
     5 IDIJSI,IDIJSJ,IDKLSK,IDKLSL,IABV,ICV,IRW
      COMMON/xINDD80/IMAX,JMAX,KKKMAX,LMAX
C
      PARAMETER (RLN10=2.30258D+00)
      PARAMETER (TEN=10.0D+00, ONE=1.0D+00)
      PARAMETER (TENM9=1.0D-09, TENM11=1.0D-11)
      PARAMETER (TENM20=1.0D-20, PT5=0.5D+00, TENM12=1.0D-12)
C
      DIMENSION LENSHL(5)
      DATA LENSHL/1,4,10,20,35/
      DATA CHECK,GRD2,DEBUG/8HCHECK   ,8HGRD2    ,8HDEBUG   /
      DATA RHF,UHF,ROHF,GVB
     *    /8HRHF     ,8HUHF     ,8HROHF    ,8HGVB     /
      DATA RMC/8HMCSCF    /
      DATA ALDET,GUGA,GENCI,ORMAS
     *       /8HALDET   ,8HGUGA    ,8HGENCI   ,8HORMAS   /
      DATA CIS/8HCIS     /
C
C     ----- THIS IS THE MAIN 2E- GRADIENT DRIVER -----
C
cwc      DBG = EXETYP.EQ.DEBUG
cwc      OUT = EXETYP.EQ.GRD2.OR.NPRINT.EQ.-4
cwc1
      dbg=.false.
      out=.false.
cwc0
      SOME = MASWRK  .AND.  NPRINT.NE.-5
cwc      IF(SOME) WRITE (IW,9008)
C
C        INITIALIZE PARALLEL
C
      NXT = IBTYP.EQ.1
      IPCOUNT = ME - 1
      NEXT = -1
      MINE = -1
cwc1
      nloop=-1
cwc0
C
C     ----- SET STARTING PARAMETERS -----
C
      HONDO = .TRUE.
      POPLE = IPOPLE.EQ.1
C
C      CUTOFF IS THE SCHWARZ SCREENING CUT OFF
C      DABCUT IS THE TWO PARTICLE DENSITY CUT OFF
C
      CUTOFF=TENM9
      IF(.NOT.POPLE) CUTOFF=CUTOFF/TEN
      CUTOFF2=CUTOFF/2.0D+00
C
      ZBIG = 0.0D+00
      DO ISH=1,NSHELL
         I1=KSTART(ISH)
         I2=I1+KNG(ISH)-1
         DO IG=I1,I2
            IF(EX(IG).GT.ZBIG) ZBIG = EX(IG)
         ENDDO
      ENDDO
      DABCUT=TENM11
      IF(ZBIG.GT.1.0D+06) DABCUT = DABCUT/TEN
      IF(ZBIG.GT.1.0D+07) DABCUT = DABCUT/TEN
C
C      VTOLS ARE CUT OFFS USED BY THE POPLE PACKAGE
C      CURRENT VALUES ARE FROM HONDO 8, SEE G92 FOR OTHER POSSIBILITIES
C
      VTOL1 = TENM12
      VTOL2 = TENM12
      VTOLS = TENM20
      DTOL = TEN**(-ITOL)
      RTOL = RLN10*ITOL
C
C      INITIALIZE THE INTEGRAL BLOCK COUNTERS TO ZERO
C
      IISKIP1= 0
      IISKIP2= 0
      IDID = 0
C
      LCIS = .FALSE.
cwc      IF(CITYP.EQ.CIS) LCIS = .TRUE.
C
      MC    = SCFTYP.EQ.RMC
      MP2   = SCFTYP.EQ.RHF .AND.  MPLEVL.EQ.2
      UMP2  = SCFTYP.EQ.UHF .AND.  MPLEVL.EQ.2
      HFSCF = MPLEVL.EQ.0   .AND.   CITYP.NE.CIS .AND.
     *        (SCFTYP.EQ.RHF .OR. SCFTYP.EQ.UHF .OR. SCFTYP.EQ.ROHF)
      UHFTYP=                     SCFTYP.EQ.UHF .OR. SCFTYP.EQ.ROHF
      ROGVB = SCFTYP.EQ.GVB
C
C     ----- SET POINTERS FOR PARTITIONING OF MEMORY -----
C
      L1 = NUM
      L2 = (NUM*NUM+NUM)/2
      L3 = NUM*NUM
      NSH2=(NSHELL*NSHELL+NSHELL)/2
C
      DO 100 I = 1,NUM
         IA(I) = (I*I-I)/2
  100 CONTINUE
C
C     ----- READ IN 1E-GRADIENT -----
C
      CALL yDAREAD(IDAF,IODA,DE,3*NAT,3,0)
cwc      IF (GOPARR) CALL XSCAL(3*NAT,ONE/NPROC,DE,1)
cwc1
      call XSCAL (3*NAT,one/nprocs,DE,1)
cwc0
C
C     ----- GRID-DFT EXCHANGE-CORRELATION -----
C     -----    DERIVATIVE CONTRIBUTION    -----
C
cwc      IF (NDFTFG.EQ.1) CALL DFTDER
C
C              CALCULATE THE LARGEST SHELL TYPE
C
      CALL xBASCHK(MAXTYP)
      MAXSHL = LENSHL(MAXTYP+1)
C              DO AT LEAST AN L SHELL
      IF (MAXSHL.LT.4) MAXSHL=4
C
C       IF WE ARE USING THE POPLE PACKAGE AND DO NOT HAVE ANY SHELLS
C       LARGER THAN AN L-SHELL THEN SKIP THE SETUP FOR THE RYS PACKAGE
C
      IF (POPLE.AND.MAXTYP.LT.2) HONDO = .FALSE.
C
C      CALCULATE THE NUMBER OF MCSCF CORE AND ACTIVE ORBITALS
C
      NCI = 0
cwc      IF(MC) THEN
cwc         IF(CISTEP.EQ.ALDET) NCI = NACT
cwc         IF(CISTEP.EQ.GENCI) NCI = NACT
cwc         IF(CISTEP.EQ.ORMAS) NCI = NACT
cwc         IF(CISTEP.EQ.GUGA)  NCI = NDOC+NAOS+NBOS+NALP+NVAL+NNEXT
cwc         IF(CISTEP.EQ.ALDET) NMCCOR = NCORSV
cwc         IF(CISTEP.EQ.GENCI) NMCCOR = NCORSV
cwc         IF(CISTEP.EQ.ORMAS) NMCCOR = NCORSV
cwc         IF(CISTEP.EQ.GUGA)  NMCCOR = NMCC
cwc      END IF
C
C     FIGURE OUT THE MEMORY WE NEED FOR STORING DENSITY MATRIX
C     AND OTHER WAVEFUNCTION INFORMATION. -JKDMEM- ALLOCATES
C     MEMORY FOR DERIVATIVE COMPUTATION AND 2ND ORDER DENSITY
C     AFTER -LENGTH- WORDS.
C
      LENRHF = L2
      LENUHF = L2+L2
cwc      LENGVB = L2+L3
cwc      LENMC  = L2+L2+L3+(NCI*NCI*NCI*NCI)
cwc      LENMC1 = NINTMX+NINTMX+(NCI*NCI+NCI)/2
cwc      LENMC2 = MAXSHL*NCI   *NCI   *NCI
cwc     *       + MAXSHL*MAXSHL*NCI   *NCI
cwc     *       + MAXSHL*MAXSHL*MAXSHL*NCI
cwc     *       + MAXSHL*MAXSHL*MAXSHL*MAXSHL
cwc      LENMC  = LENMC+MAX(LENMC1,LENMC2)
      LENMP  = L2+L2+L3
      LENUMP = 2*LENMP
      LENCIS = L2+L2+L2+L3
C
                  LENGTH=LENRHF
      IF (UHFTYP) LENGTH=LENUHF
cwc      IF (ROGVB)  LENGTH=LENGVB
cwc      IF (MC)     LENGTH=LENMC
      IF (MP2)    LENGTH=LENMP
      IF (UMP2)   LENGTH=LENUMP
cwc      IF (LCIS)   LENGTH=LENCIS
cwc      IF (NTMO.GT.0) LENGTH = LENGTH + L2
C
C       CALCULATE THE AMOUNT OF MEMORY NEEDED AND SET THE POINTERS
C       FOR BOTH PACKAGES
C
      CALL yVALFM(LOADFM)
      CALL xJKDMEM(1,LOADFM,IADDR,LENGTH,MINXYZ,MAXXYZ,MINVEC,POPLE,
     *            MP2.OR.UMP2)
C
C     ----- CARRY OUT SET UP TASKS -----
C
      NOCORE = .FALSE.
      NOOPEN = .FALSE.
C
C          MEMORY ALLOCATION FOR MCSCF FUNCTIONS
C
cwc      IF (MC) THEN
cwc         NEED = IADDR - LOADFM -1
cwc         CALL yGETFM(NEED)
cwcC
cwc         K10  = IWFN
cwc         K20  = K10 + L2
cwc         K30  = K20 + L2
cwc         K40  = K30 + L3
cwc         K51  = K40 + NCI*NCI*NCI*NCI
cwc         K61  = K51 + (NCI*NCI+NCI)/2
cwc         K71  = K61 + NINTMX
cwc         LAST = K71 + NINTMX
cwcC
cwc         IF (EXETYP .EQ. CHECK) GO TO 600
cwcC
cwcC     ----- READ MCSCF WAVEFUNCTION INFORMATION -----
cwcC
cwc         CALL JKDMC(X(K10),X(K20),X(K30),X(K40),
cwc     *              X(K51),X(K61),X(K71),IA,NUM,NMCCOR,NCI)
cwcC
cwc         K50  = K40 + NCI   *NCI   *NCI   *NCI
cwc         K60  = K50 + MAXSHL*NCI   *NCI   *NCI
cwc         K70  = K60 + MAXSHL*MAXSHL*NCI   *NCI
cwc         K80  = K70 + MAXSHL*MAXSHL*MAXSHL*NCI
cwc         LAST = K80 + MAXSHL*MAXSHL*MAXSHL*MAXSHL
cwc         CALL VCLR(X(K50),1,MAXSHL*NCI   *NCI   *NCI   )
cwc         CALL VCLR(X(K60),1,MAXSHL*MAXSHL*NCI   *NCI   )
cwc         CALL VCLR(X(K70),1,MAXSHL*MAXSHL*MAXSHL*NCI   )
cwc         CALL VCLR(X(K80),1,MAXSHL*MAXSHL*MAXSHL*MAXSHL)
cwc         CALL BASCHK(MAXANG)
cwc         NSHLMX=((MAXANG+1)*(MAXANG+2))/2
cwc         IF(NSHLMX.LE.4) NSHLMX=4
cwc         NGIJKL = NSHLMX**4
cwc         CALL VCLR(X(IDAB),1,NGIJKL)
cwc      END IF
C
C     ----- CLOSED SHELL -MP2- WAVEFUNCTIONS -----
C
cwc      IF (.NOT.MP2) GO TO 355
cwcC
cwc      NOASAV=NOA
cwc      NOA=NOA-NACORE
cwc      NOC1 = 1
cwc      ISTEP= 0
cwc      LPHF   = IWFN
cwc      LPMP2  = LPHF  + L2
cwc      LVEC   = LPMP2 + L2
cwc      LAST   = LVEC  + L3
cwc      LMPDM2 = IADDR
cwcC
cwc      NEED1= IADDR - LOADFM -1
cwc      CALL GETFM(NEED1)
cwcC
cwcC        THE FINAL ARRAY IS THE MP2 DENSITY, WHICH IS PROPORTIONAL
cwcC        TO N**4 WORDS.  THIS MUST BE SLICED UP INTO PASSES OVER
cwcC        THE OCCUPIED ORBITALS, CONTAINING -NUMCO- EACH.
cwcC
cwc      CALL VALFM(LOADFM)
cwc      CALL GOTFM(NGOTMX)
cwc      NMIN = L2*L1*1
cwc      IF(NMIN.GT.NGOTMX) THEN
cwc         WRITE(IW,9997) NMIN-NGOTMX
cwc         CALL ABRT
cwc      END IF
cwc      NUMCO=MIN(NOA,NGOTMX/(L2*L1))
cwcC
cwc      NMIN  = NEED1 + L2*L1*1
cwc      NMAX  = NEED1 + L2*L1*NOA
cwc      NUSED = NEED1 + L2*L1*NUMCO
cwc      NPASS = (NOA-1)/NUMCO + 1
cwc      IF(SOME) WRITE(IW,9998) NOA,NMIN,NMAX,NUSED,NPASS,NUMCO
cwc      IF(NPASS.GT.1) THEN
cwc         NINCR=L2*L1
cwc         IF(SOME) WRITE(IW,9996) NINCR
cwc      END IF
cwcC
cwc      CALL VALFM(LOADFM)
cwc      LAST = LMPDM2 + L2*L1*NUMCO
cwc      NEED2= LAST - LOADFM -1
cwc      CALL GETFM(NEED2)
cwc      NEED = NEED1 + NEED2
cwc      IF(EXETYP.EQ.CHECK) GO TO 600
cwcC
cwcC     ----- READ WAVEFUNCTION INFORMATION -----
cwcC     THIS POINT IS RETURNED TO FOR EACH PASS OVER OCCUPIED MO-S
cwcC
cwc  350 CONTINUE
cwc      ISTEP = ISTEP+1
cwc      NOC2 = NOC1+NUMCO-1
cwc      NOC2 = MIN(NOA,NOC2)
cwc      CALL JKDMP2(X(LPHF),X(LPMP2),X(LVEC),X(LMPDM2),L3,L2,L1)
cwcC
cwc      GO TO 375
cwcC
cwc  355 CONTINUE
cwcC
cwcC     ----- OPEN SHELL -MP2- WAVEFUNCTIONS -----
cwcC
cwc      IF (.NOT.UMP2) GO TO 365
cwcC
cwc      NOCA=NOA-NACORE
cwc      NOCB=NOB-NBCORE
cwc      NOC1 = 1
cwc      NOC3 = 1
cwc      ISTEP= 0
cwc      LPHFA  = IWFN
cwc      LPHFB  = LPHFA  + L2
cwc      LPMP2A = LPHFB  + L2
cwc      LPMP2B = LPMP2A + L2
cwc      LVECA  = LPMP2B + L2
cwc      LVECB  = LVECA  + L3
cwc      LAST   = LVECB  + L3
cwc      LMPDM2 = IADDR
cwcC
cwc      NEED1= IADDR - LOADFM -1
cwc      CALL GETFM(NEED1)
cwcC
cwcC        LMPDM2 CONTAINS THE MP2 2 E- DENSITY, WHICH IS PROPORTIONAL
cwcC        TO N**4 WORDS.  THIS MUST BE SLICED UP INTO PASSES OVER
cwcC        THE OCCUPIED ORBITALS, CONTAINING -NUMCO- EACH.
cwcC
cwc      CALL VALFM(LOADFM)
cwc      CALL GOTFM(NGOTMX)
cwc      NMIN = L2*L1*2
cwc      IF(NMIN.GT.NGOTMX) THEN
cwc         WRITE(IW,9997) NMIN-NGOTMX
cwc         CALL ABRT
cwc      END IF
cwc      NOCX=MAX(NOCA,NOCB)
cwc      NUMCO=MIN(NOCX,NGOTMX/(L2*L1*2))
cwcC
cwc      NMIN  = NEED1 + L2*L1*2
cwc      NMAX  = NEED1 + L2*L1*2*NOCX
cwc      NUSED = NEED1 + L2*L1*2*NUMCO
cwc      NPASS = (NOCX-1)/NUMCO + 1
cwc      IF(SOME) WRITE(IW,9995) NOCA,NOCB,NMIN,NMAX,NUSED,NPASS,NUMCO
cwc      IF(NPASS.GT.1) THEN
cwc         NINCR=L2*L1*2
cwc         IF(SOME) WRITE(IW,9996) NINCR
cwc      END IF
cwcC
cwc      CALL VALFM(LOADFM)
cwc      LMPDM2B = LMPDM2 + L2*L1*NUMCO
cwc      LAST    = LMPDM2B+ L2*L1*NUMCO
cwc      NEED2   = LAST - LOADFM -1
cwc      CALL GETFM(NEED2)
cwc      NEED    = NEED1 + NEED2
cwc      IF(EXETYP.EQ.CHECK) GO TO 600
cwcC
cwcC     ----- READ WAVEFUNCTION INFORMATION -----
cwcC     THIS POINT IS RETURNED TO FOR EACH PASS OVER OCCUPIED MO-S
cwcC
cwc  360 CONTINUE
cwc      ISTEP = ISTEP+1
cwc      NOC2  = MIN(NOCA,NOC1+NUMCO-1)
cwc      NOC4  = MIN(NOCB,NOC3+NUMCO-1)
cwcC
cwc      CALL JKDUMP(X(LPHFA),X(LPMP2A),X(LVECA),X(LPHFB),X(LPMP2B),
cwc     *            X(LVECB),X(LMPDM2),X(LMPDM2B),L3,L2,L1)
cwcC
cwc      GO TO 375
cwcC
cwc  365 CONTINUE
cwcC
cwcC     ----- CLOSED SHELL CIS WAVEFUNCTIONS -----
cwcC
cwc      IF (.NOT.LCIS) GO TO 375
cwcC
cwc      NOASAV=NOA
cwc      NOA=NOA-NACORE
cwc      NOC1 = 1
cwc      ISTEP= 0
cwc      LPHF   = IWFN
cwc      LPCIS  = LPHF  + L2
cwc      LTCIS  = LPCIS + L2
cwc      LAST   = LTCIS + L3
cwcC
cwc      NEED= LAST - LOADFM -1
cwc      CALL GETFM(NEED)
cwcC
cwc      IF(EXETYP.EQ.CHECK) GO TO 600
cwcC
cwcC     ----- READ WAVEFUNCTION INFORMATION -----
cwcC
cwc      ISTEP = ISTEP+1
cwc      NOC2 = NOC1+NUMCO-1
cwc      NOC2 = MIN(NOA,NOC2)
cwc      CALL JKDCIS(X(LPHF),X(LPCIS),X(LTCIS),L3,L2)
C
C     NOW HANDLE SCF (RHF, UHF, ROHF, GVB) WAVEFUNCTIONS
C
  375 CONTINUE
      IF (HFSCF  .OR.  ROGVB) THEN
         LDEN = IWFN
         LVEC = LDEN + L2
         LAST = LVEC + MAX(L2,L3)
C
C     ADD IN DFT MEMORY
C
cwc1
         nauxfun=0
cwc0
         L4=NUM+NAUXFUN
         L6=L4*L4
cwc         IF((IADDR-LOADFM).LT.(5*L6+L4+3*MXATM)  .AND.
cwc     *       DFTTYP(1).NE.0.0D+00) THEN
cwcC
cwcC          GENERALLY, THE MEMORY THAT WE NEED FOR DFT IS LESS THAN THAT
cwcC          NEEDED FOR NON-DFT RUNS, SO ALLOCATING IT IS NOT A PROBLEM
cwcC          UNLESS IT'S A RUN WITH AUX FUNCTIONS
cwcC
cwc           IADDR=LOADFM+5*L6+L4+3*MXATM
cwc           IDFT1=LOADFM+1
cwc           IDFT2=IDFT1+L6
cwc           IDFT3=IDFT2+L6
cwc           IDFT4=IDFT3+L6
cwc           IDFT5=IDFT4+L6
cwc           IDFT6=IDFT5+L6
cwc           IDFT7=IDFT6+L4
cwc         ELSE
           IDFT1=LOADFM+1
           IDFT2=IDFT1+L6
           IDFT3=IDFT2+L6
           IDFT4=IDFT3+L6
           IDFT5=IDFT4+L6
           IDFT6=IDFT5+L6
           IDFT7=IDFT6+L4
cwc         END IF
C
C        TOTAL USED FOR DFT = IDFT7+3*MXATM
C
         NEED=IADDR-LOADFM
         CALL yGETFM(NEED)
cwc         IF (EXETYP .EQ. CHECK) GO TO 600
C
C     ----- READ WAVEFUNCTION INFORMATION -----
C
         CALL xDDEBUT(IA,X(LDEN),X(LVEC),L1,L2,L3,NOCORE,NOOPEN,
     *               HFSCF,UHFTYP,ROGVB)
C
cwc         IF (NTMO.GT.0) THEN
cwc            CALL DAREAD(IDAF,IODA,X(LVEC),L2,79,0)
cwc         END IF
      END IF
C
C        READ IN THE EXCHANGE INTEGRALS FROM DISK. IF THEY WERE NOT
C        PREVIOUSLY COMPUTED, THEN JUST SET THE ARRAY TO ONE, WHICH
C        EFFECTIVELY DEACTIVATES THE SCHWARZ SCREENING
C
      IF(ISCHWZ.EQ.1) THEN
         CALL yDAREAD(IDAF,IODA,X(IXCH),NSH2,54,0)
      ELSE
         DO 400 I=0,NSH2-1
            X(IXCH+I) = ONE
  400    CONTINUE
      END IF
C
C     ----- GET SYMMETRY MAPPING OF SHELLS -----
C
      CALL xJKDSET
C
C     ----- PREPARE FOR USE OF G80 DERIVATIVE ROUTINES -----
C
      IF (POPLE) CALL xGAMGEN(1)
C
C        SET UP THE 1-ELECTRON CHARGE DISTRIBUTION
C
      IF (HONDO) CALL xOEDHND(X(INIJG),X(ICHRG))
C
C        SQUARE DTOL FOR USE IN JKDSPD
C
      DTOL = DTOL*DTOL
      NC=1
      LDF=1
C
C     ----- I SHELL -----
C
cwc1
      ist=1
cwc0
      DO 560 II = IST,NSHELL
         DO 160 IT = 1,NT
            ID = MAPSHL(II,IT)
            IF (ID .GT. II) GO TO 560
            M0(IT) = ID
  160    CONTINUE
cwc         IF(MC.AND.NCI.GT.0) CALL DABMCI(X(K30),X(K40),X(K50),NCI,
cwc     *                                   II,NUM,MAXSHL,NMCCOR)
C
C     ----- J SHELL -----
C
        J0 = JST
        DO 540 JJ = J0,II
          JST = 1
          DO 220 IT = 1,NT
           JD = MAPSHL(JJ,IT)
           IF (JD .GT. II) GO TO 540
           ID = M0(IT)
           IF (ID .GE. JD) GO TO 200
           ND = ID
           ID = JD
           JD = ND
  200      IF (ID .EQ. II .AND. JD .GT. JJ) GO TO 540
           M1(IT) = ID
           M2(IT) = JD
  220     CONTINUE
C
C     ----- GO PARALLEL! -----
C
          IF (NXT .AND. GOPARR) THEN
             MINE = MINE + 1
cwc             IF (MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
cwc             IF (NEXT.NE.MINE) GO TO 540
          END IF
C
cwc          IF (MC.AND.NCI.GT.0)
cwc     1         CALL DABMCJ(X(K30),X(K50),X(K60),NCI,
cwc     2                     II,JJ,NUM,MAXSHL,NMCCOR)
C
C        GET IJ CHARGE DISTRIBUTION
C        REALLY WE JUST SET THE POINTERS TO THE CHARGE DISTRIBUTION
C
          IF (HONDO) THEN
            IIJJ=IA(MAX0(II,JJ))+MIN0(II,JJ)
            CALL xOEDRD(X(INIJG),NIJ,NIJ0,IIJJ)
            IF(NIJ.EQ.0) GO TO 540
          END IF
C
C     ----- K SHELL -----
C
          K0 = KST
          DO 520 KK = K0,II
            KST = 1
            DO 260 IT = 1,NT
             KD = MAPSHL(KK,IT)
             IF (KD .GT. II) GO TO 520
             M3(IT) = KD
  260       CONTINUE
cwc            IF(MC.AND.NCI.GT.0)
cwc     1         CALL DABMCK(X(K30),X(K60),X(K70),NCI,
cwc     2                     II,JJ,KK,NUM,MAXSHL,NMCCOR)
C
C     ----- L SHELL -----
C
            L0 = LST
            MAXLL = KK
            IF (KK .EQ. II) MAXLL = JJ
            DO 500 LL = L0,MAXLL
              LST = 1
              N4 = 0
              DO 340 IT = 1,NT
               LD = MAPSHL(LL,IT)
               IF (LD .GT. II) GO TO 500
               KD = M3(IT)
               IF (KD .GE. LD) GO TO 300
               ND = KD
               KD = LD
               LD = ND
  300          ID = M1(IT)
               JD = M2(IT)
               IF (ID .NE. II .AND. KD .NE. II) GO TO 340
               IF (KD .LT. ID) GO TO 320
               IF (KD .EQ. ID .AND. LD .LE. JD) GO TO 320
               ND = ID
               ID = KD
               KD = ND
               ND = JD
               JD = LD
               LD = ND
  320          IF (JD .LT. JJ) GO TO 340
               IF (JD .GT. JJ) GO TO 500
               IF (KD .LT. KK) GO TO 340
               IF (KD .GT. KK) GO TO 500
               IF (LD .LT. LL) GO TO 340
               IF (LD .GT. LL) GO TO 500
               N4 = N4+1
  340         CONTINUE
C
C     ----- GO PARALLEL! -----
C
cwc1
cwc              IF ((.NOT.NXT) .AND. GOPARR) THEN
cwc                 IPCOUNT = IPCOUNT + 1
cwc                 IF (MOD(IPCOUNT,NPROC).NE.0) GO TO 500
cwc              END IF
              nloop=nloop+1
              nloop=mod(nloop,nprocs)
              if (nloop.ne.irank) goto 500
cwc0

C
C     ----- CALCULATE Q4 FACTOR FOR THIS GROUP OF SHELLS -----
C
              Q4 = NT
              Q4 = Q4 / N4
C
C     ----- DECIDE ON DERIVATIVE INTEGRAL METHOD -----
C     ANY PURELY SP SET OF SHELLS CAN BE DONE WITH THE FASTER
C     POPLE/SCHLEGEL ROTATION ALGORITHM.  INTEGRALS INVOLVING
C     D AND HIGHER FUNCTIONS MUST USE RYS POLYNOMIAL CODE.
C
              POPLE = .TRUE.
              IF(IPOPLE.EQ.0) POPLE=.FALSE.
              IF(KTYPE(II).GT.2) POPLE=.FALSE.
              IF(KTYPE(JJ).GT.2) POPLE=.FALSE.
              IF(KTYPE(KK).GT.2) POPLE=.FALSE.
              IF(KTYPE(LL).GT.2) POPLE=.FALSE.
C
C         IMPLEMENT INTEGRAL SCREENING HERE USING EXCHANGE INTEGRALS
C         SEE H.HORN, H.WEISS, M.HAESER, M.EHRIG, R.AHLRICHS
C             J.COMPUT.CHEM. 12, 1058-1064(1991)
C         REGARDING THE ESTIMATION FORMULA (31) THAT IS USED HERE.
C
              IJIJ=IA(MAX0(II,JJ))+MIN0(II,JJ)
              KLKL=IA(MAX0(KK,LL))+MIN0(KK,LL)
              GMAX=(X(IXCH+IJIJ-1)*X(IXCH+KLKL-1))
C
C                COARSE SCREENING, ON JUST THE INTEGRAL VALUE
C
              IF (GMAX.LT.CUTOFF) THEN
                 IISKIP1 = IISKIP1+1
                 GO TO 500
              END IF
C
           ISH=II
           JSH=JJ
           KSH=KK
           LSH=LL
C
           IF (POPLE) THEN
             AX1=PT5
             IF(ISH.NE.JSH) AX1=AX1+AX1
             IF(KSH.NE.LSH) AX1=AX1+AX1
             IF(ISH.NE.KSH.OR.JSH.NE.LSH) AX1=AX1+AX1
             Q4 = Q4*AX1
             INEW=ISH
             JNEW=JSH
             KNEW=KSH
             LNEW=LSH
             IMAX=KTYPE(INEW)-1
             JMAX=KTYPE(JNEW)-1
             KKKMAX=KTYPE(KNEW)-1
             LMAX=KTYPE(LNEW)-1
             PER1=.FALSE.
             PER2=.FALSE.
             PER3=.FALSE.
             IF (IMAX.LT.JMAX) THEN
               INEW = JSH
               JNEW = ISH
               PER1 = .TRUE.
             END IF
             IF (KKKMAX.LT.LMAX) THEN
               KNEW = LSH
               LNEW = KSH
               PER2 = .TRUE.
             END IF
             IF ((IMAX+JMAX).LT.(KKKMAX+LMAX)) THEN
               ID = INEW
               INEW = KNEW
               KNEW = ID
               ID = JNEW
               JNEW = LNEW
               LNEW = ID
               PER3 = .TRUE.
             END IF
             IMAX=3*(KTYPE(INEW)-1)+1
             JMAX=3*(KTYPE(JNEW)-1)+1
             KKKMAX=3*(KTYPE(KNEW)-1)+1
             LMAX=3*(KTYPE(LNEW)-1)+1
             JTYPE=(IMAX+JMAX+KKKMAX+KKKMAX+LMAX-2)/3
             IAT = KATOM(INEW)
             JAT = KATOM(JNEW)
             KAT = KATOM(KNEW)
             LAT = KATOM(LNEW)
             IF ((IAT.EQ.JAT).AND.(IAT.EQ.KAT).AND.(IAT.EQ.LAT))
     1           GO TO 500
           ELSE
C
C     ----- GET -KL- CHARGE DISTRIBUTION -----
C       ACTUALLY JUST THE POINTERS
C
              KKLL=IA(MAX0(KK,LL))+MIN0(KK,LL)
              CALL xOEDRD(X(INIJG),NKL,NKL0,KKLL)
              IF(NKL.EQ.0) GO TO 500
C
C     ----- SELECT CENTERS FOR DERIVATIVES -----
C
              CALL xJKDATM(ISH,JSH,KSH,LSH)
              IF(SKIPI.AND.SKIPJ.AND.SKIPK.AND.SKIPL) GO TO 500
C
C     ----- SET INDICES FOR SHELL BLOCK -----
C
              CALL xJKDSHL(ISH,JSH,KSH,LSH)
              CALL xJKDNDX(X(IIJKLG))
              INEW = ISH
              JNEW = JSH
              KNEW = KSH
              LNEW = LSH
           END IF
C
C     ----- OBTAIN 2 BODY DENSITY FOR THIS SHELL BLOCK -----
C
cwc         IF (NTMO.GT.0) THEN
cwc            CALL DABPAU(INEW,JNEW,KNEW,LNEW,UHFTYP,X(LDEN),
cwc     *                  X(LVEC),X(IDAB),DABMAX,Q4,POPLE)
cwc            GO TO 111
cwc         END IF
cwc         IF(HFSCF .AND. DFTTYP(3).NE.1.0D+00)CALL DABDFT
cwc     *                        (INEW,JNEW,KNEW,LNEW,UHFTYP,X(LDEN),
cwc     *                         X(LVEC),X(IDAB),DABMAX,Q4,POPLE)
cwc         IF(HFSCF .AND. DFTTYP(3).EQ.1.0D+00) CALL xDABCLU
cwc     *                        (INEW,JNEW,KNEW,LNEW,UHFTYP,X(LDEN),
cwc     *                         X(LVEC),X(IDAB),DABMAX,Q4,POPLE)
cwc1
         CALL xDABCLU
     *                        (INEW,JNEW,KNEW,LNEW,UHFTYP,X(LDEN),
     *                         X(LVEC),X(IDAB),DABMAX,Q4,POPLE)
cwc0
cwc         IF(ROGVB) CALL DABGVB(INEW,JNEW,KNEW,LNEW,NOCORE,NOOPEN,
cwc     *                         X(LDEN),X(LVEC),NUM,X(IDAB),
cwc     *                         DABMAX,Q4,POPLE)
cwc         IF(MC) THEN
cwc            IF(NCI.GT.0) CALL DABMCL(X(K30),X(K70),X(K80),NCI,
cwc     *                               ISH,JSH,KSH,LSH,NUM,MAXSHL,NMCCOR)
cwc            CALL DABMC(INEW,JNEW,KNEW,LNEW,X(K10),X(K20),X(K80),NMCCOR,
cwc     *                 NCI,MAXSHL,X(IDAB),DABMAX,Q4,POPLE,
cwc     *                 PER1,PER2,PER3)
cwc         END IF
cwc         IF(MP2) THEN
cwc            CALL DABMP2(INEW,JNEW,KNEW,LNEW,X(LMPDM2),X(LVEC),X(LPHF),
cwc     *                  X(LPMP2),X(IDAB),DABMAX,L1,L2,Q4,POPLE)
cwc         END IF
cwc         IF(UMP2) THEN
cwc            CALL DABUMP(INEW,JNEW,KNEW,LNEW,X(LMPDM2),X(LMPDM2B),
cwc     *                  X(LVECA),X(LVECB),X(LPHFA),X(LPHFB),X(LPMP2A),
cwc     *                  X(LPMP2B),X(IDAB),DABMAX,L1,L2,Q4,POPLE)
cwc         END IF
cwc         IF(LCIS) THEN
cwc            CALL DABCIS(INEW,JNEW,KNEW,LNEW,X(LPHF),X(LPCIS),
cwc     *                  X(LTCIS),X(IDAB),DABMAX,L1,L2,Q4,POPLE)
cwc         END IF
  111    CONTINUE
C
C                FINE SCREENING, ON INTEGRAL VALUE TIMES DENSITY FACTOR
C
         IF(DABMAX*GMAX.LT.CUTOFF2) THEN
            IISKIP2 = IISKIP2+1
            GO TO 500
         END IF
C
C     ----- EVALUATE DERIVATIVE INTEGRAL, AND ADD TO THE GRADIENT -----
C
         IDID = IDID+1
         IF(POPLE) THEN
            CALL xJKDG80(DABMAX,INEW,JNEW,KNEW,LNEW,
     *                  JTYPE,IAT,JAT,KAT,LAT)
         ELSE
c            write (96,'(a,4i3)') 'ish,jsh,ksh,lsh:',ish,jsh,ksh,lsh
            CALL xJKDSPD(NIJ0,NKL,NKL0,X(ICHRG),
     *      X(IGINT),X(IFINT),X(ISINT),X(IIJKLG),X(IGIJKL),
     1      X(IGNKL),X(IGNM),X(IXY),X(IXZ),X(IYZ),X(IX),X(IY),X(IZ),
     2      X(ISJ),X(ISK),X(ISL),X(IB00),X(IB01),X(IB10),X(IC00),
     3      X(ID00),X(IF00),X(IDIJ),X(IDKL),X(IDIJSI),X(IDIJSJ),
     4      X(IDKLSK),X(IDKLSL),X(IABV),X(ICV),X(IRW),X(IAAI),X(IAAJ),
     5      X(IBBK),X(IBBL),X(IFI),X(IFJ),X(IFK),X(IFL),X(ISII),X(ISJJ),
     6      X(ISKK),X(ISLL),X(ISIJ),X(ISIK),X(ISIL),X(ISJK),X(ISJL),
     7      X(ISKL),X(IDAB),MAXXYZ,FC,NC,DF,LDF,NBF,DDA,Q4,MINVEC,
     8      DABCUT,DABMAX)
         END IF
C
C     ----- END OF *SHELL* LOOPS -----
C
  500 CONTINUE
  520 CONTINUE
  540 CONTINUE
cwc      IF (TIM .GE. TIMLIM) GO TO 600
  560 CONTINUE
C
C     ----- GO BACK FOR ANOTHER MP2 PASS, IF NECESSARY -----
C
cwc      IF(MP2) THEN
cwc         NOC1 = NOC2+1
cwc         IF(NOC1.LE.NOA) GO TO 350
cwc      END IF
cwc      IF(UMP2) THEN
cwc         NOC1 = NOC2+1
cwc         NOC3 = NOC4+1
cwc         IF(NOC1.LE.NOCA) GO TO 360
cwc         IF(NOC3.LE.NOCB) GO TO 360
cwc      END IF
C
C     ----- FINISH UP THE FINAL GRADIENT -----
C
C           GLOBAL SUM CONTRIBUTIONS FROM EACH NODE
C
      IF (GOPARR) THEN
cwc         IF (NXT) CALL DDI_DLBRESET
cwc         CALL DDI_GSUMF(1600,DE,3*NAT)
cwc         CALL DDI_GSUMI(1601,IISKIP1,1)
cwc         CALL DDI_GSUMI(1601,IISKIP2,1)
cwc         CALL DDI_GSUMI(1602,IDID,1)
      END IF
cwc1
      iskp1i=IISKIP1
      iskp2i=IISKIP2
      ididi=IDID
c      call timer(0)
      call MPIAR (DE,buf,3*nat)
      call MPIAR1i (iskp1i,IISKIP1)
      call MPIAR1i (iskp2i,IISKIP2)
      call MPIAR1i (ididi,IDID)
c      call timer(1)
c      write (96,'(a,f10.2)') '2-el MPIAR time:     ',timenew
cwc0
C
C           ADD IN THE REST OF THE DFT GRADIENT
C
cwc      IF(DFTTYP(1) .NE. 0.0D+00) THEN
cwc         CALL DFTGRD(X(IDFT1),X(IDFT2),X(IDFT3),X(IDFT4),X(IDFT5),
cwc     *               X(IDFT6),X(IDFT7),NAUXFUN)
cwc         CALL VADD(X(IDFT7),1,DE,1,DE,1,3*NAT)
cwc      END IF
C
C           SYMMETRIZE SKELETON GRADIENT VECTOR
C
cwc      CALL SYMEG(DE)
C
C           PROJECT ROTATIONAL CONTAMINANT FROM DFT GRADIENTS
C
cwc      IF(NDFTFG.EQ.1  .OR.  DFTTYP(1).NE.0.0D+00) THEN
cwc         NCCF = 3*NAT + 6*NFRG
cwc         CALL VALFM(LOADFM)
cwc         LX   = LOADFM+ 1
cwc         LGTOT= LX    + NCCF
cwc         LRM  = LGTOT + NCCF
cwc         LP   = LRM   + 3*NAT + 21*NFRG
cwc         LAST = LP    + NCCF*NCCF
cwc         NEEDP= LAST-LOADFM-1
cwc         CALL GETFM(NEEDP)
cwc         CALL PRJGRD(DE,X(LX),X(LGTOT),X(LRM),X(LP),3*NAT,NCCF)
cwc         CALL RETFM(NEEDP)
cwc      END IF
C
cwc1
c      write (96,*)
c      write (96,'(a)') 'FINAL GRADIENT'
c      do 7777 i=1,nat
c 7777 write (96,'(i5,f29.9,2f18.9)') i,de(1,i),de(2,i),de(3,i)
cwc0
      CALL xDFINAL(1)
C
C     ----- DEALLOCATE MEMORY -----
C
  600 CONTINUE
      IF(LAST.GT.0) LAST=0
      CALL yRETFM(NEED)
cwc1
      nprtgo=0
cwc0

cwc      IF (MASWRK .AND. NPRTGO.NE.2) THEN
cwc        IF(SOME) WRITE(IW,9999) IISKIP1,IISKIP2,IDID
cwc        WRITE(IW,FMT='(/'' ...... END OF 2-ELECTRON GRADIENT ......'')')
cwc      END IF
cwc      CALL TEXIT(1,4)
cwc      IF(MP2) NOA=NOASAV
cwc1
        IF (irank.eq.0) WRITE(IW,9999) IISKIP1,IISKIP2,IDID
cwc0
      RETURN
C
 9008 FORMAT(/10X,22(1H-)/10X,22HGRADIENT OF THE ENERGY/10X,22(1H-))
 9995 FORMAT(1X,'UMP2 MEMORY ALLOCATION FOR 2E- GRADIENT INTEGRALS.'/
     *       1X,'THERE ARE',I5,' ALPHA OCCUPIED MO-S AND',I5,
     *          ' BETA OCCUPIED MO-S'/
     *       1X,'MINIMUM MEMORY =',I15,
     *          ' (ONE OCCUPIED MO/DERIVATIVE INTEGRAL PASS)'/
     *       1X,'MAXIMUM MEMORY =',I15,
     *          ' (ALL OCCUPIED MO-S IN A SINGLE PASS)'/
     *       1X,'   MEMORY USED =',I15,', WITH ',I5,
     *          ' DERIVATIVE PASSES,',I5,' MO-S/PASS')
 9996 FORMAT(1X,'INCREASE MEMORY ',I15,
     *          ' WORDS TO REDUCE PASS COUNT BY ONE.')
 9997 FORMAT(/' NOT ENOUGH MEMORY FOR -2DM- CONTRIBUTION'/
     *        ' INCREASE MEMORY SIZE BY AT LEAST ',I10,' WORDS')
 9998 FORMAT(1X,'MP2 MEMORY ALLOCATION FOR 2E- GRADIENT INTEGRALS.  (',
     *          I5,' OCCUPIED MO-S)'/
     *       1X,'MINIMUM MEMORY =',I15,
     *          ' (ONE OCCUPIED MO/DERIVATIVE INTEGRAL PASS)'/
     *       1X,'MAXIMUM MEMORY =',I15,
     *          ' (ALL OCCUPIED MO-S IN A SINGLE PASS)'/
     *       1X,'   MEMORY USED =',I15,', WITH ',I5,
     *          ' DERIVATIVE PASSES,',I5,' MO-S/PASS')
 9999 FORMAT(1X,'THE COARSE/FINE SCHWARZ SCREENINGS SKIPPED ',I10,'/',
     *          I10,' BLOCKS.'/
     *     1X,'THE NUMBER OF GRADIENT INTEGRAL BLOCKS COMPUTED WAS',I10)
      END
c------------------------------------------------------------------------
C*MODULE GRD2A   *DECK OEDSHL
      SUBROUTINE xOEDSHL(ISH,JSH,DCHRG,NIJ0)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "gamess.par"
C
      DIMENSION GA(MXGSH),CCA(MXGSH),CCAS(MXGSH),
     *          GB(MXGSH),CCB(MXGSH),CCBS(MXGSH)
      DIMENSION CSPDFG(MXGTOT,5)
      DIMENSION DCHRG(15,*)
C
      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL, use_symm
      LOGICAL SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      LOGICAL EXPNDI,EXPNDK,OUT,DBG
C
      COMMON /xDSHLNO/ LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     *                NIJ,IJD,KLD,IJ,KL
      COMMON /xDSHLT / RTOL,DTOL,VTOL1,VTOL2,VTOLS,OUT,DBG
      COMMON /xINFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM)
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON/xSHLEQU/IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL,use_symm
      COMMON/xSHLNUM/NUMI,NUMJ,NUMK,NUML,IJKL
      COMMON/xSHLTYP/SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      COMMON/xSHLXPN/EXPNDI,EXPNDK
C
      EQUIVALENCE (CSPDFG(1,1),CS(1))
C
      PARAMETER (ONE=1.0D+00)
C
      if (use_symm) then
         IIEQJJ=ISH.EQ.JSH
      else
         IIEQJJ=.FALSE. !VFL 
      endif
C
C     ----- ISHELL -----
C
      I=KATOM(ISH)
      XI=C(1,I)
      YI=C(2,I)
      ZI=C(3,I)
      I1=KSTART(ISH)
      I2=I1+KNG(ISH)-1
      LIT=KTYPE(ISH)
      MINI=KMIN(ISH)
      MAXI=KMAX(ISH)
      NUMI=MAXI-MINI+1
      LOCI=KLOC(ISH)-MINI
      SPI=LIT.EQ.2.AND.MINI.EQ.1
      NGA=0
      DO 10 I=I1,I2
        NGA=NGA+1
        GA(NGA)=EX(I)
        CCA(NGA)=CSPDFG(I,LIT)
        IF(SPI) CCAS(NGA)=CSPDFG(I,1)/CSPDFG(I,2)
   10 CONTINUE
C
C     ----- JSHELL -----
C
      J=KATOM(JSH)
      XJ=C(1,J)
      YJ=C(2,J)
      ZJ=C(3,J)
      J1=KSTART(JSH)
      J2=J1+KNG(JSH)-1
      LJT=KTYPE(JSH)
      MINJ=KMIN(JSH)
      MAXJ=KMAX(JSH)
      NUMJ=MAXJ-MINJ+1
      LOCJ=KLOC(JSH)-MINJ
      SPJ=LJT.EQ.2.AND.MINJ.EQ.1
      NGB=0
      DO 20 J=J1,J2
        NGB=NGB+1
        GB(NGB)=EX(J)
        CCB(NGB)=CSPDFG(J,LJT)
        IF(SPJ) CCBS(NGB)=CSPDFG(J,1)/CSPDFG(J,2)
   20 CONTINUE
      RRI=((XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2)
      SPIJ=SPI.OR.SPJ
      EXPNDI=LIT.GE.LJT
C
C     ----- -IJ- CHARGE DISTRIBUTION -----
C
      XC=XI
      YC=YI
      ZC=ZI
      DXIJ=XI-XJ
      DYIJ=YI-YJ
      DZIJ=ZI-ZJ
      IF(EXPNDI) GO TO 100
      XC=XJ
      YC=YJ
      ZC=ZJ
      DXIJ=XJ-XI
      DYIJ=YJ-YI
      DZIJ=ZJ-ZI
  100 CONTINUE
C
C     ----- - I- PRIMITIVE           -----
C
      NIJ=0
      DO 300 IA=1,NGA
        AI=GA(IA)
        ARRI=AI*RRI
        AXI=AI*XI
        AYI=AI*YI
        AZI=AI*ZI
        CCI=CCA(IA)
C
C     ----- - J- PRIMITIVE           -----
C
        DO 200 JB=1,NGB
          AJ=GB(JB)
          AA=AI+AJ
          AA1=ONE/AA
          DUM=AJ*ARRI*AA1
          IF(DUM.GT.RTOL) GO TO 200
          DAEXPA=CCI*CCB(JB)* EXP(-DUM)*AA1
          DUM=  ABS(DAEXPA)
          IF(DUM.LE.DTOL) GO TO 200
C
          NIJ=NIJ+1
          DCHRG( 1,NIJ+NIJ0)= DAEXPA
          DCHRG( 2,NIJ+NIJ0)= AA
          DCHRG( 3,NIJ+NIJ0)=(AXI+AJ*XJ)*AA1
          DCHRG( 4,NIJ+NIJ0)=(AYI+AJ*YJ)*AA1
          DCHRG( 5,NIJ+NIJ0)=(AZI+AJ*ZJ)*AA1
          DCHRG( 6,NIJ+NIJ0)= XC
          DCHRG( 7,NIJ+NIJ0)= YC
          DCHRG( 8,NIJ+NIJ0)= ZC
          DCHRG( 9,NIJ+NIJ0)= DXIJ
          DCHRG(10,NIJ+NIJ0)= DYIJ
          DCHRG(11,NIJ+NIJ0)= DZIJ
          DCHRG(12,NIJ+NIJ0)= AI+AI
          DCHRG(13,NIJ+NIJ0)= AJ+AJ
          IF(SPI) DCHRG(14,NIJ+NIJ0)=CCAS(IA)
          IF(SPJ) DCHRG(15,NIJ+NIJ0)=CCBS(JB)
C
  200   CONTINUE
  300 CONTINUE
      RETURN
      END
c------------------------------------------------------------------------
C*MODULE GRD2A   *DECK OEDRD
      SUBROUTINE xOEDRD(NIJG,NIJ,NIJ0,IIJJ)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION NIJG(2,*)
C
      NIJ0=NIJG(1,IIJJ)
      NIJ =NIJG(2,IIJJ)
      RETURN
      END
c--------------------------------------------------------------------------
C*MODULE GRD2A   *DECK JKDATM
      SUBROUTINE xJKDATM(II,JJ,KK,LL)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "gamess.par"
C
      LOGICAL DBUG
      LOGICAL OUT
      LOGICAL SKIPI,SKIPJ,SKIPK,SKIPL
      LOGICAL IANDJ,IANDK,IANDL,JANDK,JANDL,KANDL
      COMMON /xDSHLT / RTOL,DTOL,VTOL1,VTOL2,VTOLS,OUT,DBUG
      COMMON /xIOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON/xDERSKP/IIAT,JJAT,KKAT,LLAT,SKIPI,SKIPJ,SKIPK,SKIPL
      COMMON/xDERINV/INVTYP
cwc1
      common /gms_deriv_ctrl/ want_integrals, all_centers
      logical want_integrals, all_centers
cwc0
C
      SKIPI=.TRUE.
      SKIPJ=.TRUE.
      SKIPK=.TRUE.
      SKIPL=.TRUE.
      IIAT=KATOM(II)
      JJAT=KATOM(JJ)
      KKAT=KATOM(KK)
      LLAT=KATOM(LL)
      IANDJ=IIAT.EQ.JJAT
      IANDK=IIAT.EQ.KKAT
      IANDL=IIAT.EQ.LLAT
      JANDK=JJAT.EQ.KKAT
      JANDL=JJAT.EQ.LLAT
      KANDL=KKAT.EQ.LLAT
C
      IF(.NOT.IANDJ) GO TO 500
      IF(.NOT.IANDK) GO TO 200
      IF(.NOT.IANDL) GO TO 100
C
C     ----- IAT = JAT = KAT = LAT ----- (IAT,IAT/IAT,IAT) -----
C
      INVTYP=1
      GO TO 1500
  100 CONTINUE
C
C     ----- IAT = JAT = KAT ; LAT ----- (IAT,IAT/IAT,LAT) -----
C
      SKIPL=.FALSE.
      INVTYP=2
      GO TO 1500
  200 IF(.NOT.IANDL) GO TO 300
C
C     ----- IAT = JAT = LAT ; KAT ----- (IAT,IAT/KAT,IAT) -----
C
      SKIPK=.FALSE.
      INVTYP=3
      GO TO 1500
  300 IF(.NOT.KANDL) GO TO 400
C
C     ----- IAT = JAT ; KAT = LAT ----- (IAT,IAT/KAT,KAT) -----
C
      SKIPK=.FALSE.
      SKIPL=.FALSE.
      INVTYP=4
      GO TO 1500
  400 CONTINUE
C
C     ----- IAT = JAT ; KAT ; LAT ----- (IAT,IAT/KAT,LAT) -----
C
      SKIPK=.FALSE.
      SKIPL=.FALSE.
      INVTYP=5
      GO TO 1500
  500 IF(.NOT.IANDK) GO TO 800
      IF(.NOT.IANDL) GO TO 600
C
C     ----- IAT = KAT = LAT ; JAT ----- (IAT,JAT/IAT,IAT) -----
C
      SKIPJ=.FALSE.
      INVTYP=6
      GO TO 1500
  600 IF(.NOT.JANDL) GO TO 700
C
C     ----- IAT = KAT ; JAT = LAT ----- (IAT,JAT/IAT,JAT) -----
C
      SKIPJ=.FALSE.
      SKIPL=.FALSE.
      INVTYP=7
      GO TO 1500
  700 CONTINUE
C
C     ----- IAT = KAT ; JAT ; LAT ----- (IAT,JAT/IAT,LAT) -----
C
      SKIPJ=.FALSE.
      SKIPL=.FALSE.
      INVTYP=8
      GO TO 1500
  800 IF(.NOT.IANDL) GO TO 1000
      IF(.NOT.JANDK) GO TO 900
C
C     ----- IAT = LAT ; JAT = KAT ----- (IAT,JAT/JAT,IAT) -----
C
      SKIPJ=.FALSE.
      SKIPK=.FALSE.
      INVTYP=9
      GO TO 1500
  900 CONTINUE
C
C     ----- IAT = LAT ; JAT , KAT ----- (IAT,JAT/KAT,IAT) -----
C
      SKIPJ=.FALSE.
      SKIPK=.FALSE.
      INVTYP=10
      GO TO 1500
 1000 IF(.NOT.JANDK) GO TO 1200
      IF(.NOT.JANDL) GO TO 1100
C
C     ----- IAT ; JAT = JAT = JAT ----- (IAT,JAT/JAT,JAT) -----
C
      SKIPI=.FALSE.
      INVTYP=11
      GO TO 1500
 1100 CONTINUE
C
C     ----- IAT ; JAT = KAT ; LAT ----- (IAT,JAT/JAT,LAT) -----
C
      SKIPI=.FALSE.
      SKIPL=.FALSE.
      INVTYP=12
      GO TO 1500
 1200 IF(.NOT.JANDL) GO TO 1300
C
C     ----- JAT = LAT ; IAT ; KAT ----- (IAT,JAT/KAT,JAT) -----
C
      SKIPI=.FALSE.
      SKIPK=.FALSE.
      INVTYP=13
      GO TO 1500
 1300 IF(.NOT.KANDL) GO TO 1400
C
C     ----- KAT = LAT ; IAT ; JAT ----- (IAT,JAT/KAT,KAT) -----
C
      SKIPI=.FALSE.
      SKIPJ=.FALSE.
      INVTYP=14
      GO TO 1500
 1400 CONTINUE
C
C     ----- IAT ; JAT ; KAT ; LAT ----- (IAT,JAT/KAT,LAT) -----
C
      SKIPI=.FALSE.
      SKIPJ=.FALSE.
      SKIPK=.FALSE.
      INVTYP=15
C
 1500 CONTINUE
C
C     ----- FOR DEBUGGING PURPOSES CALCULATE ALL TERMS -----
C
cwc      IF(DBUG) THEN
      IF(DBUG.or.all_centers) THEN
        SKIPI=.FALSE.
        SKIPJ=.FALSE.
        SKIPK=.FALSE.
        SKIPL=.FALSE.
        INVTYP=16
      END IF
C
      IF(OUT) WRITE(IW,9999) II,JJ,KK,LL,SKIPI,SKIPJ,SKIPK,SKIPL
      RETURN
 9999 FORMAT(/,' -- II,JJ,KK,LL =',4I3,' SKIPI,J,K,L =',4L3)
      END
c--------------------------------------------------------------------------
C*MODULE GRD2A   *DECK DABCLU
      SUBROUTINE xDABCLU(II,JJ,KK,LL,UHFTYP,DA,DB,DAB,DABMAX,Q4,POPLE)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION DA(*),DB(*),DAB(*)
      LOGICAL UHFTYP,POPLE
C
#include "gamess.par"
C
      LOGICAL OUT,DBG
      COMMON /xDSHLT / RTOL,DTOL,VTOL1,VTOL2,VTOLS,OUT,DBG
      COMMON /xIJPAIR/ IA(MXAO)
      COMMON /xIOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /xINDD80/ IMAX,JMAX,KKMAX,LMAX
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /xSHLLMN/ IGXYZ(4,35),JGXYZ(4,35),KGXYZ(4,35),LGXYZ(4,35)
      COMMON /xSHLNRM/ PNRM(35)
      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL, use_symm
      COMMON /xSHLEQU/ IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL, use_symm
C
CJMS  LABELLED COMMON GSPG80 DEFINED FOR COMPUTATIONAL EFFICIENCY.
CJMS  FOR SP BASES ONLY, IT CONTAINS THE E ARRAY WHICH IS THE DAB
CJMS  ARRAY WITH INDICES IN REVERSE ORDER: E(I,J,K,L)= DAB(L,K,J,I)
CJMS  AND IS USED IN SUB JKDG80 (MODULE GRD2B). IT ORIGINATES IN:
CJMS
CJMS     1. SUBS DABCLU, DABDFT, DABGVB, DABMC, DABMP2 AND DABUMP
CJMS        (MODULE GRD2A) AND SUB DABPAU (MODULE EFPAUL) WHICH ARE
CJMS        ALL CALLED BY SUB JKDER (MODULE GRD2A)
CJMS
CJMS     2. SUB DABCLU (MODULE GRD2A) WHICH IS CALLED BY SUB EFDEN OF
CJMS        MODULE EFGRD2
CJMS
CJMS     3. SUB PAR2PDM (MODULE MP2DDI) WHICH IS CALLED BY SUB PJKDMP2
CJMS        OF MODULE MP2DDI
C
      COMMON /xGSPG80/ E(4,4,4,4)
C
      PARAMETER (ZER=0.0D+00, PT5=0.5D+00, P25=.25D+00, F04=4.0D+00)
C
C     ----- THIS ROUTINE FORMS THE PRODUCT OF DENSITY       -----
C           MATRICES FOR USE IN FORMING THE TWO ELECTRON
C           GRADIENT.  VALID FOR CLOSED AND OPEN SHELL SCF.
C
      DABMAX= ZER
      MINI= KMIN(II)
      MINJ= KMIN(JJ)
      MINK= KMIN(KK)
      MINL= KMIN(LL)
      LOCI= KLOC(II)-MINI
      LOCJ= KLOC(JJ)-MINJ
      LOCK= KLOC(KK)-MINK
      LOCL= KLOC(LL)-MINL
      IF(POPLE) THEN
         DO 110 L=1,LMAX
            DO 110 K=1,KKMAX
               KAL= MAX0(LOCK+K,LOCL+L)
               KIL= MIN0(LOCK+K,LOCL+L)
               DO 110 J=1,JMAX
                  DO 110 I=1,IMAX
                     IAJ= MAX0(LOCI+I,LOCJ+J)
                     IIJ= MIN0(LOCI+I,LOCJ+J)
                     IN = IAJ
                     JN = IIJ
                     KN = KAL
                     LN = KIL
                     IF(IN.LT.KN .OR.(IN.EQ.KN .AND. JN.LT.LN)) THEN
                        IN = KAL
                        JN = KIL
                        KN = IAJ
                        LN = IIJ
                     ENDIF
C        -----  NOW CALCULATE DENSITY MATRIX CONTRIBUTION.      -----
C        -----  (IJ/KL),(IJ/IL),(IJ/JL),(IJ/KJ)                 -----
C        -----  EIGHT DISTINCT INTEGRALS                        -----
C        -----  CONTRIBUTION TO THE ENERGY                      -----
C        -----  4(IJ)(KL) - (IK)(JL) - (JK)(IL)                 -----
                     IJ = IA(IN)+JN
                     IK = IA(IN)+KN
                     IL = IA(IN)+LN
                     JK = IA(MAX0(JN,KN))+MIN0(JN,KN)
                     JL = IA(JN)+LN
                     IF(JN.LT.KN) JL = IA(MAX0(JN,LN))+MIN0(JN,LN)
                     KL = IA(KN)+LN
                     DF1= DA(IJ)*DA(KL)
                     DQ1= DA(IK)*DA(JL)+DA(IL)*DA(JK)
                     IF(UHFTYP) DQ1=DQ1+DB(IK)*DB(JL)+DB(IL)*DB(JK)
                     DF1=(DF1-P25*DQ1)*Q4
                     IF(DABMAX.LT. ABS(DF1)) DABMAX= ABS(DF1)
                     E(I,J,K,L)= DF1
                     IF(OUT) WRITE(IW,9010) II,JJ,KK,LL,I,J,K,L,DF1
  110    CONTINUE
      ELSE
C
C D AND HIGHER FUNCTIONS OR HONDO ONLY RUNS
C
         MAXI= KMAX(II)
         MAXJ= KMAX(JJ)
         MAXK= KMAX(KK)
         MAXL= KMAX(LL)
         DO 210 I=MINI,MAXI
            P1I = PNRM(I)
            JMAX= MAXJ
            IF(use_symm .and. IIEQJJ) JMAX= I
            DO 210 J=MINJ,JMAX
               P2J = P1I*PNRM(J)
               IAJ= MAX0(LOCI+I,LOCJ+J)
               IIJ= MIN0(LOCI+I,LOCJ+J)
               KMMAX=MAXK
               IF(use_symm .and. IJEQKL) KMMAX= I
               DO 210 K=MINK,KMMAX
                  P3K = P2J*PNRM(K)
                  LMAX= MAXL
                  if (use_symm) then
                     IF(KKEQLL) LMAX= K
                     IF(IJEQKL .AND. K.EQ.I) LMAX= J
                  endif

                  DO 210 L=MINL,LMAX
                     P4L= P3K*PNRM(L)
                     KAL= MAX0(LOCK+K,LOCL+L)
                     KIL= MIN0(LOCK+K,LOCL+L)
                     IN = IAJ
                     JN = IIJ
                     KN = KAL
                     LN = KIL
                     IF(IN.LT.KN .OR.(IN.EQ.KN .AND. JN.LT.LN)) THEN
                        IN = KAL
                        JN = KIL
                        KN = IAJ
                        LN = IIJ
                     ENDIF
                     IJ = IA(IN)+JN
                     IK = IA(IN)+KN
                     IL = IA(IN)+LN
                     JK = IA(MAX0(JN,KN))+MIN0(JN,KN)
                     JL = IA(JN)+LN
                     IF(JN.LT.KN) JL = IA(MAX0(JN,LN))+MIN0(JN,LN)
                     KL = IA(KN)+LN
                     DF1= DA(IJ)*DA(KL)
                     DQ1= DA(IK)*DA(JL)+DA(IL)*DA(JK)
                     IF(UHFTYP) DQ1=DQ1+DB(IK)*DB(JL)+DB(IL)*DB(JK)
                     DF1= DF1*F04-DQ1
                     IF(JN.EQ.IN               ) DF1= DF1*PT5
                     IF(LN.EQ.KN               ) DF1= DF1*PT5
                     IF(KN.EQ.IN .AND. LN.EQ.JN) DF1= DF1*PT5
                     IF(DABMAX.LT. ABS(DF1)) DABMAX= ABS(DF1)
C
C IGXYZ AND J, K, AND L ARE SET UP IN JKDNDX
C
                     IJKL=IGXYZ(1,I)+JGXYZ(1,J)+KGXYZ(1,K)+LGXYZ(1,L)
                     DAB(IJKL)= DF1*P4L
                     IF(OUT) WRITE(IW,9020) II,JJ,KK,LL,I,J,K,L,IJKL,DF1
  210    CONTINUE
      ENDIF
      RETURN
 9010 FORMAT(' -DABCLU,POPLE- ',4I4,4I3,D20.12)
 9020 FORMAT(' -DABCLU,HONDO- ',4I4,4I3,I5,D20.12)
      END
c--------------------------------------------------------------------------
C*MODULE GRD2A   *DECK JKDNDX
      SUBROUTINE xJKDNDX(IJKLG)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL SKIPI,SKIPJ,SKIPK,SKIPL
      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL, use_symm
C--   LOGICAL IEQJ0,KEQL0,SAME0,FIRST
      COMMON /xDSHLNO/ LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     *                NIJ,IJD,KLD,IJ,KL
      COMMON /xROOT  / XX,U(13),W(13),NROOTS
      COMMON/xSHLNUM/NUMI,NUMJ,NUMK,NUML,IJKL
      COMMON/xSHLLMN/IGXYZ(4,35),JGXYZ(4,35),KGXYZ(4,35),LGXYZ(4,35)
      COMMON/xSHLEQU/IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL, use_symm
C--   COMMON/xSHLOLD/NUMI0,NUMJ0,NUMK0,NUML0,IEQJ0,KEQL0,SAME0,FIRST
      COMMON/xDERPAR/IDER,JDER,KDER,LDER,NDER
      COMMON/xDERSKP/IIAT,JJAT,KKAT,LLAT,SKIPI,SKIPJ,SKIPK,SKIPL
      DIMENSION IJKLX(35),IJKLY(35),IJKLZ(35)
      DIMENSION IJKLN(5)
      DIMENSION IJKLG(4,*)
C
      DATA IJKLN /   1,  4, 10, 20, 35/
      DATA IJKLX /   0,  1,  0,  0,  2,  0,  0,  1,  1,  0,
     1               3,  0,  0,  2,  2,  1,  0,  1,  0,  1,
     2               4,  0,  0,  3,  3,  1,  0,  1,  0,  2,
     3               2,  0,  2,  1,  1/
      DATA IJKLY /   0,  0,  1,  0,  0,  2,  0,  1,  0,  1,
     1               0,  3,  0,  1,  0,  2,  2,  0,  1,  1,
     2               0,  4,  0,  1,  0,  3,  3,  0,  1,  2,
     3               0,  2,  1,  2,  1/
      DATA IJKLZ /   0,  0,  0,  1,  0,  0,  2,  0,  1,  1,
     1               0,  0,  3,  0,  1,  0,  1,  2,  2,  1,
     2               0,  0,  4,  0,  1,  0,  1,  3,  3,  0,
     3               2,  2,  1,  1,  2/
C
C...  SAMTYP=NUMI.EQ.NUMI0.AND.NUMJ.EQ.NUMJ0.AND.
C... 1       NUMK.EQ.NUMK0.AND.NUML.EQ.NUML0          .AND.(.NOT.FIRST)
C...  SAMSYM=(IIEQJJ.EQV.IEQJ0).AND.
C... 1       (KKEQLL.EQV.KEQL0).AND.(IJEQKL.EQV.SAME0).AND.(.NOT.FIRST)
C--   NUMI0=NUMI
C--   NUMJ0=NUMJ
C--   NUMK0=NUMK
C--   NUML0=NUML
C--   IEQJ0=IIEQJJ
C--   KEQL0=KKEQLL
C--   SAME0=IJEQKL
C--   FIRST=.FALSE.
C.... IF(SAMTYP.AND.SAMSYM) RETURN
C
      IDER=NDER
      JDER=NDER
      KDER=NDER
      LDER=NDER
      IF(SKIPI) IDER=0
      IF(SKIPJ) JDER=0
      IF(SKIPK) KDER=0
      IF(SKIPL) LDER=0
      LJTMOD=LJT + JDER
      LKTMOD=LKT + KDER
      LLTMOD=LLT + LDER
C
C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS -----
C
      NI=NUML*NUMK*NUMJ
      DO 10 I=MINI,MAXI
   10 IGXYZ(1,I)=NI*(I-MINI)+1
      LLKJT=LLTMOD*LKTMOD*LJTMOD
      DO 20 I=1,IJKLN(LIT)
        IGXYZ(2,I)=IJKLX(I)*LLKJT+1
        IGXYZ(3,I)=IJKLY(I)*LLKJT+1
   20   IGXYZ(4,I)=IJKLZ(I)*LLKJT+1
C
      NJ=NUML*NUMK
      DO 30 J=MINJ,MAXJ
   30 JGXYZ(1,J)=NJ*(J-MINJ)
      LLKT=LLTMOD*LKTMOD
      DO 40 J=1,IJKLN(LJT)
        JGXYZ(2,J)=IJKLX(J)*LLKT
        JGXYZ(3,J)=IJKLY(J)*LLKT
   40   JGXYZ(4,J)=IJKLZ(J)*LLKT
C
C     ----- PREPARE INDICES FOR PAIRS OF (K,L) FUNCTIONS -----
C
      NK=NUML
      DO 110 K=MINK,MAXK
  110 KGXYZ(1,K)=NK*(K-MINK)
      DO 120 K=1,IJKLN(LKT)
        KGXYZ(2,K)=IJKLX(K)*LLTMOD
        KGXYZ(3,K)=IJKLY(K)*LLTMOD
  120   KGXYZ(4,K)=IJKLZ(K)*LLTMOD
C
      NL=1
      DO L=MINL,MAXL
         LGXYZ(1,L)=NL*(L-MINL)
      ENDDO
      DO L=1,IJKLN(LLT)
         LGXYZ(2,L)=IJKLX(L)
         LGXYZ(3,L)=IJKLY(L)
         LGXYZ(4,L)=IJKLZ(L)
      ENDDO
C
C     ----- PREPARE INDICES FOR (IJ/KL) -----
C
      IJKL=0
      DO 240 I=MINI,MAXI
        JMAX=MAXJ
        IF(use_symm .and. IIEQJJ) JMAX=I
        DO 230 J=MINJ,JMAX
          KMAX=MAXK
          IF(use_symm .and. IJEQKL) KMAX=I
          DO 220 K=MINK,KMAX
            LMAX=MAXL
            if (use_symm) then
               IF(KKEQLL           ) LMAX=K
               IF(IJEQKL.AND.K.EQ.I) LMAX=J
            endif

            DO 210 L=MINL,LMAX
              IJKL=IJKL+1
              NN=((IGXYZ(1,I)+JGXYZ(1,J))+KGXYZ(1,K))+LGXYZ(1,L)
              NX=((IGXYZ(2,I)+JGXYZ(2,J))+KGXYZ(2,K))+LGXYZ(2,L)
              NY=((IGXYZ(3,I)+JGXYZ(3,J))+KGXYZ(3,K))+LGXYZ(3,L)
              NZ=((IGXYZ(4,I)+JGXYZ(4,J))+KGXYZ(4,K))+LGXYZ(4,L)
              IJKLG(1,IJKL)=   NN
              IJKLG(2,IJKL)=3*(NX-1)+1
              IJKLG(3,IJKL)=3*(NY-1)+2
              IJKLG(4,IJKL)=3*(NZ-1)+3
  210       CONTINUE
  220     CONTINUE
  230   CONTINUE
  240 CONTINUE
C
C     ----- SET NUMBER OF QUADRATURE POINTS -----
C
      NROOTS=(LIT+LJT+LKT+LLT-2 + NDER )/2
C
      RETURN
      END
c--------------------------------------------------------------------------
C*MODULE GRD2A   *DECK JKDSHL
      SUBROUTINE xJKDSHL(ISH,JSH,KSH,LSH)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "gamess.par"
C
      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL, use_symm
      LOGICAL SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      LOGICAL EXPNDI,EXPNDK
      COMMON /xDSHLNO/ LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     *                NIJ,IJD,KLD,IJ,KL
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON/xSHLTYP/SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      COMMON/xSHLXPN/EXPNDI,EXPNDK
      COMMON/xSHLNUM/NUMI,NUMJ,NUMK,NUML,IJKL
      COMMON/xSHLEQU/IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL, use_symm
C
      if (use_symm) then
         IIEQJJ=ISH.EQ.JSH
         KKEQLL=KSH.EQ.LSH
         IJEQKL=ISH.EQ.KSH.AND.JSH.EQ.LSH
      else
         IIEQJJ=.FALSE. ! VFL
         KKEQLL=.FALSE. ! VFL 
         IJEQKL=.FALSE. ! VFL 
      endif
      IJGTKL=MAX0(ISH,JSH).GT.MAX0(KSH,LSH)
      IJLTKL=MAX0(ISH,JSH).LT.MAX0(KSH,LSH)
C
C     ----- ISHELL -----
C
      LIT=KTYPE(ISH)
      MINI=KMIN(ISH)
      MAXI=KMAX(ISH)
      NUMI=MAXI-MINI+1
      LOCI=KLOC(ISH)-MINI
      SPI=LIT.EQ.2.AND.MINI.EQ.1
C
C     ----- JSHELL -----
C
      LJT=KTYPE(JSH)
      MINJ=KMIN(JSH)
      MAXJ=KMAX(JSH)
      NUMJ=MAXJ-MINJ+1
      LOCJ=KLOC(JSH)-MINJ
      SPJ=LJT.EQ.2.AND.MINJ.EQ.1
      SPIJ=SPI.OR.SPJ
      EXPNDI=LIT.GE.LJT
C
C     ----- KSHELL -----
C
      LKT=KTYPE(KSH)
      MINK=KMIN(KSH)
      MAXK=KMAX(KSH)
      NUMK=MAXK-MINK+1
      LOCK=KLOC(KSH)-MINK
      SPK=LKT.EQ.2.AND.MINK.EQ.1
C
C     ----- LSHELL -----
C
      LLT=KTYPE(LSH)
      MINL=KMIN(LSH)
      MAXL=KMAX(LSH)
      NUML=MAXL-MINL+1
      LOCL=KLOC(LSH)-MINL
      SPL=LLT.EQ.2.AND.MINL.EQ.1
      SPKL=SPK.OR.SPL
      SPIJKL=SPIJ.OR.SPKL
      EXPNDK=LKT.GE.LLT
C
      RETURN
      END
c--------------------------------------------------------------------------
C*MODULE GRD2A   *DECK JKDMEM
      SUBROUTINE xJKDMEM(MDER,LOADFM,IADDR,LENGTH,MINXYZ,MAXXYZ,MINVEC,
     *                  POPLE,MP2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "gamess.par"
C
      LOGICAL SOME,OUT,DBG,NORM,POPLE,MP2,GOPARR,DSKWRK,MASWRK,
     *        VECTOR,ASKVEC
C
      DIMENSION LENSHL(5)
C
      COMMON /xDSHLT / RTOL,DTOL,VTOL1,VTOL2,VTOLS,OUT,DBG
      COMMON /xIOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /xMACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /xOUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /xPAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON/xSHLBAS/MAXTYP,MAXNUM
      COMMON/xSHLNRM/PNRM(35)
      COMMON/xDERPAR/IDER,JDER,KDER,LDER,NDER
      COMMON/xDERMEM/IWFN,IXCH,INIJG,IGINT,IFINT,ISINT,IIJKLG,
     1 IDAB,ICHRG,IXY,IXZ,IYZ,IX,IY,IZ,ISJ,ISK,ISL,IGIJKL,IGNKL,IGNM,
     2 IDIJ,IDKL,IB00,IB01,IB10,IC00,ID00,IFXX,
     3 IAAI,IAAJ,IBBK,IBBL,IFI,IFJ,IFK,IFL,
     4 ISII,ISJJ,ISKK,ISLL,ISIJ,ISIK,ISIL,ISJK,ISJL,ISKL,
     5 IDIJSI,IDIJSJ,IDKLSK,IDKLSL,IABV,ICV,IRW
C
cwc1
      common /gms_deriv_ctrl/ want_integrals, all_centers
      logical want_integrals, all_centers
cwc0
      PARAMETER (LENVEC=255)
      PARAMETER(ONE=1.0D+00,SQRT3=1.73205080756888D+00)
      PARAMETER(SQRT5=2.23606797749979D+00,SQRT7=2.64575131106459D+00)
      DATA LENSHL /1,4,10,20,35/
C
      SOME = MASWRK  .AND.  NPRINT.NE.-5
C
C      FIRST FIND OUT WHETHER OR NOT THIS IS A VECTOR MACHINE, IF SO
C      THE VECTOR CODE WILL BE USED FOR LONG VECTORS
C
cwc      VECTOR=ASKVEC()
cwc1
      vector=.false.
cwc0
C
C       FIND OUT HOW MUCH MEMORY IS AVAILABLE
C
      CALL yGOTFM(NGOTMX)
C
      MAXVEC=LENVEC/3
      IF(     VECTOR) MINVEC=24
      IF(.NOT.VECTOR) MINVEC=MAXVEC*3+1
C
C     ----- THIS PROGRAM FOR DERIVATIVES -----
C
      NDER=MDER
      IF(NDER.NE.1.AND.NDER.NE.2) NDER=0
      NDER0=0
      NDER1=0
      NDER2=0
cwc      IF(OUT              ) NDER0=1
      IF(OUT.or.want_integrals              ) NDER0=1
cwc      IF(OUT.AND.NDER.GE.1) NDER1=12
      IF((OUT.or.want_integrals).AND.NDER.GE.1) NDER1=12
cwc      IF(OUT.AND.NDER.EQ.2) NDER2=78
      IF((OUT.or.want_integrals).AND.NDER.EQ.2) NDER2=78
      IF(        NDER.EQ.2) NDER1=12
C
C     ----- CHECK MAXIMUM ANGULAR MOMENTUM -----
C
      CALL xBASCHK(LMAX)
      MAXTYP=LMAX+1
cwc      IF (MAXTYP.GT.5) THEN
cwc         WRITE(IW,9999)
cwc         CALL ABRT
cwc      END IF
      MINXYZ=(4*MAXTYP -2 +NDER)/2
C
C     ----- GET NUMBER OF PRIMITIVE CHARGE DISTRIBUTIONS -----
C
      NIJG=0
      DO 20 II=1,NSHELL
        DO 20 JJ=1,II
          NIJG=NIJG+KNG(II)*KNG(JJ)
   20 CONTINUE
C
C     -----  AT THIS POINT IT IS GOOD TO REMEMBER THAT    -----
C            -MAXTYP- = HIGHEST SHELL ANGULAR MOMENTUM
C            -MAXFUN- = NUMBER OF FUNCTIONS WITH ANGULAR
C                       MOMENTUM LESS OR EQUAL TO -MAXTYP-
C            -MAXNUM- = NUMBER OF FUNCTIONS WITH ANGULAR
C                       MOMENTUM         EQUAL TO -MAXTYP-
C            -MAXXYZ- = MAXIMUM NUMBER OF PRIMITIVE INTEGRALS
C                       THAT CAN BE HANDLED IN ONE -VECTOR-
C            -NUMXYZ- = ACTUAL MAXIMUM LENGTH OF ONE -VECTOR-
C            -MAXXYZ- = IT IS NUMXYZ/3 . SINCE THE X, Y, AND Z
C                       COMPONENTS ARE TREATED AS A SINGLE VECTOR,
C                       -MAXXYZ- IS THE NUMBER OF (PRIMITIVE-ROOTS)
C                       COMBINATIONS WHICH CAN BE TREATED IN ONE
C                       VECTOR. FOR -SSSS- INTEGRALS WHICH REQUIRE
C                       ONE RYS ROOT, MAXXYZ HAPPENS TO COINCIDE WITH
C                       THE NUMBER OF PRIMITIVE INTEGRALS TREATED IN
C                       ONE VECTOR. FOR -DDDD- INTEGRALS WHICH
C                       REQUIRE FIVE RYS ROOTS, THE NUMBER OF PRIMITIVE
C                       INTEGRALS TREATED IN ONE VECTOR IS -MAXXYZ-/5 .
C
C
C     ----- SET NORMALIZATION CONSTANTS -----
C
      MAXFUN=LENSHL(MAXTYP)
      DO 100 I=1,MAXFUN
  100 PNRM(I)=ONE
      NORM=NORMF.NE.1.OR.NORMP.NE.1
      IF(.NOT.NORM) GO TO 180
C
      SQRT53=SQRT5/SQRT3
      DO 170 I=1,MAXFUN
        IF (I.EQ.1.OR.I.EQ.2.OR.I.EQ.5.OR.I.EQ.11.OR.I.EQ.21) THEN
           FI = ONE
        ELSE IF (I.EQ.8.OR.I.EQ.20.OR.I.EQ.33) THEN
           FI=FI*SQRT3
        ELSE IF (I.EQ.14) THEN
           FI=FI*SQRT5
        ELSE IF (I.EQ.24) THEN
           FI=FI*SQRT7
        ELSE IF (I.EQ.30) THEN
           FI=FI*SQRT53
        END IF
        PNRM(I)=FI
  170 CONTINUE
C
C        ALWAYS ALLOCATE -L- SHELLS IF -P- IS HIGHEST ANGULAR MOMENTUM
C
  180 CONTINUE
      IDUM=MAXTYP-1
      MAXNUM=((IDUM+1)*(IDUM+2))/2
      IF(MAXNUM.LT.4) MAXNUM=4
      NGIJKL=(MAXNUM**4)
C
C     ----- FOR DERIVATIVES -----
C
      MODTYP=MAXTYP+NDER
C
C     ----- CALCULATE VECTOR LENGTH AND SET CORE POINTERS -----
C
      LVAR=0
      LFIX=LENGTH
      LFIX=LFIX  +(NSHELL*(NSHELL+1))/2
C
C     ----- -SP- FUNCTIONS FOR FIRST DERIVATIVES ARE SPECIAL -----
C
      IWFN  = LOADFM + 1
      IXCH  = IWFN + LENGTH
      IDAB  = IXCH + (NSHELL*(NSHELL+1))/2
      ILAST = IDAB + NGIJKL
      INEED = ILAST- LOADFM
      IADDR = ILAST
      IF (POPLE.AND.MAXTYP.LT.3.AND.NDER.EQ.1) GO TO 300
C
      LFIX=LFIX+( (NSHELL*(NSHELL+1))/2 )*2
      LFIX=LFIX+NGIJKL*(NDER0+NDER1+NDER2)
      LFIX=LFIX+NGIJKL* 4
      LFIX=LFIX+NGIJKL
      LFIX=LFIX+NIJG*15
      LVAR=     ( MODTYP**2       * MODTYP**2       )*3
      LVAR=LVAR+( MODTYP**2       *(MODTYP+MODTYP-1))*3
      LVAR=LVAR+((MODTYP+MODTYP-1)*(MODTYP+MODTYP-1))*3
      LVAR=LVAR+( MODTYP**2                         )*3
      LVAR=LVAR+((MODTYP+MODTYP-1)                  )*3
      LVAR=LVAR+((MODTYP+MODTYP-1)* 3               )*3
      LVAR=LVAR+(  3                                )*3
      LVAR=LVAR+(  9                                )
      LVAR=LVAR+(  4                                )
      LVAR=LVAR+(  5                                )
      LVAR=LVAR+( 18                                )
      LVAR=LVAR+(  2                                )
      LVAR=LVAR+(  4                                )*3
      LVAR=LVAR+( MODTYP**2       * MODTYP**2       )*3*14
C
      MAXXYZ=(NGOTMX-LFIX-1)/LVAR
      IF(MP2) MAXXYZ = MIN(MAXXYZ,2*MINXYZ)
      MINMEM = (MINXYZ*LVAR)+1+LFIX
cwc1
      write (96,'(a,i9,a)') 
     >     'minimum memory for 2-el gradient:',minmem,' words'
cwc0
      IF(MAXXYZ.LT.MINXYZ) THEN
cwc1
         write (*,'(a)') '@XJKDMEM: not enough memory'
         call aces_exit(1)
cwc0
cwc         IF(MASWRK) WRITE(IW,9998) MINMEM
cwc         CALL ABRT
      END IF
cwc      IF(.NOT.MP2  .AND.  SOME) WRITE(IW,9994) MINMEM
      IF(MAXXYZ.GT.MAXVEC) MAXXYZ=MAXVEC
C
C     X(IWFN  ) = WAVEFUNCTION DATA
C     X(IXCH  ) = EXCHANGE INTEGRAL THRESHOLD
C     X(INIJG ) = CHARGE DISTRIBUTION POINTERS
C     X(IGINT ) = ELECTRON REPULSION INTEGRALS
C     X(IFINT ) = FIRST DERIVATIVE INTEGRALS
C     X(ISINT ) = SECOND DERIVATIVE INTEGRALS
C     X(IIJKLG) = INDICES
C     X(IDAB  ) = DENSITY ARRAY -DAB-
C     X(ICHRG ) = CHARGE DISTRIBUTION PARAMETERS
C     X(IXY   ) = TEMPORARY ARRAY WHEN FORMING DERIVATIVE INTEGRALS
C     X(IXZ   ) = TEMPORARY ARRAY WHEN FORMING DERIVATIVE INTEGRALS
C     X(IYZ   ) = TEMPORARY ARRAY WHEN FORMING DERIVATIVE INTEGRALS
C     X(IX    ) = TEMPORARY ARRAY WHEN FORMING DERIVATIVE INTEGRALS
C     X(IY    ) = TEMPORARY ARRAY WHEN FORMING DERIVATIVE INTEGRALS
C     X(IZ    ) = TEMPORARY ARRAY WHEN FORMING DERIVATIVE INTEGRALS
C     X(ISJ   ) = TEMPORARY ARRAY WHEN -SP- SHELLS
C     X(ISK   ) = TEMPORARY ARRAY WHEN -SP- SHELLS
C     X(ISL   ) = TEMPORARY ARRAY WHEN -SP- SHELLS
C     X(IGIJKL) = ( 2-D , 4 CENTERS ) INTEGRALS
C     X(IGNKL ) = ( 2-D , 3 CENTERS ) INTEGRALS
C     X(IGNM  ) = ( 2-D , 2 CENTERS ) INTEGRALS
C     X(IDIJ  ) = CONTRACTION DENSITY FOR -IJ- CHARGE DISTRIBUTION
C     X(IDKL  ) = CONTRACTION DENSITY FOR -KL- CHARGE DISTRIBUTION
C     X(IB00  ) = -B00-
C     X(IB01  ) = -B01-
C     X(IB10  ) = -B10-
C     X(IC00  ) = -C00-
C     X(ID00  ) = -D00-
C     X(IFXX  ) = -F00-
C     X(IDIJSI) = SCALING FACTOR FOR -S- FUNCTION OF AN -SP- II SHELL
C     X(IDIJSJ) = SCALING FACTOR FOR -S- FUNCTION OF AN -SP- JJ SHELL
C     X(IDKLSK) = SCALING FACTOR FOR -S- FUNCTION OF AN -SP- KK SHELL
C     X(IDKLSL) = SCALING FACTOR FOR -S- FUNCTION OF AN -SP- LL SHELL
C     X(IABV  ) = -AB- VECTOR FOR PRIMITIVE INTEGRALS
C     X(ICV   ) = -CV- VECTOR FOR PRIMITIVE INTEGRALS
C     X(IRW   ) = -RW- VECTOR FOR RYS ROOTS AND WEIGHTS
C     X(IAAI  ) = EXPONENT FOR DERIVATIVE OF II SHELL
C     X(IAAJ  ) = EXPONENT FOR DERIVATIVE OF JJ SHELL
C     X(IBBK  ) = EXPONENT FOR DERIVATIVE OF KK SHELL
C     X(IBBL  ) = EXPONENT FOR DERIVATIVE OF LL SHELL
C     X(IFI   ) = FIRST DERIVATIVE WRT. II OF ( 2-D , 4 CENTERS ) INT.
C     X(IFJ   ) = FIRST DERIVATIVE WRT. JJ OF ( 2-D , 4 CENTERS ) INT.
C     X(IFK   ) = FIRST DERIVATIVE WRT. KK OF ( 2-D , 4 CENTERS ) INT.
C     X(IFL   ) = FIRST DERIVATIVE WRT. LL OF ( 2-D , 4 CENTERS ) INT.
C     X(ISII  ) = SECOND DER. WRT II AND II OF ( 2-D , 4 CENTERS ) INT.
C     X(ISJJ  ) = SECOND DER. WRT JJ AND JJ OF ( 2-D , 4 CENTERS ) INT.
C     X(ISKK  ) = SECOND DER. WRT KK AND KK OF ( 2-D , 4 CENTERS ) INT.
C     X(ISLL  ) = SECOND DER. WRT LL AND LL OF ( 2-D , 4 CENTERS ) INT.
C     X(ISIJ  ) = SECOND DER. WRT II AND JJ OF ( 2-D , 4 CENTERS ) INT.
C     X(ISIK  ) = SECOND DER. WRT II AND KK OF ( 2-D , 4 CENTERS ) INT.
C     X(ISIL  ) = SECOND DER. WRT II AND LL OF ( 2-D , 4 CENTERS ) INT.
C     X(ISJK  ) = SECOND DER. WRT JJ AND KK OF ( 2-D , 4 CENTERS ) INT.
C     X(ISJL  ) = SECOND DER. WRT JJ AND LL OF ( 2-D , 4 CENTERS ) INT.
C     X(ISKL  ) = SECOND DER. WRT KK AND LL OF ( 2-D , 4 CENTERS ) INT.
C
      INIJG =IXCH  + (NSHELL*(NSHELL+1))/2
      IGINT =INIJG +(((NSHELL*(NSHELL+1))/2)*2)/NWDVAR
      IFINT =IGINT +  NGIJKL*NDER0
      ISINT =IFINT +  NGIJKL*NDER1
      IIJKLG=ISINT +  NGIJKL*NDER2
      IDAB  =IIJKLG+  NGIJKL*4/NWDVAR
      ICHRG =IDAB  +  NGIJKL
      IXY   =ICHRG +  NIJG  *15
      IXZ   =IXY   +(  1                                )*MAXXYZ
      IYZ   =IXZ   +(  1                                )*MAXXYZ
      IX    =IYZ   +(  1                                )*MAXXYZ
      IY    =IX    +(  1                                )*MAXXYZ
      IZ    =IY    +(  1                                )*MAXXYZ
      ISJ   =IZ    +(  1                                )*MAXXYZ
      ISK   =ISJ   +(  1                                )*MAXXYZ
      ISL   =ISK   +(  1                                )*MAXXYZ
      IGIJKL=ISL   +(  1                                )*MAXXYZ
      IGNKL =IGIJKL+( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      IGNM  =IGNKL +( MODTYP**2       *(MODTYP+MODTYP-1))*MAXXYZ*3
      IDIJ  =IGNM  +((MODTYP+MODTYP-1)*(MODTYP+MODTYP-1))*MAXXYZ*3
      IDKL  =IDIJ  +( MODTYP**2                         )*MAXXYZ*3
      IB00  =IDKL  +((MODTYP+MODTYP-1)                  )*MAXXYZ*3
      IB01  =IB00  +((MODTYP+MODTYP-1)                  )*MAXXYZ*3
      IB10  =IB01  +((MODTYP+MODTYP-1)                  )*MAXXYZ*3
      IC00  =IB10  +((MODTYP+MODTYP-1)                  )*MAXXYZ*3
      ID00  =IC00  +(  1                                )*MAXXYZ*3
      IFXX  =ID00  +(  1                                )*MAXXYZ*3
      IDIJSI=IFXX  +(  1                                )*MAXXYZ*3
      IDIJSJ=IDIJSI+(  1                                )*MAXXYZ
      IDKLSK=IDIJSJ+(  1                                )*MAXXYZ
      IDKLSL=IDKLSK+(  1                                )*MAXXYZ
      IABV  =IDKLSL+(  1                                )*MAXXYZ
      ICV   =IABV  +(  5                                )*MAXXYZ
      IRW   =ICV   +( 18                                )*MAXXYZ
      IAAI  =IRW   +(  2                                )*MAXXYZ
      IAAJ  =IAAI  +(  1                                )*MAXXYZ*3
      IBBK  =IAAJ  +(  1                                )*MAXXYZ*3
      IBBL  =IBBK  +(  1                                )*MAXXYZ*3
      IFI   =IBBL  +(  1                                )*MAXXYZ*3
      IFJ   =IFI   +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      IFK   =IFJ   +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      IFL   =IFK   +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      ISII  =IFL   +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      ISJJ  =ISII  +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      ISKK  =ISJJ  +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      ISLL  =ISKK  +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      ISIJ  =ISLL  +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      ISIK  =ISIJ  +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      ISIL  =ISIK  +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      ISJK  =ISIL  +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      ISJL  =ISJK  +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      ISKL  =ISJL  +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      ILAST =ISKL  +( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      INEED =ILAST - LOADFM
      IADDR =ILAST

C
cwc      IF(.NOT.MP2  .AND.  SOME) WRITE(IW,9997) INEED
      IF(OUT) WRITE(IW,9996) IWFN,IXCH,INIJG,IGINT,IFINT,ISINT,IIJKLG,
     1 IDAB,ICHRG,IXY,IXZ,IYZ,IX,IY,IZ,ISJ,ISK,ISL,IGIJKL,IGNKL,IGNM,
     2 IDIJ,IDKL,IB00,IB01,IB10,IC00,ID00,IFXX,
     3 IAAI,IAAJ,IBBK,IBBL,IFI,IFJ,IFK,IFL,
     4 ISII,ISJJ,ISKK,ISLL,ISIJ,ISIK,ISIL,ISJK,ISJL,ISKL,
     5 IDIJSI,IDIJSJ,IDKLSK,IDKLSL,IABV,ICV,IRW
      RETURN
C
  300 CONTINUE
      IF(OUT) WRITE(IW,9995) IWFN,IXCH,INEED,NGOTMX
      RETURN
C
 9999 FORMAT(' GRADIENTS ARE LIMITED TO G AND LOWER FUNCTIONS ')
 9998 FORMAT(/,' NOT ENOUGH MEMORY FOR THE TWO-ELECTRON GRADIENT',/,
     1         ' YOU WILL NEED AT LEAST ',I10,' WORDS.')
 9997 FORMAT(' USING ',I10,' WORDS OF MEMORY.')
 9996 FORMAT(
     1 ' IWFN  ',I8,' IXCH  ',I8,' INIJG ',I8,' IGINT ',I8,' IFINT ',I8,
     1 ' ISINT ',I8,' IIJKLG',I8,/,
     1 ' IDAB  ',I8,' ICHRG ',I8,' IXY   ',I8,' IXZ   ',I8,' IYZ   ',I8,
     1 ' IX    ',I8,' IY    ',I8,' IZ    ',I8,/,
     1 ' ISJ   ',I8,' ISK   ',I8,' ISL   ',I8,
     1 ' IGIJKL',I8,' IGNKL ',I8,' IGNM  ',I8,/,
     2 ' IDIJ  ',I8,' IDKL  ',I8,' IB00  ',I8,' IB01  ',I8,' IB10  ',I8,
     3 ' IFI   ',I8,' IFJ   ',I8,' IFK   ',I8,' IFL   ',I8,/,
     4 ' ISII  ',I8,' ISJJ  ',I8,' ISKK  ',I8,' ISLL  ',I8,
     4 ' ISIJ  ',I8,' ISIK  ',I8,' ISIL  ',I8,' ISJK  ',I8,/
     4 ' ISJL  ',I8,' ISKL  ',I8,/,
     5 ' IDIJSI',I8,' IDIJSJ',I8,' IDKLSK',I8,' IDKLSL',I8,/,
     5 ' IABV  ',I8,' ICV   ',I8,' IRW   ',I8)
 9995 FORMAT(' SPECIAL -SP- ROUTINES USED.',
     1       ' IWFN,IXCH,INEED,MAXFM = ',4I8)
 9994 FORMAT(/,' THE MINIMUM MEMORY REQUIRED FOR THIS RUN IS ',I10,
     1         ' WORDS.')
      END
c--------------------------------------------------------------------------
C*MODULE GRD2A   *DECK DDEBUT
      SUBROUTINE xDDEBUT(IA,DA,DB,L1,L2,L3,NOCORE,NOOPEN,HFSCF,
     *                  UHFTYP,ROGVB)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "gamess.par"
C
      LOGICAL NOCORE,NOOPEN
      LOGICAL HFSCF,UHFTYP,ROGVB
C
      DIMENSION DA(L2),DB(L2),IA(*)
C
      COMMON /xIOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
cwc      COMMON /xGVBWFN/ CICOEF(2,12),F(25),ALPHA(325),BETA(325),NO(10),
cwc     *                NCO,NSETO,NOPEN,NPAIR,NMOGVB,NCONF(MXAO),NHAM
C
C     ---- READ IN 1 BODY DENSITY MATRICES FOR RHF,UHF,ROHF ----
C     DA = TOTAL (ALPHA+BETA) DENSITY, DB = SPIN (ALPHA-BETA) DENSITY
C
      IF(HFSCF) THEN
         CALL yDAREAD(IDAF,IODA,DA,L2,16,0)
         IF(UHFTYP) THEN
            CALL yDAREAD(IDAF,IODA,DB,L2,20,0)
            DO 120 I = 1,L2
               DUMA = DA(I)
               DUMB = DB(I)
               DA(I) = DUMA+DUMB
               DB(I) = DUMA-DUMB
  120       CONTINUE
         END IF
      END IF
C
C     ----- SET UP CORE DENSITY MATRIX IN -DA- FOR GVB -----
C     SNEAKY, THE MO VECTORS ARE STORED IN SQUARE DB.
C
cwc      IF (ROGVB) THEN
cwc         CALL DAREAD(IDAF,IODA,DB,L3,15,0)
cwc         CALL DENCOR(DA,DB,IA,L1)
cwc         NOCORE = NCO .EQ. 0
cwc         NOOPEN = NSETO .EQ. 0 .AND. NPAIR .EQ. 0
cwc      END IF
C
      RETURN
      END
c--------------------------------------------------------------------------
C*MODULE GRD2A   *DECK JKDSET
      SUBROUTINE xJKDSET
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "gamess.par"
C
      LOGICAL SP
C
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /xSYMTRY/ MAPSHL(MXSH,48),MAPCTR(MXATM,48),
     *                T(432),INVT(48),NT
C
      DIMENSION MI(48)
      DATA NTYP  /6/
C
      DO 270 IITYP=1,NTYP
        DO 260 II=1,NSHELL
          LIT=KTYPE(II)
          SP =KTYPE(II).EQ.2.AND.KMIN(II).EQ.1
          IF(SP) LIT=NTYP
          IF(LIT.NE.IITYP) GO TO 260
          DO 210 IT=1,NT
            ID=MAPSHL(II,IT)
            IF( ID.GT.II   ) GO TO 260
  210     MI(IT)=ID
C
          IF(NT.EQ.1) GO TO 260
          DO 230 IT=2,NT
            MAX=IT-1
            DO 220 JT=1,MAX
              IF(MI(JT).NE.MI(IT)) GO TO 220
              MI(IT)=0
  220       CONTINUE
  230     CONTINUE
C
  260   CONTINUE
  270 CONTINUE
C
      RETURN
      END
c--------------------------------------------------------------------------
C*MODULE GRD2A   *DECK OEDHND
      SUBROUTINE xOEDHND(NIJG,DCHRG)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "gamess.par"
C
      COMMON /xDSHLNO/ LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     *                NIJ,IJD,KLD,IJ,KL
      COMMON /xIJPAIR/ IA(MXAO)
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      DIMENSION NIJG(2,*), DCHRG(15,*)
C
C     ----- ONE-ELECTRON CHARGE DISTRIBUTION AND DERIVATIVES -----
C
      NIJ0=0
      DO 9000 II=1,NSHELL
        DO 8000 JJ=1,II
          ISH=II
          JSH=JJ
          CALL xOEDSHL(ISH,JSH,DCHRG,NIJ0)
          IIJJ=IA(MAX0(II,JJ))+MIN0(II,JJ)
          NIJG(1,IIJJ)=NIJ0
          NIJG(2,IIJJ)=NIJ
          NIJ0=NIJ0+NIJ
 8000   CONTINUE
 9000 CONTINUE
C
      RETURN
      END
c--------------------------------------------------------------------------
C*MODULE GRD2A   *DECK DFINAL
      SUBROUTINE xDFINAL(INDEX)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
#include "gamess.par"
C
      COMMON /xFUNCT / EHF,EG(3*MXATM)
      COMMON /xGRAD  / DE(3*MXATM)
      COMMON /xINFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM)
      COMMON /xIOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /xPAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /xRESTAR/ TIMLIM,IREST,NREC,INTLOC,IST,JST,KST,LST
cwc      COMMON /xTMVALS/ TI,TX,TIM
C
C     ----- SAVE INCOMPLETE GRADIENT -----
C
      IF (INDEX .EQ. 1) GO TO 100
      IREST = 4
      CALL yDAWRIT(IDAF,IODA,DE,3*NAT,3,0)
cwc      CALL TSECND(TIM)
cwc      IF (TIM .GE. TIMLIM .AND. MASWRK) WRITE (IW,9008)
cwc     *    IREST,IST,JST,KST,LST
      RETURN
C
C     ----- GRADIENT VECTOR HAS BEEN COMPLETED -----
C
  100 CONTINUE
      NCOORD = 3*NAT
      CALL XCOPY(NCOORD,DE,1,EG,1)
      CALL yDAWRIT(IDAF,IODA,EG,NCOORD,3,0)
      IREST = 5
      IST = 1
      JST = 1
      KST = 1
      LST = 1
      RETURN
C
 9008 FORMAT(1X,'WARNING.......THIS JOB MUST BE RESTARTED .....'/
     *       1X,'ENTER IREST=',I5,' IN $CONTRL'/
     *       1X,'ENTER IST,JST,KST,LST=',4I5,' IN $INTGRL')
      END
c--------------------------------------------------------------------------
c---------------- Rys routines --------------------------------------------
c--------------------------------------------------------------------------
C*MODULE GRD2C   *DECK JKDSPD
      SUBROUTINE xJKDSPD(NIJ0,NKL,NKL0,DCHRG,GINT,FINT,SINT,IJKLG,
     1                  GIJKL,GNKL,GNM,XY,XZ,YZ,X,Y,Z,SJ,SK,SL,
     2                  B00,B01,B10,C00,D00,F00,DIJ,DKL,DIJSI,DIJSJ,
     3                  DKLSK,DKLSL,ABV,CV,RWV,AAI,AAJ,BBK,BBL,FI,FJ,
     4                  FK,FL,SII,SJJ,SKK,SLL,SIJ,SIK,SIL,SJK,SJL,SKL,
     5                  DAB,MAXXYZ,FC,NC,DF,LDF,NBF,DEN,QQ4,MINVEC,
     6                  DABCUT, DABMAX)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SKIPI,SKIPJ,SKIPK,SKIPL
      LOGICAL SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      LOGICAL NMAXS,NMAXP,MMAXS,MMAXP
      LOGICAL EXPNDI,EXPNDK
      LOGICAL FIRST
      LOGICAL LAST
      LOGICAL OUT,DBG
      COMMON /xDSHLNO/ LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     *                NIJ,IJD,KLD,IJ,KL
      COMMON /xDSHLT / RTOL,DTOL,VTOL1,VTOL2,VTOLS,OUT,DBG
      COMMON /xROOT  / XX,U(13),W(13),NROOTS
      COMMON/xSHLGNM/NMAXS,NMAXP,MMAXS,MMAXP
      COMMON/xSHLTYP/SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      COMMON/xSHLXPN/EXPNDI,EXPNDK
      COMMON/xDERPAR/IDER,JDER,KDER,LDER,NDER
      COMMON/xDERSKP/IIAT,JJAT,KKAT,LLAT,SKIPI,SKIPJ,SKIPK,SKIPL
      DIMENSION DCHRG(15,*)
      DIMENSION GINT(1),FINT(12,1),SINT(78,1),IJKLG(4,1)
      DIMENSION XY(1),XZ(1),YZ(1),X(1),Y(1),Z(1),SJ(1),SK(1),SL(1)
      DIMENSION GIJKL(1),GNKL(1),GNM(1)
      DIMENSION DIJ(1),DKL(1),DIJSI(1),DIJSJ(1),DKLSK(1),DKLSL(1)
      DIMENSION FI(1),FJ(1),FK(1),FL(1)
      DIMENSION SII(1),SJJ(1),SKK(1),SLL(1)
      DIMENSION SIJ(1),SIK(1),SIL(1),SJK(1),SJL(1),SKL(1)
      DIMENSION B00(1),B01(1),B10(1),C00(1),D00(1),F00(1)
      DIMENSION AAI(1),AAJ(1),BBK(1),BBL(1)
      DIMENSION ABV(5,1),CV(18,1),RWV(2,1)
      DIMENSION DAB(1)
      DIMENSION FC(NC,1),DF(LDF,1),DEN(1)
      PARAMETER(PI252=34.986836655250D+00, ONE=1.0D+00)
C
      Q4=PI252*QQ4
C
      NIMAX=LIT + IDER
      NJMAX=LJT + JDER
      NKMAX=LKT + KDER
      NLMAX=LLT + LDER
      NMAX=LIT+LJT-1 + MIN0(IDER+JDER,NDER)
      MMAX=LKT+LLT-1 + MIN0(KDER+LDER,NDER)
      NMAXS=NMAX.EQ.1
      NMAXP=NMAX.LE.2
      MMAXS=MMAX.EQ.1
      MMAXP=MMAX.LE.2
C
      MAXG=MAXXYZ/NROOTS
C
C     ----- PAIR OF K,L PRIMITIVES -----
C
      FIRST=.TRUE.
      NG=0
      KLG=0
  100 KLG=KLG+1
      IF(KLG.GT.NKL) GO TO 300
      DB=DCHRG( 1,KLG+NKL0)
      BB=DCHRG( 2,KLG+NKL0)
      XB=DCHRG( 3,KLG+NKL0)
      YB=DCHRG( 4,KLG+NKL0)
      ZB=DCHRG( 5,KLG+NKL0)
      XD=DCHRG( 6,KLG+NKL0)
      YD=DCHRG( 7,KLG+NKL0)
      ZD=DCHRG( 8,KLG+NKL0)
      DXKL=DCHRG( 9,KLG+NKL0)
      DYKL=DCHRG(10,KLG+NKL0)
      DZKL=DCHRG(11,KLG+NKL0)
      Q4DB=Q4*DB
C
C     ----- PAIR OF I,J PRIMITIVES -----
C
      IJG=0
  200 IJG=IJG+1
      IF(IJG.GT.NIJ) GO TO 100
      DA=DCHRG( 1,IJG+NIJ0)
      AA=DCHRG( 2,IJG+NIJ0)
      XA=DCHRG( 3,IJG+NIJ0)
      YA=DCHRG( 4,IJG+NIJ0)
      ZA=DCHRG( 5,IJG+NIJ0)
      AANDB1=ONE/(AA+BB)
      Q4DBDA=Q4DB*DA
      DUM   =Q4DBDA*Q4DBDA*AANDB1
      IF(DUM.LE.DTOL) GO TO 200
      Q4DBDA=Q4DBDA* SQRT(AANDB1)
      IF( ABS(Q4DBDA*DABMAX).LT.DABCUT) GO TO 200
      RHO   =AA*BB*AANDB1
      XX=RHO*((XA-XB)**2+(YA-YB)**2+(ZA-ZB)**2)
C
      NG=NG+1
      ABV(1,NG)=AA
      ABV(2,NG)=BB
      ABV(3,NG)=RHO
      ABV(4,NG)=Q4DBDA
      ABV(5,NG)=XX
C
      XC=DCHRG( 6,IJG+NIJ0)
      YC=DCHRG( 7,IJG+NIJ0)
      ZC=DCHRG( 8,IJG+NIJ0)
      DXIJ=DCHRG( 9,IJG+NIJ0)
      DYIJ=DCHRG(10,IJG+NIJ0)
      DZIJ=DCHRG(11,IJG+NIJ0)
C
      AAI(NG)=DCHRG(12,IJG+NIJ0)
      AAJ(NG)=DCHRG(13,IJG+NIJ0)
      BBK(NG)=DCHRG(12,KLG+NKL0)
      BBL(NG)=DCHRG(13,KLG+NKL0)
C
      IF(MMAXS) GO TO 210
      CV( 1,NG)=AA*(XA-XD)
      CV( 2,NG)=BB*(XB-XD)
      CV( 3,NG)=AA*(YA-YD)
      CV( 4,NG)=BB*(YB-YD)
      CV( 5,NG)=AA*(ZA-ZD)
      CV( 6,NG)=BB*(ZB-ZD)
  210 IF(NMAXS) GO TO 220
      CV( 7,NG)=AA*(XA-XC)
      CV( 8,NG)=BB*(XB-XC)
      CV( 9,NG)=AA*(YA-YC)
      CV(10,NG)=BB*(YB-YC)
      CV(11,NG)=AA*(ZA-ZC)
      CV(12,NG)=BB*(ZB-ZC)
  220 CONTINUE
      CV(13,NG)=DXIJ
      CV(14,NG)=DYIJ
      CV(15,NG)=DZIJ
      CV(16,NG)=DXKL
      CV(17,NG)=DYKL
      CV(18,NG)=DZKL
      IF(SPI) DIJSI(NG)=DCHRG(14,IJG+NIJ0)
      IF(SPJ) DIJSJ(NG)=DCHRG(15,IJG+NIJ0)
      IF(SPK) DKLSK(NG)=DCHRG(14,KLG+NKL0)
      IF(SPL) DKLSL(NG)=DCHRG(15,KLG+NKL0)
C
      IF(NG.LT.MAXG) GO TO 200
      LAST=.FALSE.
      GO TO 310
C
  300 CONTINUE
      LAST=.TRUE.
  310 CONTINUE
      NUMG=NG
      IF(NUMG.EQ.0) GO TO 1000
C
      IF(NROOTS.EQ.1) GO TO 480
      IF (SPI) THEN
          DO 410 IROOT=2,NROOTS
          DO 410 IG=1,NUMG
              DIJSI(IG+NUMG*(IROOT-1))=DIJSI(IG)
  410     CONTINUE
      ENDIF
      IF (SPJ) THEN
         DO 430 IROOT=2,NROOTS
         DO 430 IG=1,NUMG
            DIJSJ(IG+NUMG*(IROOT-1))=DIJSJ(IG)
  430    CONTINUE
      ENDIF
      IF (SPK) THEN
         DO 450 IROOT=2,NROOTS
         DO 450 IG=1,NUMG
            DKLSK(IG+NUMG*(IROOT-1))=DKLSK(IG)
  450    CONTINUE
      ENDIF
      IF (SPL) THEN
         DO 470 IROOT=2,NROOTS
         DO 470 IG=1,NUMG
            DKLSL(IG+NUMG*(IROOT-1))=DKLSL(IG)
  470    CONTINUE
      ENDIF
C
  480 IF(SKIPI) GO TO 520
      DO 510 IRXYZ=2,NROOTS*3
      DO 510 IG=1,NUMG
      AAI(IG+NUMG*(IRXYZ-1))=AAI(IG)
  510 CONTINUE
  520 IF(SKIPJ) GO TO 540
      DO 530 IRXYZ=2,NROOTS*3
      DO 530 IG=1,NUMG
      AAJ(IG+NUMG*(IRXYZ-1))=AAJ(IG)
  530 CONTINUE
  540 IF(SKIPK) GO TO 560
      DO 550 IRXYZ=2,NROOTS*3
      DO 550 IG=1,NUMG
      BBK(IG+NUMG*(IRXYZ-1))=BBK(IG)
  550 CONTINUE
  560 IF(SKIPL) GO TO 580
      DO 570 IRXYZ=2,NROOTS*3
      DO 570 IG=1,NUMG
      BBL(IG+NUMG*(IRXYZ-1))=BBL(IG)
  570 CONTINUE
  580 CONTINUE


C
C     ----- COMPUTE ROOTS AND WEIGHTS FOR QUADRATURE -----
C
      CALL xJKWRYS(RWV,ABV,NUMG)
C
C     ----- COMPUTE COEFFICIENTS FOR RECURSION FORMULAE -----
C
      CALL xJKBCDF(B00,B01,B10,C00,D00,F00,DIJ,DKL,
     1            ABV,CV,RWV,NUMG,NROOTS)
C
C     ----- COMPUTE -X- , -Y- , -Z- INTEGRALS ( 2 CENTERS, 2-D ) -----
C
      CALL xJKGNMS(GNM,NUMG*NROOTS*3,NMAX,MMAX,
     2             B00,B01,B10,C00,D00,F00)
C
C     ----- COMPUTE -X- , -Y- , -Z- INTEGRALS ( 4 CENTERS, 2-D ) -----
C
      CALL xJKXYZS(GIJKL,GIJKL,GNKL,GNKL,GNKL,GNM,GNM,
     2 NUMG*NROOTS*3,NMAX,MMAX,NIMAX,NJMAX,NKMAX,NLMAX,DIJ,DKL,
     3 EXPNDI,EXPNDK)
C
C     ----- COMPUTE -X- , -Y- , -Z- INTEGRALS FOR DERIVATIVES -----
C
      CALL xJDXYZS(GIJKL,GIJKL,GIJKL,GIJKL,GIJKL,GIJKL,GIJKL,GIJKL,
     2 NUMG*NROOTS*3,NIMAX,NJMAX,NKMAX,NLMAX,LIT,LJT,LKT,LLT,
     3 AAI,AAJ,BBK,BBL,
     4 FI,FJ,FK,FL,SII,SJJ,SKK,SLL,SIJ,SIK,SIL,SJK,SJL,SKL)
C
C     ----- ZERO OUT FIRST TIME AROUND -----
C
      IF(FIRST) THEN
         CALL xJKDZER(GINT,FINT,SINT,IJKLG)
         FIRST=.FALSE.
      ENDIF
C
C     ----- COMPUTE DERIVATIVE INTEGRALS -----
C
      CALL xDSPDFS(NUMG,NROOTS,IJKLG,GINT,FINT,SINT,
     2 GIJKL,FI,FJ,FK,FL,SII,SJJ,SKK,SLL,SIJ,SIK,SIL,SJK,SJL,SKL,
     3 XY,XZ,YZ,X,Y,Z,SJ,SK,SL,DIJSI,DIJSJ,DKLSK,DKLSL,DAB)

C
      IF(LAST) GO TO 1000
      NG=0
      GO TO 200
 1000 IF(NUMG.EQ.0.AND.FIRST) RETURN
C
C     ----- PROCESS DERIVATIVE INTEGRALS -----
C
      CALL xJKDINV(FINT,SINT,IJKLG,FC,NC,DF,LDF,NBF,DEN)
C
      IF(OUT) CALL xJKDOUT(GINT,FINT,SINT,IJKLG,DAB)
      RETURN
      END
c--------------------------------------------------------------------------
C*MODULE GRD2C   *DECK JKXYZS
      SUBROUTINE xJKXYZS(GIJKL,HIJKL,GNKL,HNKL,FNKL,GNM,HNM,
     1 NG,NMAX,MMAX,NIMAX,NJMAX,NKMAX,NLMAX,DIJ,DKL,EXPNDI,EXPNDK)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL EXPNDI,EXPNDK
      DIMENSION GIJKL(NG*NLMAX*NKMAX,NJMAX,NIMAX)
      DIMENSION HIJKL(NG*NLMAX*NKMAX*NJMAX,NIMAX)
      DIMENSION  GNKL(NG,NLMAX,NKMAX,NMAX)
      DIMENSION  HNKL(NG*NLMAX*NKMAX,NMAX)
      DIMENSION  FNKL(NG*NLMAX*NKMAX*NMAX)
      DIMENSION   GNM(NG,NMAX,MMAX)
      DIMENSION   DIJ(NG)
      DIMENSION   DKL(NG)
C
C     ----- G(N,K,L) -----
C
      IF (.NOT.EXPNDK) THEN
C
        DO 30 NK=1,NKMAX
          DO 10 NL=1,NLMAX
            DO 10  N=1,NMAX
              DO 10 IG=1,NG
   10           GNKL(IG,NL,NK,N)=GNM(IG,N,NL)
          IF(NK.EQ.NKMAX) GO TO 30
          MAX=MMAX-NK
          DO 20  M=1,MAX
            DO 20  N=1,NMAX
            DO 20 IG=1,NG
   20         GNM(IG,N,M)=DKL(IG)*GNM(IG,N,M)+GNM(IG,N,M+1)
   30   CONTINUE
C
      ELSE
C
        DO 70 NL=1,NLMAX
          DO 50 NK=1,NKMAX
            DO 50  N=1,NMAX
              DO 50 IG=1,NG
   50           GNKL(IG,NL,NK,N)=GNM(IG,N,NK)
          IF(NL.EQ.NLMAX) GO TO 70
          MAX=MMAX-NL
          DO 60  M=1,MAX
            DO 60  N=1,NMAX
              DO 60 IG=1,NG
   60           GNM(IG,N,M)=DKL(IG)*GNM(IG,N,M)+GNM(IG,N,M+1)
   70   CONTINUE
C
      ENDIF
C
C     ----- G(I,J,K,L) -----
C
      IF (.NOT.EXPNDI) THEN
C
        DO 130 NI=1,NIMAX
          call xcopy(NG*NLMAX*NKMAX*NJMAX, FNKL, 1, HIJKL(1,NI), 1)
          IF(NI.EQ.NIMAX) GO TO 130
          MAX=NMAX-NI
          DO 120  N=1,MAX
            DO 120 NK=1,NKMAX
              DO 120 NL=1,NLMAX
                DO 120 IG=1,NG
  120             GNKL(IG,NL,NK,N)=DIJ(IG)*GNKL(IG,NL,NK,N)+
     1                             GNKL(IG,NL,NK,N+1)
  130   CONTINUE
C
      ELSE
C
        DO 170 NJ=1,NJMAX
          DO NI = 1, NIMAX
             call xcopy(NG*NLMAX*NKMAX, GNKL(1,1,1,NI), 1, 
     *               GIJKL(1,NJ,NI), 1)      
          ENDDO

          IF (NJ.EQ.NJMAX) GO TO 170
          MAX=NMAX-NJ
          DO 160  N=1,MAX
            DO 160 NK=1,NKMAX
              DO 160 NL=1,NLMAX
                DO 160 IG=1,NG
  160             GNKL(IG,NL,NK,N)=DIJ(IG)*GNKL(IG,NL,NK,N)+
     1                             GNKL(IG,NL,NK,N+1)
  170   CONTINUE

      ENDIF
C
      RETURN
      END
c--------------------------------------------------------------------------
C*MODULE GRD2C   *DECK JKBCDF
      SUBROUTINE xJKBCDF(B00,B01,B10,C00,D00,F00,DIJ,DKL,
     1                  ABV,CV,RWV,NUMG,NROOTS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL NMAXS,NMAXP,MMAXS,MMAXP
      COMMON/xSHLGNM/NMAXS,NMAXP,MMAXS,MMAXP
      DIMENSION B00(NUMG,NROOTS,3),B01(NUMG,NROOTS,3),B10(NUMG,NROOTS,3)
      DIMENSION C00(NUMG,NROOTS,3)
      DIMENSION D00(NUMG,NROOTS,3)
      DIMENSION F00(NUMG,NROOTS,3),DIJ(NUMG,NROOTS,3),DKL(NUMG,NROOTS,3)
      DIMENSION ABV(5,1),CV(18,1)
      DIMENSION RWV(2,NUMG,NROOTS)
      PARAMETER (PT5=0.5D+00, ONE=1.0D+00)
C
      DO 40 NR=1,NROOTS
        DO 30 NG=1,NUMG
          AA =ABV(1,NG)
          BB =ABV(2,NG)
          RHO=ABV(3,NG)
          QAB=ABV(4,NG)
          UU =RHO*RWV(1,NG,NR)
          WW =    RWV(2,NG,NR)
          AAUU=AA+UU
          BBUU=BB+UU
          F00(NG,NR,1)=WW*QAB
          F00(NG,NR,2)=ONE
          F00(NG,NR,3)=ONE
          DUM2=PT5/(AA*BB+UU*(AA+BB))
          AUDUM=AAUU*DUM2
          BUDUM=BBUU*DUM2
           UDUM=  UU*DUM2
          B00(NG,NR,1)= UDUM
          B00(NG,NR,2)= UDUM
          B00(NG,NR,3)= UDUM
          B01(NG,NR,1)=AUDUM
          B01(NG,NR,2)=AUDUM
          B01(NG,NR,3)=AUDUM
          B10(NG,NR,1)=BUDUM
          B10(NG,NR,2)=BUDUM
          B10(NG,NR,3)=BUDUM
           UDUM= UDUM+ UDUM
          IF(MMAXS) GO TO 10
          AUDUM=AUDUM+AUDUM
          D00(NG,NR,1)= UDUM*CV( 1,NG) + AUDUM*CV( 2,NG)
          D00(NG,NR,2)= UDUM*CV( 3,NG) + AUDUM*CV( 4,NG)
          D00(NG,NR,3)= UDUM*CV( 5,NG) + AUDUM*CV( 6,NG)
   10     IF(NMAXS) GO TO 20
          BUDUM=BUDUM+BUDUM
          C00(NG,NR,1)= UDUM*CV( 8,NG) + BUDUM*CV( 7,NG)
          C00(NG,NR,2)= UDUM*CV(10,NG) + BUDUM*CV( 9,NG)
          C00(NG,NR,3)= UDUM*CV(12,NG) + BUDUM*CV(11,NG)
   20     CONTINUE
C
   30   CONTINUE
   40 CONTINUE
C
      DO 60 NR=1,NROOTS
        DO 50 NG=1,NUMG
          DIJ(NG,NR,1)=CV(13,NG)
          DIJ(NG,NR,2)=CV(14,NG)
          DIJ(NG,NR,3)=CV(15,NG)
          DKL(NG,NR,1)=CV(16,NG)
          DKL(NG,NR,2)=CV(17,NG)
          DKL(NG,NR,3)=CV(18,NG)
   50   CONTINUE
   60 CONTINUE
C
      RETURN
      END
c--------------------------------------------------------------------------
C*MODULE GRD2C   *DECK JKGNMS
      SUBROUTINE xJKGNMS(GNM,NG,NMAX,MMAX,B00,B01,B10,C00,D00,F00)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL NMAXS,NMAXP,MMAXS,MMAXP
C
      DIMENSION GNM(NG,NMAX,MMAX)
      DIMENSION C00(NG),D00(NG),F00(NG)
      DIMENSION B00(NG,1),B01(NG,1),B10(NG,1)
C
      COMMON /xSHLGNM/NMAXS,NMAXP,MMAXS,MMAXP
C
C     ----- G(0,0) -----
C
      DO 10 IG=1,NG
   10   GNM(IG,1,1)=F00(IG)
      IF(NMAXS.AND.MMAXS) RETURN
      IF(NMAXS) GO TO 30
C
C     ----- G(1,0) = C00 * G(0,0) -----
C
      DO 20 IG=1,NG
   20   GNM(IG,2,1)=C00(IG)*GNM(IG,1,1)
C
   30 IF(MMAXS) GO TO 60
C
C     ----- G(0,1) = D00 * G(0,0) -----
C
      DO 40 IG=1,NG
   40   GNM(IG,1,2)=D00(IG)*GNM(IG,1,1)
      IF(NMAXS) GO TO 60
C
C     ----- G(1,1) = B00 * G(0,0) + D00 * G(1,0) -----
C
      DO 50 IG=1,NG
   50   GNM(IG,2,2)=B00(IG,1)*GNM(IG,1,1)+D00(IG)*GNM(IG,2,1)
C
   60 MAX=MAX0(NMAX-1,MMAX-1)
      DO 70 M=2,MAX
        DO 70 IG=1,NG
   70     B00(IG,M)=B00(IG,M-1)+B00(IG,1)
C
      IF(NMAXP) GO TO 120
C
C     ----- G(N+1,0) = N * B10 * G(N-1,0) + C00 * G(N,0) -----
C
      DO 80 N=2,NMAX-1
        DO 80 IG=1,NG
   80     B10(IG,N)=B10(IG,N-1)+B10(IG,1)
      DO 90 N=2,NMAX-1
        DO 90 IG=1,NG
   90     GNM(IG,N+1,1)=B10(IG,N-1)*GNM(IG,N-1,1)+C00(IG)*GNM(IG,N,1)
      IF(MMAXS) GO TO 120
C
C     ----- G(N,1) = N * B00 * G(N-1,0) + D00 * G(N,0) -----
C
      DO 100 N=2,NMAX-1
        DO 100 IG=1,NG
  100     GNM(IG,N+1,2)=B00(IG,N)*GNM(IG,N,1)+D00(IG)*GNM(IG,N+1,1)
C
  120 IF(MMAXP) GO TO 170
C
C     ----- G(0,M+1) = M * B01 * G(0,M-1) + D00 * G(O,M) -----
C
      DO 130 M=2,MMAX-1
        DO 130 IG=1,NG
  130     B01(IG,M)=B01(IG,M-1)+B01(IG,1)
      DO 140 M=2,MMAX-1
        DO 140 IG=1,NG
  140     GNM(IG,1,M+1)=B01(IG,M-1)*GNM(IG,1,M-1)+D00(IG)*GNM(IG,1,M)
      IF(NMAXS) GO TO 170
C
C     ----- G(1,M) = M * B00 * G(0,M-1) + C00 * G(0,M) -----
C
      DO 150 M=2,MMAX-1
        DO 150 IG=1,NG
  150     GNM(IG,2,M+1)=B00(IG,M)*GNM(IG,1,M)+C00(IG)*GNM(IG,1,M+1)
C
  170 IF(NMAXP.OR.MMAXP) RETURN
C
C     ----- G(N+1,M) = N * B10 * G(N-1,M  ) -----
C                    +     C00 * G(N  ,M  )
C                    + M * B00 * G(N  ,M-1)
C
      DO 180 M=2,MMAX-1
        DO 180 N=2,NMAX-1
          DO 180 IG=1,NG
            GNM(IG,N+1,M+1)=B10(IG,N-1)*GNM(IG,N-1,M+1)+
     1                C00(IG    )*GNM(IG,N  ,M+1)+
     2                B00(IG,M  )*GNM(IG,N  ,M  )
  180 CONTINUE
C
      RETURN
      END
c--------------------------------------------------------------------------
C*MODULE GRD2C   *DECK JKWRYS
      SUBROUTINE xJKWRYS(RWV,ABV,NUMG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /xROOT  / XX,U(13),W(13),NROOTS
      DIMENSION RWV(2,NUMG,1),ABV(5,1)
C
      DO 595 NG=1,NUMG
         XX=ABV(5,NG)
         IF(NROOTS.LE.3) CALL xRT123
         IF(NROOTS.EQ.4) CALL xROOT4
         IF(NROOTS.EQ.5) CALL xROOT5
         IF(NROOTS.GE.6) CALL xROOT6
         DO 590 NR=1,NROOTS
           RWV(1,NG,NR)=U(NR)
           RWV(2,NG,NR)=W(NR)
 590     CONTINUE
 595  CONTINUE
C
      RETURN
      END
c--------------------------------------------------------------------------
C*MODULE GRD2C   *DECK JDXYZS
      SUBROUTINE xJDXYZS(GI,GIJ,GIJK,GIJKL,HIJK,HIJKL,FIJK,FIJKL,
     1 NG,NIMAX,NJMAX,NKMAX,NLMAX,NI,NJ,NK,NL,AAI,AAJ,AAK,AAL,
     2 FI,FJ,FK,FL,SII,SJJ,SKK,SLL,SIJ,SIK,SIL,SJK,SJL,SKL)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SKIPI,SKIPJ,SKIPK,SKIPL
      LOGICAL IS,JS,KS,LS
      COMMON/xDERSKP/IIAT,JJAT,KKAT,LLAT,SKIPI,SKIPJ,SKIPK,SKIPL
      DIMENSION GIJKL(NG,NLMAX,NKMAX*NJMAX*NIMAX)
      DIMENSION  GIJK(NG,NLMAX,NKMAX,NJMAX*NIMAX)
      DIMENSION   GIJ(NG,NLMAX*NKMAX,NJMAX,NIMAX)
      DIMENSION    GI(NG,NLMAX*NKMAX*NJMAX,NIMAX)
      DIMENSION HIJKL(NG,NLMAX,NKMAX*NJMAX,NIMAX)
      DIMENSION  HIJK(NG,NLMAX,NKMAX,NJMAX,NIMAX)
      DIMENSION FIJKL(NG,NLMAX,NKMAX,NJMAX,NIMAX)
      DIMENSION  FIJK(NG,NLMAX,NKMAX,NJMAX*NIMAX)
      DIMENSION   AAI(NG)
      DIMENSION   AAJ(NG)
      DIMENSION   AAK(NG)
      DIMENSION   AAL(NG)
      DIMENSION    FI(NG,NLMAX*NKMAX*NJMAX,NIMAX)
      DIMENSION    FJ(NG,NLMAX*NKMAX,NJMAX,NIMAX)
      DIMENSION    FK(NG,NLMAX,NKMAX,NJMAX*NIMAX)
      DIMENSION    FL(NG,NLMAX,NKMAX*NJMAX*NIMAX)
      DIMENSION   SII(NG,NLMAX*NKMAX*NJMAX,NIMAX)
      DIMENSION   SJJ(NG,NLMAX*NKMAX,NJMAX,NIMAX)
      DIMENSION   SKK(NG,NLMAX,NKMAX,NJMAX*NIMAX)
      DIMENSION   SLL(NG,NLMAX,NKMAX*NJMAX*NIMAX)
      DIMENSION   SIJ(NG,NLMAX*NKMAX,NJMAX,NIMAX)
      DIMENSION   SIK(NG,NLMAX,NKMAX,NJMAX,NIMAX)
      DIMENSION   SIL(NG,NLMAX,NKMAX*NJMAX,NIMAX)
      DIMENSION   SJK(NG,NLMAX,NKMAX,NJMAX,NIMAX)
      DIMENSION   SJL(NG,NLMAX,NKMAX,NJMAX,NIMAX)
      DIMENSION   SKL(NG,NLMAX,NKMAX,NJMAX*NIMAX)
C
      IS=NI.EQ.1
      JS=NJ.EQ.1
      KS=NK.EQ.1
      LS=NL.EQ.1
C
C     ----- FIRST DERIVATIVES ONLY -----
C
      IF(SKIPI) GO TO 1030
C
C     ----- -FI- ONLY -----
C
      DO 1010 LKJ=1,NLMAX*NKMAX*NJMAX
        DO 1010 IG =1,NG
           FI(IG,LKJ,1)=  GI(IG,LKJ,2)*AAI(IG)
 1010 CONTINUE
      IF(IS) GO TO 1030
      DO 1020 I  =2,NI
        DO 1020 LKJ=1,NLMAX*NKMAX*NJMAX
          DO 1020 IG =1,NG
             FI(IG,LKJ,I)= GI(IG,LKJ,I+1)*AAI(IG)
     1                    -GI(IG,LKJ,I-1)*(I-1)
 1020 CONTINUE
C
 1030 IF(SKIPJ) GO TO 1130
C
C     ----- -FJ- ONLY -----
C
      DO 1110 I =1,NIMAX
        DO 1110 LK=1,NLMAX*NKMAX
          DO 1110 IG=1,NG
             FJ(IG,LK,1,I)=  GIJ(IG,LK,2,I)*AAJ(IG)
 1110 CONTINUE
      IF(JS) GO TO 1130
      DO 1120 I =1,NIMAX
        DO 1120 J =2,NJ
          DO 1120 LK=1,NLMAX*NKMAX
            DO 1120 IG=1,NG
               FJ(IG,LK,J,I)= GIJ(IG,LK,J+1,I)*AAJ(IG)
     1                       -GIJ(IG,LK,J-1,I)*(J-1)
 1120 CONTINUE
C
 1130 IF(SKIPK) GO TO 1230
C
C     ----- -FK- ONLY -----
C
      DO 1210 JI=1,NJMAX*NIMAX
        DO 1210 L =1,NLMAX
          DO 1210 IG=1,NG
             FK(IG,L,1,JI)=  GIJK(IG,L,2,JI)*AAK(IG)
 1210 CONTINUE
      IF(KS) GO TO 1230
      DO 1220 JI=1,NJMAX*NIMAX
        DO 1220 K =2,NK
          DO 1220 L =1,NLMAX
            DO 1220 IG=1,NG
               FK(IG,L,K,JI)= GIJK(IG,L,K+1,JI)*AAK(IG)
     1                       -GIJK(IG,L,K-1,JI)*(K-1)
 1220 CONTINUE
C
 1230 IF(SKIPL) GO TO 1330
C
C     ----- -FL- AND -SLL- -----
C
      DO 1310 KJI=1,NKMAX*NJMAX*NIMAX
        DO 1310 IG =1,NG
           FL(IG,1,KJI)=  GIJKL(IG,2,KJI)*AAL(IG)
 1310 CONTINUE
      IF(LS) GO TO 1330
      DO 1320 KJI=1,NKMAX*NJMAX*NIMAX
        DO 1320 L  =2,NL
          DO 1320 IG =1,NG
             FL(IG,L,KJI)= GIJKL(IG,L+1,KJI)*AAL(IG)
     1                    -GIJKL(IG,L-1,KJI)*(L-1)
 1320 CONTINUE
C
 1330 CONTINUE
      RETURN
      END
c--------------------------------------------------------------------------
c--------------------------------------------------------------------------
C*MODULE GRD2C   *DECK JKDZER
      SUBROUTINE xJKDZER(GIJKL,FIJKL,SIJKL,IJKLG)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OUT,DBG
C
      COMMON /xDSHLT / RTOL,DTOL,VTOL1,VTOL2,VTOLS,OUT,DBG
      COMMON/xDERPAR/IDER,JDER,KDER,LDER,NDER
      COMMON/xDERSHL/FD(12),SD(78)
      COMMON/xSHLNUM/NUMI,NUMJ,NUMK,NUML,IJKL
      DIMENSION GIJKL(*),FIJKL(12,*),SIJKL(78,*),IJKLG(4,*)
      PARAMETER (ZERO=0.0D+00)
cwc1
      common /gms_deriv_ctrl/ want_integrals, all_centers
      logical want_integrals, all_centers
cwc0
C
C     ----- ZERO OUT -GIJKL- -----
C
cwc      IF (OUT) THEN
      IF (OUT.or.want_integrals) THEN
        DO 10 N=1,IJKL
 10     GIJKL(IJKLG(1,N))=ZERO
      ENDIF
      IF(NDER.EQ.0) RETURN
C
C     ----- ZERO OUT -FIJKL- -----
C
cwc      IF(OUT.OR.NDER.EQ.2) THEN
      IF(OUT.OR.NDER.EQ.2.or.want_integrals) THEN
        DO 30 N=1,IJKL
          DO 20 M=1,12
 20       FIJKL(M,IJKLG(1,N))=ZERO
 30     CONTINUE
      ENDIF
C
      DO 40 I=1,12
   40 FD(I)=ZERO
      IF(NDER.EQ.1) RETURN
C
C     ----- ZERO OUT -SIJKL- -----
C
cwc      IF (OUT) THEN
      IF (OUT.or.want_integrals) THEN
        DO 60 N=1,IJKL
          DO 50 M=1,78
 50       SIJKL(M,IJKLG(1,N))=ZERO
 60     CONTINUE
      ENDIF
C
      DO 70 I=1,78
   70 SD(I)=ZERO
      RETURN
      END
c--------------------------------------------------------------------------
C*MODULE GRD2A   *DECK DSPDFS
      SUBROUTINE xDSPDFS(NG,NR,IJKLG,GIJKL,FIJKL,SIJKL,XYZ,
     1 FIXYZ,FJXYZ,FKXYZ,FLXYZ,SIIXYZ,SJJXYZ,SKKXYZ,SLLXYZ,
     2 SIJXYZ,SIKXYZ,SILXYZ,SJKXYZ,SJLXYZ,SKLXYZ,
     3 XY,XZ,YZ,X,Y,Z,SJ,SK,SL,DIJSI,DIJSJ,DKLSK,DKLSL,DAB)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL IS,JS,KS,LS
      LOGICAL IJS,IJKS,IJKLS
      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL, use_symm
      LOGICAL SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      LOGICAL SKIPI,SKIPJ,SKIPK,SKIPL
      LOGICAL OUT,DBG
      COMMON /xDSHLNO/ LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     *                NIJ,IJD,KLD,IJ,KL
      COMMON /xDSHLT / RTOL,DTOL,VTOL1,VTOL2,VTOLS,OUT,DBG
      COMMON/xDERSKP/IIAT,JJAT,KKAT,LLAT,SKIPI,SKIPJ,SKIPK,SKIPL
      COMMON/xDERPAR/IDER,JDER,KDER,LDER,NDER
      COMMON/xDERSHL/FD(12),SD(78)
      COMMON/xSHLEQU/IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL, use_symm
      COMMON/xSHLTYP/SPI,SPJ,SPK,SPL,SPIJ,SPKL,SPIJKL
      DIMENSION GIJKL(1),FIJKL(12,1),SIJKL(78,1),IJKLG(4,1)
      DIMENSION DAB(1)
      DIMENSION    XYZ(NG*NR,1)
      DIMENSION  FIXYZ(NG*NR,1), FJXYZ(NG*NR,1)
      DIMENSION  FKXYZ(NG*NR,1), FLXYZ(NG*NR,1)
      DIMENSION SIIXYZ(NG*NR,1),SJJXYZ(NG*NR,1)
      DIMENSION SKKXYZ(NG*NR,1),SLLXYZ(NG*NR,1)
      DIMENSION SIJXYZ(NG*NR,1),SIKXYZ(NG*NR,1)
      DIMENSION SILXYZ(NG*NR,1),SJKXYZ(NG*NR,1)
      DIMENSION SJLXYZ(NG*NR,1),SKLXYZ(NG*NR,1)
      DIMENSION XY(NG*NR),XZ(NG*NR),YZ(NG*NR)
      DIMENSION  X(NG*NR), Y(NG*NR), Z(NG*NR)
      DIMENSION SJ(NG*NR),SK(NG*NR),SL(NG*NR)
      DIMENSION DIJSI(NG*NR),DIJSJ(NG*NR),DKLSK(NG*NR),DKLSL(NG*NR)
      PARAMETER (ZERO=0.0D+00)
cwc1
      common /gms_deriv_ctrl/ want_integrals, all_centers
      logical want_integrals, all_centers
cwc0
C
      IF(SPIJKL) GO TO 1000
C
C     ----- NO SHARED EXPONENTS ; SUM UP ( IX * IY * IZ ) -----
C
      IF(NDER.LT.1) GO TO 700
C
C     ----- GRADIENT -----
C
      IJKLN=0
      DO 640 I=MINI,MAXI
        JMAX=MAXJ
        IF(use_symm .and. IIEQJJ) JMAX=I
        DO 630 J=MINJ,JMAX
          KMAX=MAXK
          IF(use_symm .and. IJEQKL) KMAX=I
          DO 620 K=MINK,KMAX
            LMAX=MAXL
            if (use_symm) then
               IF(KKEQLL           ) LMAX=K
               IF(IJEQKL.AND.K.EQ.I) LMAX=J
            endif

            DO 610 L=MINL,LMAX
              IJKLN=IJKLN+1
              NN=IJKLG(1,IJKLN)
              NX=IJKLG(2,IJKLN)
              NY=IJKLG(3,IJKLN)
              NZ=IJKLG(4,IJKLN)
C
              DO 510 IGR=1,NG*NR
                XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)
                XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)
                YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)
  510         CONTINUE
      
C
              IF(SKIPI) GO TO 530
              DUMFX =ZERO
              DUMFY =ZERO
              DUMFZ =ZERO
              DO 520 IGR=1,NG*NR
                DUMFX =DUMFX + FIXYZ(IGR,NX)*YZ(IGR)
                DUMFY =DUMFY + FIXYZ(IGR,NY)*XZ(IGR)
                DUMFZ =DUMFZ + FIXYZ(IGR,NZ)*XY(IGR)
  520         CONTINUE
              FD( 1)=FD( 1)+DAB(NN)*DUMFX
              FD( 2)=FD( 2)+DAB(NN)*DUMFY
              FD( 3)=FD( 3)+DAB(NN)*DUMFZ
  530         IF(SKIPJ) GO TO 550
              DUMFX =ZERO
              DUMFY =ZERO
              DUMFZ =ZERO
              DO 540 IGR=1,NG*NR
                DUMFX =DUMFX + FJXYZ(IGR,NX)*YZ(IGR)
                DUMFY =DUMFY + FJXYZ(IGR,NY)*XZ(IGR)
  540           DUMFZ =DUMFZ + FJXYZ(IGR,NZ)*XY(IGR)
              FD( 4)=FD( 4)+DAB(NN)*DUMFX
              FD( 5)=FD( 5)+DAB(NN)*DUMFY
              FD( 6)=FD( 6)+DAB(NN)*DUMFZ
  550         IF(SKIPK) GO TO 570
              DUMFX =ZERO
              DUMFY =ZERO
              DUMFZ =ZERO
              DO 560 IGR=1,NG*NR
                DUMFX =DUMFX + FKXYZ(IGR,NX)*YZ(IGR)
                DUMFY =DUMFY + FKXYZ(IGR,NY)*XZ(IGR)
  560           DUMFZ =DUMFZ + FKXYZ(IGR,NZ)*XY(IGR)
              FD( 7)=FD( 7)+DAB(NN)*DUMFX
              FD( 8)=FD( 8)+DAB(NN)*DUMFY
              FD( 9)=FD( 9)+DAB(NN)*DUMFZ
  570         IF(SKIPL) GO TO 600
              DUMFX =ZERO
              DUMFY =ZERO
              DUMFZ =ZERO
              DO 580 IGR=1,NG*NR
                DUMFX =DUMFX + FLXYZ(IGR,NX)*YZ(IGR)
                DUMFY =DUMFY + FLXYZ(IGR,NY)*XZ(IGR)
  580           DUMFZ =DUMFZ + FLXYZ(IGR,NZ)*XY(IGR)
              FD(10)=FD(10)+DAB(NN)*DUMFX
              FD(11)=FD(11)+DAB(NN)*DUMFY
              FD(12)=FD(12)+DAB(NN)*DUMFZ
  600         CONTINUE
C
  610       CONTINUE
  620     CONTINUE
  630   CONTINUE
  640 CONTINUE
C
  700 CONTINUE
C
C     ----- ZEROTH AND FIRST DERIVATIVE INTEGRALS -----
C
cwc      IF(.NOT.OUT) RETURN
      IF(.NOT.OUT.and.(.not.want_integrals)) RETURN
C
      IJKLN=0
      DO 940 I=MINI,MAXI
        JMAX=MAXJ
        IF(use_symm .and. IIEQJJ) JMAX=I
        DO 930 J=MINJ,JMAX
          KMAX=MAXK
          IF(use_symm .and. IJEQKL) KMAX=I
          DO 920 K=MINK,KMAX
            LMAX=MAXL
            if (use_symm) then
               IF(KKEQLL           ) LMAX=K
               IF(IJEQKL.AND.K.EQ.I) LMAX=J
            endif

            DO 910 L=MINL,LMAX
            IJKLN=IJKLN+1
            NN=IJKLG(1,IJKLN)
            NX=IJKLG(2,IJKLN)
            NY=IJKLG(3,IJKLN)
            NZ=IJKLG(4,IJKLN)
C
            DO 710 IGR=1,NG*NR
              XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)
              XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)
  710         YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)
            DUM=ZERO
            DO 720 IGR=1,NG*NR
  720         DUM=DUM+XYZ(IGR,NX)*YZ(IGR)
            GIJKL(NN)=GIJKL(NN)+DUM
C
            IF(SKIPI) GO TO 830
            DUMFX =ZERO
            DUMFY =ZERO
            DUMFZ =ZERO
            DO 820 IGR=1,NG*NR
              DUMFX =DUMFX + FIXYZ(IGR,NX)*YZ(IGR)
              DUMFY =DUMFY + FIXYZ(IGR,NY)*XZ(IGR)
  820         DUMFZ =DUMFZ + FIXYZ(IGR,NZ)*XY(IGR)
            FIJKL( 1,NN)=FIJKL( 1,NN)+DUMFX
            FIJKL( 2,NN)=FIJKL( 2,NN)+DUMFY
            FIJKL( 3,NN)=FIJKL( 3,NN)+DUMFZ

  830       IF(SKIPJ) GO TO 850
            DUMFX =ZERO
            DUMFY =ZERO
            DUMFZ =ZERO
            DO 840 IGR=1,NG*NR
              DUMFX =DUMFX + FJXYZ(IGR,NX)*YZ(IGR)
              DUMFY =DUMFY + FJXYZ(IGR,NY)*XZ(IGR)
  840         DUMFZ =DUMFZ + FJXYZ(IGR,NZ)*XY(IGR)
            FIJKL( 4,NN)=FIJKL( 4,NN)+DUMFX
            FIJKL( 5,NN)=FIJKL( 5,NN)+DUMFY
            FIJKL( 6,NN)=FIJKL( 6,NN)+DUMFZ
  850       IF(SKIPK) GO TO 870
            DUMFX =ZERO
            DUMFY =ZERO
            DUMFZ =ZERO
            DO 860 IGR=1,NG*NR
              DUMFX =DUMFX + FKXYZ(IGR,NX)*YZ(IGR)
              DUMFY =DUMFY + FKXYZ(IGR,NY)*XZ(IGR)
  860         DUMFZ =DUMFZ + FKXYZ(IGR,NZ)*XY(IGR)
            FIJKL( 7,NN)=FIJKL( 7,NN)+DUMFX
            FIJKL( 8,NN)=FIJKL( 8,NN)+DUMFY
            FIJKL( 9,NN)=FIJKL( 9,NN)+DUMFZ
  870       IF(SKIPL) GO TO 900
            DUMFX =ZERO
            DUMFY =ZERO
            DUMFZ =ZERO
            DO 880 IGR=1,NG*NR
              DUMFX =DUMFX + FLXYZ(IGR,NX)*YZ(IGR)
              DUMFY =DUMFY + FLXYZ(IGR,NY)*XZ(IGR)
  880         DUMFZ =DUMFZ + FLXYZ(IGR,NZ)*XY(IGR)
            FIJKL(10,NN)=FIJKL(10,NN)+DUMFX
            FIJKL(11,NN)=FIJKL(11,NN)+DUMFY
            FIJKL(12,NN)=FIJKL(12,NN)+DUMFZ
  900       CONTINUE
C
  910       CONTINUE
  920     CONTINUE
  930   CONTINUE
  940 CONTINUE
C
      RETURN
C
C     ----- SHARED EXPONENTS ; FORM ( IX * IY * IZ ) -----
C
 1000 CONTINUE
      IF(NDER.LT.1) GO TO 3000
C
C     ----- GRADIENT -----
C
      IJKLN=0
      DO 2640 I=MINI,MAXI
        IS=SPI.AND.I.EQ.1
C
        JMAX=MAXJ
        IF(use_symm .and. IIEQJJ) JMAX=I
        DO 2630 J=MINJ,JMAX
          JS=SPJ.AND.J.EQ.1
C
          IF(JS) THEN
             IF(IS) THEN
                DO 2110 IGR=1,NG*NR
 2110             SJ(IGR)=DIJSJ(IGR)*DIJSI(IGR)
             ELSE
                DO 2120 IGR=1,NG*NR
 2120             SJ(IGR)=DIJSJ(IGR)
             ENDIF
          ELSE
             IF(IS) THEN
                DO 2130 IGR=1,NG*NR
 2130             SJ(IGR)=DIJSI(IGR)
             ENDIF
          ENDIF
          IJS=IS.OR.JS
C
          KMAX=MAXK
          IF(use_symm .and. IJEQKL) KMAX=I
          DO 2620 K=MINK,KMAX
            KS=SPK.AND.K.EQ.1
C
            IF(KS) THEN
               IF(IJS) THEN
                  DO 2210 IGR=1,NG*NR
 2210               SK(IGR)=DKLSK(IGR)*SJ(IGR)
               ELSE
                  DO 2220 IGR=1,NG*NR
 2220               SK(IGR)=DKLSK(IGR)
               ENDIF
            ELSE
               IF(IJS) THEN
                  DO 2230 IGR=1,NG*NR
 2230               SK(IGR)=SJ(IGR)
               ENDIF
            ENDIF
            IJKS=IJS.OR.KS
C
            LMAX=MAXL
            if (use_symm) then
               IF(KKEQLL           ) LMAX=K
               IF(IJEQKL.AND.K.EQ.I) LMAX=J
            endif

            DO 2610 L=MINL,LMAX
              LS=SPL.AND.L.EQ.1
C
              IF(LS) THEN
                 IF(IJKS) THEN
                    DO 2310 IGR=1,NG*NR
 2310                 SL(IGR)=DKLSL(IGR)*SK(IGR)
                 ELSE
                    DO 2320 IGR=1,NG*NR
 2320                 SL(IGR)=DKLSL(IGR)
                 ENDIF
              ELSE
                 IF(IJKS) THEN
                    DO 2330 IGR=1,NG*NR
 2330                 SL(IGR)=SK(IGR)
                 ENDIF
              ENDIF
              IJKLS=IJKS.OR.LS
C
              IJKLN=IJKLN+1
              NN=IJKLG(1,IJKLN)
              NX=IJKLG(2,IJKLN)
              NY=IJKLG(3,IJKLN)
              NZ=IJKLG(4,IJKLN)
C
              IF(IJKLS) THEN
                 DO 2500 IGR=1,NG*NR
                   XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)*SL(IGR)
                   XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)*SL(IGR)
 2500              YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)*SL(IGR)
              ELSE
                 DO 2510 IGR=1,NG*NR
                   XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)
                   XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)
 2510              YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)
              ENDIF
C
              IF(SKIPI) GO TO 2530
              DUMFX =ZERO
              DUMFY =ZERO
              DUMFZ =ZERO
              DO 2520 IGR=1,NG*NR
                DUMFX =DUMFX + FIXYZ(IGR,NX)*YZ(IGR)
                DUMFY =DUMFY + FIXYZ(IGR,NY)*XZ(IGR)
 2520           DUMFZ =DUMFZ + FIXYZ(IGR,NZ)*XY(IGR)
              FD( 1)=FD( 1)+DAB(NN)*DUMFX
              FD( 2)=FD( 2)+DAB(NN)*DUMFY
              FD( 3)=FD( 3)+DAB(NN)*DUMFZ
 2530         IF(SKIPJ) GO TO 2550
              DUMFX =ZERO
              DUMFY =ZERO
              DUMFZ =ZERO
              DO 2540 IGR=1,NG*NR
                DUMFX =DUMFX + FJXYZ(IGR,NX)*YZ(IGR)
                DUMFY =DUMFY + FJXYZ(IGR,NY)*XZ(IGR)
 2540           DUMFZ =DUMFZ + FJXYZ(IGR,NZ)*XY(IGR)
              FD( 4)=FD( 4)+DAB(NN)*DUMFX
              FD( 5)=FD( 5)+DAB(NN)*DUMFY
              FD( 6)=FD( 6)+DAB(NN)*DUMFZ
 2550         IF(SKIPK) GO TO 2570
              DUMFX =ZERO
              DUMFY =ZERO
              DUMFZ =ZERO
              DO 2560 IGR=1,NG*NR
                DUMFX =DUMFX + FKXYZ(IGR,NX)*YZ(IGR)
                DUMFY =DUMFY + FKXYZ(IGR,NY)*XZ(IGR)
 2560           DUMFZ =DUMFZ + FKXYZ(IGR,NZ)*XY(IGR)
              FD( 7)=FD( 7)+DAB(NN)*DUMFX
              FD( 8)=FD( 8)+DAB(NN)*DUMFY
              FD( 9)=FD( 9)+DAB(NN)*DUMFZ
 2570         IF(SKIPL) GO TO 2600
              DUMFX =ZERO
              DUMFY =ZERO
              DUMFZ =ZERO
              DO 2580 IGR=1,NG*NR
                DUMFX =DUMFX + FLXYZ(IGR,NX)*YZ(IGR)
                DUMFY =DUMFY + FLXYZ(IGR,NY)*XZ(IGR)
 2580           DUMFZ =DUMFZ + FLXYZ(IGR,NZ)*XY(IGR)
              FD(10)=FD(10)+DAB(NN)*DUMFX
              FD(11)=FD(11)+DAB(NN)*DUMFY
              FD(12)=FD(12)+DAB(NN)*DUMFZ
 2600         CONTINUE
C
 2610       CONTINUE
 2620     CONTINUE
 2630   CONTINUE
 2640 CONTINUE
C
 3000 CONTINUE
C
C     ----- ZEROTH AND FIRST DERIVATIVE INTEGRALS -----
C
cwc      IF(.NOT.OUT) RETURN
      IF(.NOT.OUT.and.(.not.want_integrals)) RETURN
C
      IJKLN=0
      DO 3640 I=MINI,MAXI
        IS=SPI.AND.I.EQ.1
C
        JMAX=MAXJ
        IF(use_symm .and. IIEQJJ) JMAX=I
        DO 3630 J=MINJ,JMAX
          JS=SPJ.AND.J.EQ.1
C
          IF(JS) THEN
             IF(IS) THEN
                DO 3110 IGR=1,NG*NR
 3110             SJ(IGR)=DIJSJ(IGR)*DIJSI(IGR)
             ELSE
                DO 3120 IGR=1,NG*NR
 3120             SJ(IGR)=DIJSJ(IGR)
             ENDIF
          ELSE
             IF(IS) THEN
                DO 3130 IGR=1,NG*NR
 3130             SJ(IGR)=DIJSI(IGR)
             ENDIF
          ENDIF
          IJS=IS.OR.JS
C
          KMAX=MAXK
          IF(use_symm .and. IJEQKL) KMAX=I
          DO 3620 K=MINK,KMAX
            KS=SPK.AND.K.EQ.1
C
            IF(KS) THEN
               IF(IJS) THEN
                  DO 3210 IGR=1,NG*NR
 3210               SK(IGR)=DKLSK(IGR)*SJ(IGR)
               ELSE
                  DO 3220 IGR=1,NG*NR
 3220               SK(IGR)=DKLSK(IGR)
               ENDIF
            ELSE
               IF(IJS) THEN
                  DO 3230 IGR=1,NG*NR
 3230               SK(IGR)=SJ(IGR)
               ENDIF
            ENDIF
            IJKS=IJS.OR.KS
C
            LMAX=MAXL
            if (use_symm) then
              IF(KKEQLL           ) LMAX=K
              IF(IJEQKL.AND.K.EQ.I) LMAX=J
            endif

            DO 3610 L=MINL,LMAX
              LS=SPL.AND.L.EQ.1
C
              IF(LS) THEN
                 IF(IJKS) THEN
                    DO 3310 IGR=1,NG*NR
 3310                 SL(IGR)=DKLSL(IGR)*SK(IGR)
                 ELSE
                    DO 3320 IGR=1,NG*NR
 3320                 SL(IGR)=DKLSL(IGR)
                 ENDIF
              ELSE
                 IF(IJKS) THEN
                    DO 3330 IGR=1,NG*NR
 3330                 SL(IGR)=SK(IGR)
                 ENDIF
              ENDIF
              IJKLS=IJKS.OR.LS
C
              IJKLN=IJKLN+1
              NN=IJKLG(1,IJKLN)
              NX=IJKLG(2,IJKLN)
              NY=IJKLG(3,IJKLN)
              NZ=IJKLG(4,IJKLN)
C
              IF(IJKLS) THEN
                 DO 3400 IGR=1,NG*NR
                   XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)*SL(IGR)
                   XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)*SL(IGR)
 3400              YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)*SL(IGR)
              ELSE
                 DO 3410 IGR=1,NG*NR
                   XY(IGR)=XYZ(IGR,NX)*XYZ(IGR,NY)
                   XZ(IGR)=XYZ(IGR,NX)*XYZ(IGR,NZ)
 3410              YZ(IGR)=XYZ(IGR,NY)*XYZ(IGR,NZ)
              ENDIF
C
              DUM=ZERO
              DO 3420 IGR=1,NG*NR
 3420           DUM=DUM+XYZ(IGR,NX)*YZ(IGR)
              GIJKL(NN)=GIJKL(NN)+DUM
C
              IF(SKIPI) GO TO 3530
              DUMFX =ZERO
              DUMFY =ZERO
              DUMFZ =ZERO
              DO 3520 IGR=1,NG*NR
                DUMFX =DUMFX + FIXYZ(IGR,NX)*YZ(IGR)
                DUMFY =DUMFY + FIXYZ(IGR,NY)*XZ(IGR)
 3520           DUMFZ =DUMFZ + FIXYZ(IGR,NZ)*XY(IGR)
              FIJKL( 1,NN)=FIJKL( 1,NN)+DUMFX
              FIJKL( 2,NN)=FIJKL( 2,NN)+DUMFY
              FIJKL( 3,NN)=FIJKL( 3,NN)+DUMFZ
 3530         IF(SKIPJ) GO TO 3550
              DUMFX =ZERO
              DUMFY =ZERO
              DUMFZ =ZERO
              DO 3540 IGR=1,NG*NR
                DUMFX =DUMFX + FJXYZ(IGR,NX)*YZ(IGR)
                DUMFY =DUMFY + FJXYZ(IGR,NY)*XZ(IGR)
 3540           DUMFZ =DUMFZ + FJXYZ(IGR,NZ)*XY(IGR)
              FIJKL( 4,NN)=FIJKL( 4,NN)+DUMFX
              FIJKL( 5,NN)=FIJKL( 5,NN)+DUMFY
              FIJKL( 6,NN)=FIJKL( 6,NN)+DUMFZ
 3550         IF(SKIPK) GO TO 3570
              DUMFX =ZERO
              DUMFY =ZERO
              DUMFZ =ZERO
              DO 3560 IGR=1,NG*NR
                DUMFX =DUMFX + FKXYZ(IGR,NX)*YZ(IGR)
                DUMFY =DUMFY + FKXYZ(IGR,NY)*XZ(IGR)
 3560           DUMFZ =DUMFZ + FKXYZ(IGR,NZ)*XY(IGR)
              FIJKL( 7,NN)=FIJKL( 7,NN)+DUMFX
              FIJKL( 8,NN)=FIJKL( 8,NN)+DUMFY
              FIJKL( 9,NN)=FIJKL( 9,NN)+DUMFZ
 3570         IF(SKIPL) GO TO 3600
              DUMFX =ZERO
              DUMFY =ZERO
              DUMFZ =ZERO
              DO 3580 IGR=1,NG*NR
                DUMFX =DUMFX + FLXYZ(IGR,NX)*YZ(IGR)
                DUMFY =DUMFY + FLXYZ(IGR,NY)*XZ(IGR)
 3580           DUMFZ =DUMFZ + FLXYZ(IGR,NZ)*XY(IGR)
              FIJKL(10,NN)=FIJKL(10,NN)+DUMFX
              FIJKL(11,NN)=FIJKL(11,NN)+DUMFY
              FIJKL(12,NN)=FIJKL(12,NN)+DUMFZ
 3600         CONTINUE
C
 3610       CONTINUE
 3620     CONTINUE
 3630   CONTINUE
 3640 CONTINUE
C
      RETURN
      END
c--------------------------------------------------------------------------
C*MODULE GRD2A   *DECK JKDINV
      SUBROUTINE xJKDINV(FIJKL,SIJKL,IJKLG,FC,NC,DF,LDF,NBF,DEN)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "gamess.par"
C
      LOGICAL OUT,DBG
      LOGICAL SKIPI,SKIPJ,SKIPK,SKIPL
      COMMON /xDSHLT / RTOL,DTOL,VTOL1,VTOL2,VTOLS,OUT,DBG
      COMMON /xGRAD  / DE(3,MXATM)
      COMMON /xIOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON/xDERINV/INVTYP
      COMMON/xDERSHL/FD(3,4),SD(78)
      COMMON/xDERPAR/IDER,JDER,KDER,LDER,NDER
      COMMON/xDERSKP/IIAT,JJAT,KKAT,LLAT,SKIPI,SKIPJ,SKIPK,SKIPL
      DIMENSION FIJKL(3,4,1),SIJKL(78,1),IJKLG(4,1)
      DIMENSION FC(NC,1),DF(LDF,1),DEN(1)
      DIMENSION KAT(4)
      EQUIVALENCE (KAT(1),IIAT)
C
      IF(NDER.EQ.0) RETURN
C
C     ----- TRANSLATIONAL INVARIANCE FOR GRADIENT ELEMENTS -----
C
      IF (INVTYP.EQ.2) THEN
        DO 210 IXYZ=1,3
  210     FD(IXYZ,1)=- FD(IXYZ,4)
      ELSE IF (INVTYP.EQ.3) THEN
        DO 310 IXYZ=1,3
  310     FD(IXYZ,1)=- FD(IXYZ,3)
      ELSE IF (INVTYP.EQ.4.OR.INVTYP.EQ.5) THEN
        DO 410 IXYZ=1,3
  410     FD(IXYZ,1)=-(FD(IXYZ,3)+FD(IXYZ,4))
      ELSE IF (INVTYP.EQ.6) THEN
        DO 610 IXYZ=1,3
  610     FD(IXYZ,1)=- FD(IXYZ,2)
      ELSE IF (INVTYP.EQ.7) THEN
        DO 710 IXYZ=1,3
  710     FD(IXYZ,1)=-(FD(IXYZ,2)+FD(IXYZ,4))
      ELSE IF (INVTYP.EQ.8) THEN
        DO 810 IXYZ=1,3
  810     FD(IXYZ,1)=-(FD(IXYZ,2)+FD(IXYZ,4))
      ELSE IF (INVTYP.EQ.9.OR.INVTYP.EQ.10) THEN
        DO 910 IXYZ=1,3
  910     FD(IXYZ,1)=-(FD(IXYZ,2)+FD(IXYZ,3))
      ELSE IF (INVTYP.EQ.11) THEN
        DO 1110 IXYZ=1,3
 1110     FD(IXYZ,2)=- FD(IXYZ,1)
      ELSE IF (INVTYP.EQ.12) THEN
        DO 1210 IXYZ=1,3
 1210     FD(IXYZ,2)=-(FD(IXYZ,1)+FD(IXYZ,4))
      ELSE IF (INVTYP.EQ.13) THEN
        DO 1310 IXYZ=1,3
 1310     FD(IXYZ,2)=-(FD(IXYZ,1)+FD(IXYZ,3))
      ELSE IF (INVTYP.EQ.14) THEN
        DO 1410 IXYZ=1,3
 1410     FD(IXYZ,3)=-(FD(IXYZ,1)+FD(IXYZ,2))
      ELSE IF (INVTYP.EQ.15) THEN
        DO 1510 IXYZ=1,3
 1510     FD(IXYZ,4)=-(FD(IXYZ,1)+FD(IXYZ,2)+FD(IXYZ,3))
      ENDIF
C
 2000 CONTINUE
C
c      write (96,'(a,i2,4g15.6)') 'jkdinv:',invtyp,(fd(i,3),i=1,4)
      DO 2010 I=1,4
        IAT=KAT(I)
        DO 2010 IXYZ=1,3
 2010     DE(IXYZ,IAT)=DE(IXYZ,IAT)+FD(IXYZ,I)
c      write (96,'(4i3)') (kat(i),i=1,4)    
c      write (96,'(a,2x,4g15.6)') '       ',(de(3,kat(i)),i=1,4)
C
      IF(OUT) WRITE(IW,9999)(I,IIAT,DE(I,IIAT),I,JJAT,DE(I,JJAT),
     1                       I,KKAT,DE(I,KKAT),I,LLAT,DE(I,LLAT),I=1,3)
C
      RETURN
 9999 FORMAT(4(' DE(',I2,',',I3,') = ',E18.11))
      END
c--------------------------------------------------------------------------
C*MODULE GRD2C   *DECK JKDOUT
      SUBROUTINE xJKDOUT(GIJKL,FIJKL,SIJKL,IJKLG,DAB)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "gamess.par"
C
      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL, use_symm
      LOGICAL BUG1,BUG2,BUG
C
      COMMON /xDSHLNO/ LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     *                NIJ,IJD,KLD,IJ,KL
      COMMON /xIJPAIR/ IA(MXAO)
      COMMON /xIOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON/xDERPAR/IDER,JDER,KDER,LDER,NDER
      COMMON/xSHLEQU/IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL, use_symm
      COMMON/xSHLNRM/PNRM(35)
C
      DIMENSION GIJKL(1),FIJKL(12,1),SIJKL(78,1),IJKLG(4,1)
      DIMENSION DAB(1)
      DIMENSION F(12),S(12,12)
C
      PARAMETER (ZERO=0.0D+00, TOL=1.0D-10)
C
      BUG=.FALSE.
C
      IJKLN=0
      DO 400 I=MINI,MAXI
        PI=PNRM(I)
        JMAX=MAXJ
        IF(use_symm .and. IIEQJJ) JMAX=I
        DO 300 J=MINJ,JMAX
          PIJ=PNRM(J)*PI
          KMAX=MAXK
          IF(use_symm .and. IJEQKL) KMAX=I
          DO 200 K=MINK,KMAX
            PIJK=PNRM(K)*PIJ
            LMAX=MAXL
            if (use_symm) then
               IF(KKEQLL           ) LMAX=K
               IF(IJEQKL.AND.K.EQ.I) LMAX=J
            endif

            DO 100 L=MINL,LMAX
              PIJKL=PNRM(L)*PIJK
C
              IJKLN=IJKLN+1
              INUM=LOCI+I
              JNUM=LOCJ+J
              KNUM=LOCK+K
              LNUM=LOCL+L
C
              VAL=GIJKL(IJKLG(1,IJKLN))*PIJKL
              DEN=  DAB(IJKLG(1,IJKLN))
C
              IF(NDER.EQ.0) GO TO 50
C
              DO 10 N=1,12
   10           F(N)=FIJKL( N,IJKLG(1,IJKLN))*PIJKL
              ERR1=ZERO
              DO 20 N=1,12
   20           ERR1=ERR1+F(N)
              BUG1= ABS(ERR1).GT.TOL
              IF(BUG1) THEN
                BUG=.TRUE.
                WRITE(IW,9999) INUM,JNUM,KNUM,LNUM,VAL,DEN,
     1             IJKLG(1,IJKLN),IJKLG(2,IJKLN),IJKLG(3,IJKLN),
     2             IJKLG(4,IJKLN)
                WRITE(IW,9998) ERR1
                WRITE(IW,9995) (F(N),N=1,12)
              ENDIF
C
              IF(NDER.GE.1) GO TO 50
C
              DO 30 M=1,12
                DO 30 N=1,M
                  S(M,N)=SIJKL(IA(M)+N,IJKLG(1,IJKLN))*PIJKL
   30             S(N,M)=SIJKL(IA(M)+N,IJKLG(1,IJKLN))*PIJKL
              ERR2=ZERO
              DO 40 N=1,12
                DO 40 M=1,12
   40             ERR2=ERR2+S(M,N)
              BUG2= ABS(ERR2).GT.TOL
              IF(BUG2) THEN
                BUG=.TRUE.
                WRITE(IW,9999) INUM,JNUM,KNUM,LNUM,VAL,DEN,
     1             IJKLG(1,IJKLN),IJKLG(2,IJKLN),IJKLG(3,IJKLN),
     2             IJKLG(4,IJKLN)
                WRITE(IW,9997) ERR2
                WRITE(IW,9994) ((S(M,N),M=1,12),N=1,12)
              ENDIF
C
   50         CONTINUE
C
  100       CONTINUE
  200     CONTINUE
  300   CONTINUE
  400 CONTINUE
C
      IF(.NOT.BUG) WRITE(IW,9996)
C
      RETURN
 9999 FORMAT(' I,J,K,L = ',4I3,2E20.12,4I5)
 9998 FORMAT(' ERROR IN 1-DERIVATIVE INTEGRALS , ERR1 = ',F20.12)
 9997 FORMAT(' ERROR IN 2-DERIVATIVE INTEGRALS , ERR2 = ',F20.12)
 9996 FORMAT(' NO -TRANSLATION INVARIANCE ERROR DETECTED',
     1       ' FOR THIS GROUP OF SHELLS. GOOD GOING')
 9995 FORMAT(3E20.12)
 9994 FORMAT(6E15.8)
      END
c--------------------------------------------------------------------------
c------------------ Pople routines ----------------------------------------
c--------------------------------------------------------------------------
C*MODULE GRD2B   *DECK JKDG80
      SUBROUTINE xJKDG80(DABMAX,INEW,JNEW,KNEW,LNEW,
     2                  JTYPE,IAT,JAT,KAT,LAT)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ---- JKDG80 DRIVER FOR 2-E INTEGRAL DERIVATIVES FOR SP BASES ----
C     FROM GAUSSIAN-80, MICHEL DUPUIS, 6/6/85.
C     FROM HONDO8, BRETT BODE, 7/15/94.
C     REWRITTEN USING LOOPS AND INLINING FPPPP, JOSE SIERRA, 7/2001
C
#include "gamess.par"
C
      LOGICAL OUT,DBG
      COMMON /xDSHLT / RTOL,DTOL,VTOL1,VTOL2,VTOLS,OUT,DBG
      COMMON /xGRAD  / FXYZ(3,MXATM)
      COMMON /xINFOA / NAT,ICH,MUL,NUM,NX,NE,NA,NB,ZAN(MXATM),C(3,MXATM)
      COMMON /xIOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     2                CF(MXGTOT),CG(MXGTOT),
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
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
CJMS  LABELLED COMMON GSPG80 DEFINED FOR COMPUTATIONAL EFFICIENCY.
CJMS  FOR SP BASES ONLY, IT CONTAINS THE E ARRAY WHICH IS THE DAB
CJMS  ARRAY WITH INDICES IN REVERSE ORDER: E(I,J,K,L)= DAB(L,K,J,I)
CJMS  AND IS USED IN SUB JKDG80 (MODULE GRD2B). IT ORIGINATES IN:
CJMS
CJMS     1. SUBS DABCLU, DABDFT, DABGVB, DABMC, DABMP2 AND DABUMP
CJMS        (MODULE GRD2A) AND SUB DABPAU (MODULE EFPAUL) WHICH ARE
CJMS        ALL CALLED BY SUB JKDER (MODULE GRD2A)
CJMS
CJMS     2. SUB DABCLU (MODULE GRD2A) WHICH IS CALLED BY SUB EFDEN OF
CJMS        MODULE EFGRD2
CJMS
CJMS     3. SUB PAR2PDM (MODULE MP2DDI) WHICH IS CALLED BY SUB PJKDMP2
CJMS        OF MODULE MP2DDI
C
      COMMON /xGSPG80/ E(0:3,0:3,0:3,0:3)
C
CJMS  LABELLED COMMON JMSG80 IS ONLY USED IN THIS MODULE GRD2B.
C
      LOGICAL IJSAME,KLSAME,IKSMJL
      COMMON /xJMSG80/PI4,PIF,P12(3,3),R12,P34(3,3),R34,FCS(4,3),E34MAX,
     2                TX21(MXG2),TH21(MXG2),TY01(MXG2),TY1Y(MXG2),
     3                TF12(MXG2),P0(3,MXG2),TE12(MXG2),CSMAB(MXG2),
     4                TX43(MXG2),TH43(MXG2),TY03(MXG2),TY3Y(MXG2),
     5                TF34(MXG2),Q0(3,MXG2),TE34(MXG2),CSMCD(MXG2),
     6                IGBEG,JGBEG,KGBEG,LGBEG,IGEND,JGEND,KGEND,LGEND,
     7                IJSAME,KLSAME,IKSMJL
C
      PARAMETER (ZER=0.0D+00, PT5=0.5D+00, ONE=1.0D+00)
CC    PARAMETER (PIDIV4=0.785398163397448309615660845819D+00)
      PARAMETER (PIDIV4=0.785398163397448D+00)
      PARAMETER (PITO52=0.349868366552497D+02)
C
CJMS  THRESHOLD VALUES TO USE 'GRDPER' INSTEAD OF THE USUAL 'GRDG80'
C
      LOGICAL   CONCPT
      DIMENSION LIMG80(3:6)
      DATA      LIMG80/ 8, 8,20,20/
C
CC    PIDIV4= ATAN(ONE)
      PI4= PIDIV4
CC    PI = PI4+PI4+PI4+PI4
CC    PITO52=(PI +PI )*PI * SQRT(PI )
      PIF= PITO52
C
C     ----- TWO ELECTRON CONTRIBUTION TO THE GRADIENT -----
C
      IJSAME= INEW.EQ.JNEW
      KLSAME= KNEW.EQ.LNEW
      IKSMJL=(INEW.EQ.KNEW).AND.(JNEW.EQ.LNEW)
      IGBEG= KSTART(INEW)
      JGBEG= KSTART(JNEW)
      KGBEG= KSTART(KNEW)
      LGBEG= KSTART(LNEW)
      IGEND= IGBEG+KNG(INEW)-1
      JGEND= JGBEG+KNG(JNEW)-1
      KGEND= KGBEG+KNG(KNEW)-1
      LGEND= LGBEG+KNG(LNEW)-1
C
      MFQ= JTYPE-1
      IF(JTYPE.LE.3) MFQ= JTYPE
C
      R12= ZER
      R34= ZER
      DO 110 M=1,3
         P12(M,1)= C(M,IAT)
         P12(M,2)= C(M,JAT)
         P12(M,3)= P12(M,2)-P12(M,1)
      R12= R12+P12(M,3)*P12(M,3)
         P34(M,1)= C(M,KAT)
         P34(M,2)= C(M,LAT)
         P34(M,3)= P34(M,2)-P34(M,1)
  110 R34= R34+P34(M,3)*P34(M,3)
C
C     ----- PRELIMINARY SCAN OF MAGNITUDE -----
C
      LK= 0
      E34MAX=ZER
      K0= 1
      DO 125 KG=KGBEG,KGEND
         X03= EX(KG)
         CSUMC= ABS(CS(KG))+ ABS(CP(KG))
         LGND=LGEND
         IF(KLSAME) LGND=KG
         KL=K0
         DO 120 LG=LGBEG,LGND
            X04= EX(LG)
            X34= X03+X04
            X43= ONE/X34
            H43= PT5*X43
            Y03= X03*X43
            Y04= ONE-Y03
            Y34= Y03*X04
            TE34(KL)= ZER
            R34Y34= R34*Y34
            IF(R34Y34.GT.CUX) GO TO 120
            E34= X43* EXP(-R34Y34)
            IF(KLSAME.AND.(KG.NE.LG)) E34=E34+E34
            TE34(KL)= E34
            E34= E34*CSUMC*( ABS(CS(LG))+ ABS(CP(LG)))
            IF(E34MAX.LT.E34) E34MAX=E34
            CSMCD(KL)=E34*E34
            TX43(KL)= X43
            TH43(KL)= H43
            TY03(KL)= Y03
            TY3Y(KL)= Y03*Y04
            TF34(KL)= Y34+Y34
            DO 115 M=1,3
  115       Q0(M,KL)= P34(M,1)+P34(M,3)*Y04
C
            LK=LK+1
C
  120    KL=KL+1
  125 K0=K0+MXGSH
      IF(E34MAX*DABMAX.LT.VTOL1) GO TO 999
C
      JI= 0
      I0= 1
      DO 940 IG=IGBEG,IGEND
         X01= EX(IG)
         CSUMA= DABMAX*( ABS(CS(IG))+ ABS(CP(IG)))
         JGND=JGEND
         IF(IJSAME) JGND=IG
         IJ=I0
         DO 930 JG=JGBEG,JGND
            X02= EX(JG)
            X12= X01+X02
            X21= ONE/X12
            H21= PT5*X21
            Y01= X01*X21
            Y02= ONE-Y01
            Y12= Y01*X02
            TE12(IJ)= ZER
            R12Y12= R12*Y12
            IF(R12Y12.GT.CUX) GO TO 930
            E12= PIF*X21* EXP(-R12Y12)
            IF(IJSAME.AND.(IG.NE.JG)) E12=E12+E12
            CAB= E12*CSUMA*( ABS(CS(JG))+ ABS(CP(JG)))
            IF(E34MAX*CAB.LT.VTOL2) GO TO 930
            TE12(IJ)= E12
            CSMAB(IJ)=CAB*CAB
            TX21(IJ)= X21
            TH21(IJ)= H21
            TY01(IJ)= Y01
            TY1Y(IJ)= Y01*Y02
            TF12(IJ)= Y12+Y12
            DO 152 M=1,3
  152       P0(M,IJ)= P12(M,1)+P12(M,3)*Y02
C
            JI=JI+1
C
  930    IJ=IJ+1
  940 I0=I0+MXGSH
      IF(JI.EQ.0) GO TO 999
C
      DO 135 N=1,3
         DO 135 M=1,4
  135 FCS(M,N)= ZER
C
      CONCPT=.TRUE.
      IF(JTYPE.LT.3) THEN
C FOR JTYPE= 1 OR 2, USE THE 'GRDG80' SUBROUTINE
      ELSEIF(JTYPE.EQ.3) THEN
         IF(JI*LK.GE.LIMG80(JTYPE)) CONCPT=.FALSE.
      ELSEIF(JTYPE.EQ.4) THEN
         IF(JI*LK.GE.LIMG80(JTYPE)) CONCPT=.FALSE.
      ELSEIF(JTYPE.EQ.5) THEN
         IF(JI*LK.GE.LIMG80(JTYPE)) CONCPT=.FALSE.
      ELSEIF(JTYPE.EQ.6) THEN
         IF(JI*LK.GE.LIMG80(JTYPE)) CONCPT=.FALSE.
      ENDIF
      IF(CONCPT) THEN
         CALL xGRDG80(E,JTYPE,MFQ)
      ELSE
         CALL xGRDPER(E,JTYPE,MFQ)
      ENDIF
C
C SUMMATION OF THE CONTRIBUTIONS FROM THE SHELLS
C
      DO 990 M=1,3
         FCS(4,M)=-(FCS(1,M)+FCS(2,M)+FCS(3,M))
         FXYZ(M,IAT)= FXYZ(M,IAT)+FCS(1,M)
         FXYZ(M,JAT)= FXYZ(M,JAT)+FCS(2,M)
         FXYZ(M,KAT)= FXYZ(M,KAT)+FCS(3,M)
         FXYZ(M,LAT)= FXYZ(M,LAT)+FCS(4,M)
  990 CONTINUE
      IF(OUT) THEN
         WRITE(IW,9010) INEW,JNEW,KNEW,LNEW,IAT,JAT,KAT,LAT,
     2                 ((FCS(M,N),M=1,4),N=1,3)
         WRITE(IW,9020)(FXYZ(M,IAT),FXYZ(M,JAT),
     2                  FXYZ(M,KAT),FXYZ(M,LAT),M=1,3)
      ENDIF
  999 CONTINUE
C
      RETURN
 9010 FORMAT(' SHELLS AND CENTERS ',8I4/' ---- PARTIAL CONTRIBUTION'/
     2       4D15.7/4D15.7/4D15.7)
 9020 FORMAT(' ---- GRADIENT UP TO THIS POINT'/
     2       4D15.7/4D15.7/4D15.7)
      END
c--------------------------------------------------------------------------
C*MODULE GRD2B   *DECK GRDG80
      SUBROUTINE xGRDG80(E,JTYPE,MFQ)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION E(0:3,0:3,0:3,0:3)
C
C     ----- GRDG80 DRIVES 2-E INTEGRAL DERIVATIVES FOR SP BASES -----
C     FROM GAUSSIAN-80, MICHEL DUPUIS, 6/6/85.
C     FROM HONDO8, BRETT BODE, 7/15/94.
C     REWRITTEN USING LOOPS AND INLINING FPPPP, JOSE SIERRA, 7/2001
C
#include "gamess.par"
C
      LOGICAL OUT,DBG
      COMMON /xDSHLT / RTOL,DTOL,VTOL1,VTOL2,VTOLS,OUT,DBG
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     2                CF(MXGTOT),CG(MXGTOT),
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /xINDD80/ LA,LB,LC,LD
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
CJMS  LABELLED COMMON JMSG80 IS ONLY USED IN THIS MODULE GRD2B.
C
      LOGICAL IJSAME,KLSAME,IKSMJL
      COMMON /xJMSG80/PI4,PIF,P12(3,3),R12,P34(3,3),R34,FCS(4,3),E34MAX,
     2                TX21(MXG2),TH21(MXG2),TY01(MXG2),TY1Y(MXG2),
     3                TF12(MXG2),P0(3,MXG2),TE12(MXG2),CSMAB(MXG2),
     4                TX43(MXG2),TH43(MXG2),TY03(MXG2),TY3Y(MXG2),
     5                TF34(MXG2),Q0(3,MXG2),TE34(MXG2),CSMCD(MXG2),
     6                IGBEG,JGBEG,KGBEG,LGBEG,IGEND,JGEND,KGEND,LGEND,
     7                IJSAME,KLSAME,IKSMJL
C
      DIMENSION       G(0:3,0:3,0:3,0:3),V(0:3,0:3,0:3,0:3),Z(3),
     2                PIJM0(3),PIJ0N(3),PIJ(3,3),
     3                PKLM0(3),PKL0N(3),PKL(3,3),
     4                P1(3),P2(3,3),FQ(0:5),VEA(12),WEA(12),
     5                VES(9),WES(6),WS(8),WP(8)
C
      PARAMETER (ZER=0.0D+00, PT5=0.5D+00, ONE=1.0D+00, TWO=2.0D+00)
C
      DO 130 M=0,5
  130 FQ(M)= ZER
      DO 140 N=1,12
         VEA(N)= ZER
         WEA(N)= ZER
  140 CONTINUE
C
      CSP= ZER
      CPP= ZER
      DO 145 N=1,3
         DO 145 M=1,3
  145 PIJ(M,N)= ZER
C
C LOOP OVER THE UNCONTRACTED GAUSSIANS WITHIN THE SHELLS
C
C ... GAUSSIANS AT CENTER A.
C
C     ----- I PRIM -----
C
      I0= 1
      DO 940 IG=IGBEG,IGEND
         X01= EX(IG)
         CSA= CS(IG)
         CPA= CP(IG)
C
C ... GAUSSIANS AT CENTER B.
C
C     ----- J PRIM -----
C
         JGND=JGEND
         IF(IJSAME) JGND=IG
         IJ=I0
         DO 930 JG=JGBEG,JGND
            IF(TE12(IJ).EQ.ZER) GO TO 930
            X02= EX(JG)
            CSB= CS(JG)
            CPB= CP(JG)
            X12= X01+X02
            X21= TX21(IJ)
            H21= TH21(IJ)
            Y01= TY01(IJ)
            F12= TF12(IJ)
            X24= X12+X12
            CSB= CSB*TE12(IJ)
            CPB= CPB*TE12(IJ)
            CSS= CSA*CSB
            CPS= CPA*CSB
            CSP= CSA*CPB
            CPP= CPA*CPB
            DO 152 M=1,3
               PIJM0(M)= P0(M,IJ)-P12(M,1)
               PIJ0N(M)= P0(M,IJ)-P12(M,2)
  152       CONTINUE
            IF(LA.NE.1 .AND. LB.NE.1) THEN
               CSP= CSA*CPB
               CPP= CPA*CPB
               DO 156 N=1,3
                  DO 154 M=1,N
  154             PIJ(M,N)= PIJM0(M)*PIJ0N(N)
  156          PIJ(N,N)= PIJ(N,N)+H21
               PIJ(2,1)= PIJ(1,2)
               PIJ(3,1)= PIJ(1,3)
               PIJ(3,2)= PIJ(2,3)
            ENDIF
C
            VE0SUM= ZER
            DO 158 N=1,9
  158       VES(N)= ZER
C
C ... GAUSSIANS AT CENTER C.
C
C     ----- K PRIM -----
C
            KGND=KGEND
            IF(IKSMJL) KGND=IG
            K0= 1
            DO 920 KG=KGBEG,KGND
               X03= EX(KG)
               CSC= CS(KG)
               CPC= CP(KG)
               CSSS= CSS*CSC
               CPSS= CPS*CSC
               CSPS= CSP*CSC
               CPPS= CPP*CSC
               CSSP= CSS*CPC
               CPSP= CPS*CPC
               CSPP= CSP*CPC
               CPPP= CPP*CPC
C
C ... GAUSSIANS AT CENTER D.
C
C     ----- L PRIM -----
C
               LGND=LGEND
               IF(KLSAME) LGND=KG
               IF(IKSMJL.AND.(IG.EQ.KG)) LGND=JG
               KL=K0
               DO 910 LG=LGBEG,LGND
                  IF(TE34(KL).EQ.ZER) GO TO 910
                  X04= EX(LG)
                  CSD= CS(LG)
                  CPD= CP(LG)
                  X34= X03+X04
                  X14= X12+X34
                  X41= ONE/X14
                  VTEST= CSMAB(IJ)*CSMCD(KL)*X41
                  IF(VTEST.LT.VTOLS) GO TO 910
                  H41= PT5*X41
                  Z02= X12*X41
                  Z04= ONE-Z02
C
                  X43= TX43(KL)
                  H43= TH43(KL)
                  Y03= TY03(KL)
C
CJMS              RHO= X12*X34*X41
                  RHO= X12*Z04
                  XVA= ZER
                  DO 160 M=1,3
                     T00 = P0(M,IJ)-Q0(M,KL)
                     P1(M)=T00*Z04
  160             XVA= XVA+T00*T00
                  XVA= XVA*RHO
                  IF(XVA.LT.CUG) THEN
                     TMP= XVA*RGI
                     N  = INT(TMP)
                     C02= TMP-N
                     TMP= C02*(ONE-C02)
                     C03=-TMP*(TWO-C02)
                     C04=-TMP*(ONE+C02)
                     DO 165 M=0,MFQ
                        FQ(M)= R(1,M,N)+C02*R(2,M,N)+C03*R(3,M,N)+
     2                         C04*R(3,M,N+1)
  165                CONTINUE
                  ELSE
                     IF(VTEST.LT.VTOLS*XVA) GO TO 910
                     XIN= ONE/XVA
                     FQ(0)= SQRT(PI4*XIN)
                     TMP= PT5*XIN
                     DO 170 M=1,MFQ
                        FQ(M)= FQ(M-1)*TMP
  170                TMP= TMP+XIN
                  ENDIF
C
                  F34= TF34(KL)
                  E34= TE34(KL)* SQRT(X41)
                  IF(IKSMJL.AND.(IJ.NE.KL)) E34=E34+E34
                  CSD= CSD*E34
                  CPD= CPD*E34
                  WS(1)= CSSS*CSD
                  WS(2)= CPSS*CSD
                  ES1= E(0,0,0,0)*WS(1)
C
                  IF(JTYPE.LT.3) THEN
                     IF(JTYPE.EQ.1) THEN
C
C SPECIAL (SS,SS) SECTION
C
                        DO 175 I=1,6
  175                   WES(I)= ZER
                        VE0= ES1*FQ(0)
                        VE1= ES1*FQ(1)
                     ELSEIF(JTYPE.EQ.2) THEN
C
C SPECIAL (PS,SS) SECTION
C
                        T01=-Z04*FQ(1)
                        ES2= ZER
                        DO 180 I=1,3
                           TMP= E(I,0,0,0)*WS(2)
                           VES(  I)=VES(  I)+TMP*FQ(0)
                           WES(  I)=ZER
                           WES(3+I)=TMP*T01
                           ES1= ES1+TMP*PIJM0(I)
                           ES2= ES2+TMP*P1(I)
  180                   CONTINUE
                        VE0= ES1*FQ(0)-ES2*FQ(1)
                        VE1= ES1*FQ(1)-ES2*FQ(2)
                     ENDIF
                     VE0SUM= VE0SUM+VE0
                     FE0= F34*VE0
                     FE1= X24*VE1
                     GO TO 900
                  ENDIF
C
                  DO 185 N=1,3
                     PKLM0(N)= Q0(N,KL)-P34(N,1)
                     PKL0N(N)= Q0(N,KL)-P34(N,2)
                     DO 185 M=1,N
  185             P2(M,N)= P1(M)*P1(N)
                  P2(2,1)= P2(1,2)
                  P2(3,1)= P2(1,3)
                  P2(3,2)= P2(2,3)
C
                  IF(LC.NE.1 .AND. LD.NE.1) THEN
                     DO 195 N=1,3
                        DO 190 M=1,N
  190                   PKL(M,N)= PKLM0(M)*PKL0N(N)
  195                PKL(N,N)= PKL(N,N)+H43
                     PKL(2,1)= PKL(1,2)
                     PKL(3,1)= PKL(1,3)
                     PKL(3,2)= PKL(2,3)
                  ENDIF
                  XX1=-X12*X43
                  XX2= XX1*XX1
                  XXI=-X21*X34
                  HF2= H41*FQ(1)
                  IF(JTYPE.NE.4) THEN
                     WS(3)= CSPS*CSD
                     WS(4)= CPPS*CSD
                  ENDIF
                  IF(JTYPE.GT.3) THEN
                     XF2= XX1*FQ(1)
                     XF3= XX1*FQ(2)
                     XF4= XX1*FQ(3)
                     WS(5)= CSSP*CSD
                     WS(6)= CPSP*CSD
                     IF(JTYPE.GT.4) THEN
                        WS(7)= CSPP*CSD
                        WS(8)= CPPP*CSD
                        IF(JTYPE.GT.5) THEN
                           WP(1)= CSSS*CPD
                           WP(2)= CPSS*CPD
                           WP(3)= CSPS*CPD
                           WP(4)= CPPS*CPD
                           WP(5)= CSSP*CPD
                           WP(6)= CPSP*CPD
                           WP(7)= CSPP*CPD
                           WP(8)= CPPP*CPD
                        ENDIF
                     ENDIF
                  ENDIF
C
C BEGINING OF THE TWO PASSES THROUGH THE INTEGRAL CALCULATIONS
C
                  DO 800 ICP=1,2
                     HF1= HF2
                     HF2= H41*FQ(2)
                     IF(JTYPE.GT.3) THEN
                        XF1= XF2
                        XF2= XF3
                        XF3= XF4
                        IF(JTYPE.GT.4) XF4= XX1*FQ(4)
                     ENDIF
C
C INTEGRAL EVALUATION SECTION
C
C (SS,SS) SECTION
C
CC    G(0,0,0,0)= FQ(0)
CC    V(0,0,0,0)= FQ(0)
      VE0= ES1*FQ(0)
C
C (PS,SS) SECTION
C
      T20= WS(2)*FQ(0)
      DO 210 I=1,3
         G(I,0,0,0)=-P1(I)*FQ(1)
         V(I,0,0,0)= G(I,0,0,0)+PIJM0(I)*FQ(0)
         VEA(  I)  = E(I,0,0,0)*T20
  210 CONTINUE
      VE0= VE0+WS(2)*(
     1V(1,0,0,0)*E(1,0,0,0)+V(2,0,0,0)*E(2,0,0,0)+V(3,0,0,0)*E(3,0,0,0))
C
      IF(JTYPE.EQ.4) GO TO 400
C
C (PP,SS) + (SP,SS) SECTION
C
      T01= XXI*HF1
      DO 320 J=1,3
         V(0,J,0,0)= V(J,0,0,0)-P12(J,3)*FQ(0)
         T00= T01
         DO 310 I=J,3
            G(I,J,0,0)= T00+P2(I,J)*FQ(2)
            T= PIJ(I,J)*FQ(0)+G(I,J,0,0)
            IF(I.EQ.J) THEN
               V(I,J,0,0)=(PIJM0(J)+PIJ0N(J))*G(J,0,0,0)+T
            ELSE
               V(I,J,0,0)= PIJM0(I)*G(J,0,0,0)+PIJ0N(J)*G(I,0,0,0)+T
               G(J,I,0,0)= G(I,J,0,0)
               V(J,I,0,0)= PIJM0(J)*G(I,0,0,0)+PIJ0N(I)*G(J,0,0,0)+T
            ENDIF
  310    T00= ZER
  320 CONTINUE
      VX3=
     1V(0,1,0,0)*E(0,1,0,0)+V(0,2,0,0)*E(0,2,0,0)+V(0,3,0,0)*E(0,3,0,0)
      T30= WS(3)*FQ(0)
      VX4= ZER
      DO 340 I=1,3
         VEA(3+I)= E(0,I,0,0)*T30+WS(4)*(
     1V(1,0,0,0)*E(1,I,0,0)+V(2,0,0,0)*E(2,I,0,0)+V(3,0,0,0)*E(3,I,0,0))
         VEA(  I)= VEA(  I)+WS(4)*(
     1V(0,1,0,0)*E(I,1,0,0)+V(0,2,0,0)*E(I,2,0,0)+V(0,3,0,0)*E(I,3,0,0))
  340 VX4= VX4+
     1V(I,1,0,0)*E(I,1,0,0)+V(I,2,0,0)*E(I,2,0,0)+V(I,3,0,0)*E(I,3,0,0)
      VE0= VE0+WS(3)*VX3+WS(4)*VX4
C
      IF(JTYPE.EQ.3) GO TO 700
C
  400 CONTINUE
C
C (PS,PS) + (SS,PS) SECTION
C
      DO 410 J=1,3
         G(0,0,J,0)=-P1(J)*XF1
         V(0,0,J,0)= PKLM0(J)*FQ(0)+G(0,0,J,0)
  410 CONTINUE
      DO 430 K=1,3
         T00= HF1
         DO 420 I=K,3
            G(I,0,K,0)= T00+P2(I,K)*XF2
            V(I,0,K,0)= PIJM0(I)*V(0,0,K,0)+G(I,0,0,0)*PKLM0(K)+
     2                  G(I,0,K,0)
            IF(I.NE.K) THEN
CJMS  EQUIVALENCES BY JOSE M. SIERRA
               G(K,0,I,0)= G(I,0,K,0)
               V(K,0,I,0)= PIJM0(K)*V(0,0,I,0)+G(K,0,0,0)*PKLM0(I)+
     2                     G(I,0,K,0)
            ENDIF
  420    T00= ZER
  430 CONTINUE
      VX5=
     1V(0,0,1,0)*E(0,0,1,0)+V(0,0,2,0)*E(0,0,2,0)+V(0,0,3,0)*E(0,0,3,0)
      T50= WS(5)*FQ(0)
      VX6= ZER
      DO 440 I=1,3
         VEA(6+I)= E(0,0,I,0)*T50+WS(6)*(
     1V(1,0,0,0)*E(1,0,I,0)+V(2,0,0,0)*E(2,0,I,0)+V(3,0,0,0)*E(3,0,I,0))
         VEA(  I)= VEA(  I)+WS(6)*(
     1V(0,0,1,0)*E(I,0,1,0)+V(0,0,2,0)*E(I,0,2,0)+V(0,0,3,0)*E(I,0,3,0))
  440 VX6= VX6+
     1V(I,0,1,0)*E(I,0,1,0)+V(I,0,2,0)*E(I,0,2,0)+V(I,0,3,0)*E(I,0,3,0)
      VE0= VE0+WS(5)*VX5+WS(6)*VX6
C
      IF(JTYPE.EQ.4) GO TO 700
C
C (PP,PS) + (SP,PS) SECTION
C
      DO 510 K=1,3
         DO 510 J=1,3
  510 V(0,J,K,0)= V(0,J,0,0)*PKLM0(K)+PIJ0N(J)*G(0,0,K,0)+G(J,0,K,0)
      VX7= ZER
      DO 520 J=1,3
         WEA(6+J)= WEA(6+J)+
     1V(0,1,0,0)*E(0,1,J,0)+V(0,2,0,0)*E(0,2,J,0)+V(0,3,0,0)*E(0,3,J,0)
         WEA(3+J)= WEA(3+J)+
     1V(0,0,1,0)*E(0,J,1,0)+V(0,0,2,0)*E(0,J,2,0)+V(0,0,3,0)*E(0,J,3,0)
  520 VX7= VX7+
     1V(0,J,1,0)*E(0,J,1,0)+V(0,J,2,0)*E(0,J,2,0)+V(0,J,3,0)*E(0,J,3,0)
      VE0= VE0+WS(7)*VX7
      DO 530 N=4,9
         VEA(N)= VEA(N)+WS(7)*WEA(N)
         WEA(N)= ZER
  530 CONTINUE
      EQ0=-P1(3)*P2(2,1)*XF3
      DO 540 I=1,3
         TMP= P2(I,I)*XF3+HF2
         G(I,I,I,0)=-P1(I)*(TMP+HF2+HF2)
         DO 540 J=1,3
            IF(I.NE.J) THEN
               K=6-I-J
CJMS  EQUIVALENCES FROM FPPPP
               G(I,J,K,0)= EQ0
               G(I,K,J,0)= EQ0
CJMS  EQUIVALENCES FROM FPPPP
               EQ2=-P1(J)* TMP
               G(J,I,I,0)= EQ2
               G(I,J,I,0)= EQ2
               G(I,I,J,0)= EQ2
            ENDIF
  540 CONTINUE
      DO 550 K=1,3
         DO 550 J=1,3
            DO 550 I=1,3
  550 V(I,J,K,0)= V(I,J,0,0)*PKLM0(K)+PIJ(I,J)*G(0,0,K,0)+
     2            PIJM0(I)*G(J,0,K,0)+PIJ0N(J)*G(I,0,K,0)+G(I,J,K,0)
      VX8= ZER
      DO 560 J=1,3
         DO 560 I=1,3
            WEA(6+J)= WEA(6+J)+
     1V(I,1,0,0)*E(I,1,J,0)+V(I,2,0,0)*E(I,2,J,0)+V(I,3,0,0)*E(I,3,J,0)
            WEA(3+J)= WEA(3+J)+
     1V(I,0,1,0)*E(I,J,1,0)+V(I,0,2,0)*E(I,J,2,0)+V(I,0,3,0)*E(I,J,3,0)
            WEA(  I)= WEA(  I)+
     1V(0,J,1,0)*E(I,J,1,0)+V(0,J,2,0)*E(I,J,2,0)+V(0,J,3,0)*E(I,J,3,0)
  560 VX8= VX8+
     1V(I,J,1,0)*E(I,J,1,0)+V(I,J,2,0)*E(I,J,2,0)+V(I,J,3,0)*E(I,J,3,0)
      VE0= VE0+WS(8)*VX8
      DO 570 N=1,9
         VEA(N)= VEA(N)+WS(8)*WEA(N)
         WEA(N)= ZER
  570 CONTINUE
C
      IF(JTYPE.EQ.5) GO TO 700
C
C (PP,PP), (PP,SP), (PS,PP), (SP,PP), (SS,PP), (SP,SP), (PS,SP) AND
C (SS,SP) SECTION
C
      V(0,0,0,1)= V(0,0,1,0)-P34(1,3)*FQ(0)
      V(0,0,0,2)= V(0,0,2,0)-P34(2,3)*FQ(0)
      V(0,0,0,3)= V(0,0,3,0)-P34(3,3)*FQ(0)
      VE0= VE0+WP(1)*(
     1V(0,0,0,1)*E(0,0,0,1)+V(0,0,0,2)*E(0,0,0,2)+V(0,0,0,3)*E(0,0,0,3))
      DO 605 L=1,3
         DO 605 I=1,3
  605 V(I,0,0,L)= V(I,0,0,0)*PKL0N(L)+PIJM0(I)*G(0,0,L,0)+G(I,0,L,0)
      T10= WP(1)*FQ(0)
      VXE= ZER
      DO 610 L=1,3
         VEA(9+L)= E(0,0,0,L)*T10+WP(2)*(
     1V(1,0,0,0)*E(1,0,0,L)+V(2,0,0,0)*E(2,0,0,L)+V(3,0,0,0)*E(3,0,0,L))
         VEA(  L)= VEA(  L)+WP(2)*(
     1V(0,0,0,1)*E(L,0,0,1)+V(0,0,0,2)*E(L,0,0,2)+V(0,0,0,3)*E(L,0,0,3))
  610 VXE= VXE+
     1V(L,0,0,1)*E(L,0,0,1)+V(L,0,0,2)*E(L,0,0,2)+V(L,0,0,3)*E(L,0,0,3)
      VE0= VE0+WP(2)*VXE
      DO 615 L=1,3
         DO 615 J=1,3
  615 V(0,J,0,L)= V(0,J,0,0)*PKL0N(L)+PIJ0N(J)*G(0,0,L,0)+G(J,0,L,0)
      VXE= ZER
      DO 620 J=1,3
         VEA(9+J)= VEA(9+J)+WP(3)*(
     1V(0,1,0,0)*E(0,1,0,J)+V(0,2,0,0)*E(0,2,0,J)+V(0,3,0,0)*E(0,3,0,J))
         VEA(3+J)= VEA(3+J)+WP(3)*(
     1V(0,0,0,1)*E(0,J,0,1)+V(0,0,0,2)*E(0,J,0,2)+V(0,0,0,3)*E(0,J,0,3))
  620 VXE= VXE+
     1V(0,J,0,1)*E(0,J,0,1)+V(0,J,0,2)*E(0,J,0,2)+V(0,J,0,3)*E(0,J,0,3)
      VE0= VE0+WP(3)*VXE
      DO 625 L=1,3
         DO 625 J=1,3
            DO 625 I=1,3
  625 V(I,J,0,L)= V(I,J,L,0)-V(I,J,0,0)*P34(L,3)
      VXE= ZER
      DO 630 J=1,3
         DO 630 I=1,3
            WEA(9+J)= WEA(9+J)+
     1V(I,1,0,0)*E(I,1,0,J)+V(I,2,0,0)*E(I,2,0,J)+V(I,3,0,0)*E(I,3,0,J)
            WEA(3+J)= WEA(3+J)+
     1V(I,0,0,1)*E(I,J,0,1)+V(I,0,0,2)*E(I,J,0,2)+V(I,0,0,3)*E(I,J,0,3)
            WEA(  I)= WEA(  I)+
     1V(0,J,0,1)*E(I,J,0,1)+V(0,J,0,2)*E(I,J,0,2)+V(0,J,0,3)*E(I,J,0,3)
  630 VXE= VXE+
     1V(I,J,0,1)*E(I,J,0,1)+V(I,J,0,2)*E(I,J,0,2)+V(I,J,0,3)*E(I,J,0,3)
      VE0= VE0+WP(4)*VXE
      DO 635 N=1,3
         VEA(  N)= VEA(  N)+WP(4)*WEA(  N)
         WEA(  N)= ZER
         VEA(3+N)= VEA(3+N)+WP(4)*WEA(3+N)
         WEA(3+N)= ZER
         VEA(9+N)= VEA(9+N)+WP(4)*WEA(9+N)
         WEA(9+N)= ZER
  635 CONTINUE
      DO 640 L=1,3
         DO 640 K=L,3
            T= PKL(K,L)*FQ(0)+G(K,L,0,0)*XX2
            IF(K.EQ.L) THEN
               V(0,0,K,L)= T+G(0,0,K,0)*(PKLM0(K)+PKL0N(K))
            ELSE
               V(0,0,K,L)= T+PKLM0(K)*G(0,0,L,0)+PKL0N(L)*G(0,0,K,0)
               V(0,0,L,K)= T+PKLM0(L)*G(0,0,K,0)+PKL0N(K)*G(0,0,L,0)
            ENDIF
  640 CONTINUE
      VXE= ZER
      DO 645 K=1,3
         VEA(9+K)= VEA(9+K)+WP(5)*(
     1V(0,0,1,0)*E(0,0,1,K)+V(0,0,2,0)*E(0,0,2,K)+V(0,0,3,0)*E(0,0,3,K))
         VEA(6+K)= VEA(6+K)+WP(5)*(
     1V(0,0,0,1)*E(0,0,K,1)+V(0,0,0,2)*E(0,0,K,2)+V(0,0,0,3)*E(0,0,K,3))
  645 VXE= VXE+
     1V(0,0,K,1)*E(0,0,K,1)+V(0,0,K,2)*E(0,0,K,2)+V(0,0,K,3)*E(0,0,K,3)
      VE0= VE0+WP(5)*VXE
C
      T00= XX1*XF4
      T01= H41*XF3
      T02= T01+T01
      T05= T01+T02+T02
      T25= H41*HF2
      T75= T25+T25+T25
      EQ0= G(3,2,1,0)*XX1
      DO 650 I=1,3
         V(I,0,I,I)= G(I,I,I,0)*XX1
         TMP= P2(I,I)* T00+T01
         V(I,I,I,I)= P2(I,I)*(TMP+T05)+T75
         DO 650 J=1,3
            IF(I.NE.J) THEN
               K=6-I-J
CJMS  EQUIVALENCES BY JOSE M. SIERRA
               V(I,0,J,K)= EQ0
               V(I,0,K,J)= EQ0
CJMS  EQUIVALENCES BY JOSE M. SIERRA
               EQ1= G(I,I,J,0)*XX1
               V(J,0,I,I)= EQ1
               V(I,0,J,I)= EQ1
               V(I,0,I,J)= EQ1
CJMS  EQUIVALENCES FROM FPPPP
               EQ2= P2(J,I)*(TMP+T02)
               EQ3= P2(J,J)* TMP+T25 +P2(I,I)*T01
               EQ4= P2(K,J)* TMP
               V(J,I,I,I)= EQ2
               V(I,J,I,I)= EQ2
               V(I,I,J,I)= EQ2
               V(I,I,I,J)= EQ2
               V(J,J,I,I)= EQ3
               V(J,I,J,I)= EQ3
               V(I,J,J,I)= EQ3
               V(J,I,I,J)= EQ3
               V(I,J,I,J)= EQ3
               V(I,I,J,J)= EQ3
               V(K,J,I,I)= EQ4
               V(K,I,J,I)= EQ4
               V(K,I,I,J)= EQ4
               V(J,K,I,I)= EQ4
               V(I,K,J,I)= EQ4
               V(I,K,I,J)= EQ4
               V(J,I,K,I)= EQ4
               V(I,J,K,I)= EQ4
               V(I,I,K,J)= EQ4
               V(J,I,I,K)= EQ4
               V(I,J,I,K)= EQ4
               V(I,I,J,K)= EQ4
            ENDIF
  650 CONTINUE
      DO 657 L=1,3
         DO 657 K=1,3
            T01= PKL(K,L)
            T02= V(0,0,K,L)
            DO 653 I=1,3
               Z(I)= V(I,0,K,L)+G(I,0,0,0)*T01+
     2               G(I,0,K,0)*PKL0N(L)+G(I,0,L,0)*PKLM0(K)
               V(I,0,K,L)= Z(I)+PIJM0(I)*T02
               V(0,I,K,L)= Z(I)+PIJ0N(I)*T02
  653       CONTINUE
            DO 655 J=1,3
               DO 655 I=1,3
                  V(I,J,K,L)=V(I,J,K,L)+ Z(I)*PIJ0N(J)+Z(J)*PIJM0(I)+
     2                       G(I,J,0,0)*T01+PIJ(I,J)*T02+
     3                       G(I,J,K,0)*PKL0N(L)+G(I,J,L,0)*PKLM0(K)
  655       CONTINUE
  657 CONTINUE
C
      VXE= ZER
      DO 660 K=1,3
         DO 660 I=1,3
            WEA(9+K)= WEA(9+K)+
     1V(I,0,1,0)*E(I,0,1,K)+V(I,0,2,0)*E(I,0,2,K)+V(I,0,3,0)*E(I,0,3,K)
            WEA(6+K)= WEA(6+K)+
     1V(I,0,0,1)*E(I,0,K,1)+V(I,0,0,2)*E(I,0,K,2)+V(I,0,0,3)*E(I,0,K,3)
            WEA(  I)= WEA(  I)+
     1V(0,0,K,1)*E(I,0,K,1)+V(0,0,K,2)*E(I,0,K,2)+V(0,0,K,3)*E(I,0,K,3)
  660 VXE= VXE+
     1V(I,0,K,1)*E(I,0,K,1)+V(I,0,K,2)*E(I,0,K,2)+V(I,0,K,3)*E(I,0,K,3)
      VE0= VE0+WP(6)*VXE
      DO 665 N=1,3
         VEA(  N)= VEA(  N)+WP(6)*WEA(  N)
         WEA(  N)= ZER
         VEA(6+N)= VEA(6+N)+WP(6)*WEA(6+N)
         WEA(6+N)= ZER
         VEA(9+N)= VEA(9+N)+WP(6)*WEA(9+N)
         WEA(9+N)= ZER
  665 CONTINUE
      VXE= ZER
      DO 675 K=1,3
         DO 675 J=1,3
            WEA(9+K)= WEA(9+K)+
     1V(0,J,1,0)*E(0,J,1,K)+V(0,J,2,0)*E(0,J,2,K)+V(0,J,3,0)*E(0,J,3,K)
            WEA(6+K)= WEA(6+K)+
     1V(0,J,0,1)*E(0,J,K,1)+V(0,J,0,2)*E(0,J,K,2)+V(0,J,0,3)*E(0,J,K,3)
            WEA(3+J)= WEA(3+J)+
     1V(0,0,K,1)*E(0,J,K,1)+V(0,0,K,2)*E(0,J,K,2)+V(0,0,K,3)*E(0,J,K,3)
  675 VXE= VXE+
     1V(0,J,K,1)*E(0,J,K,1)+V(0,J,K,2)*E(0,J,K,2)+V(0,J,K,3)*E(0,J,K,3)
      VE0= VE0+WP(7)*VXE
      DO 680 N=4,12
         VEA(N)= VEA(N)+WP(7)*WEA(N)
         WEA(N)= ZER
  680 CONTINUE
      VXE= ZER
      DO 685 L=1,3
         DO 685 K=1,3
            DO 685 J=1,3
               DO 685 I=1,3
  685 VXE= VXE+V(I,J,K,L)*E(I,J,K,L)
      VE0= VE0+WP(8)*VXE
      DO 690 K=1,3
         DO 690 J=1,3
            DO 690 I=1,3
               WEA(9+K)= WEA(9+K)+
     1V(I,J,1,0)*E(I,J,1,K)+V(I,J,2,0)*E(I,J,2,K)+V(I,J,3,0)*E(I,J,3,K)
               WEA(6+K)= WEA(6+K)+
     1V(I,J,0,1)*E(I,J,K,1)+V(I,J,0,2)*E(I,J,K,2)+V(I,J,0,3)*E(I,J,K,3)
               WEA(3+J)= WEA(3+J)+
     1V(I,0,K,1)*E(I,J,K,1)+V(I,0,K,2)*E(I,J,K,2)+V(I,0,K,3)*E(I,J,K,3)
               WEA(  I)= WEA(  I)+
     1V(0,J,K,1)*E(I,J,K,1)+V(0,J,K,2)*E(I,J,K,2)+V(0,J,K,3)*E(I,J,K,3)
  690 CONTINUE
      DO 695 N=1,12
         VEA(N)= VEA(N)+WP(8)*WEA(N)
         WEA(N)= ZER
  695 CONTINUE
C
  700                CONTINUE
                     IF(ICP.EQ.1) THEN
                        DO 710 M=1,MFQ
  710                   FQ(M-1)= FQ(M)
                        VE0SUM= VE0SUM+VE0
                        FE0= F34*VE0
                        DO 720 M=1,3
                           VES(  M)= VES(  M)+VEA(  M)
                           VES(3+M)= VES(3+M)+VEA(3+M)
                           T01 = VEA(6+M)
                           T02 = VEA(9+M)
                           WES(  M)=-T01+(T01+T02)*Y03
  720                   CONTINUE
                     ELSE
                        FE1= X24*VE0
                        DO 730 M=1,3
                           T01 = VEA(  M)+VEA(3+M)
                           T02 = VEA(6+M)+VEA(9+M)
                           WES(3+M)=-T01+(T01+T02)*Z02
  730                   CONTINUE
                     ENDIF
  800             CONTINUE
  900             CONTINUE
C
C SUMMATION OF CONTRIBUTIONS FROM THE UNCONTRACTED GAUSSIANS
C
                  DO 905 M=1,3
                     WES(  M)= WES(  M)+P34(M,3)*FE0
                     WES(3+M)= WES(3+M)-P1(M)*FE1
                     VES(6+M)= VES(6+M)+WES(3+M)
                     FCS(3,M)= FCS(3,M)-WES(3+M)*Y03+WES(  M)
  905             CONTINUE
C
  910          KL=KL+1
  920       K0=K0+MXGSH
C
            TMP= F12*VE0SUM
            DO 925 M=1,3
               T01= VES(  M)
               T02= VES(3+M)+VES(6+M)
               T01=-T01+(T01+T02)*Y01+P12(M,3)*TMP
               FCS(1,M)= FCS(1,M)+T01
               FCS(2,M)= FCS(2,M)-T01+VES(6+M)
  925       CONTINUE
C
  930    IJ=IJ+1
  940 I0=I0+MXGSH
C
      RETURN
      END
c--------------------------------------------------------------------------
C*MODULE GRD2B   *DECK GRDPER
      SUBROUTINE xGRDPER(E,JTYPE,MFQ)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION E(0:3,0:3,0:3,0:3)
C
C     ----- GRDPER DRIVES 2-E INTEGRAL DERIVATIVES FOR SP BASES -----
C     PERFORMANCE VERSION, REWRITTEN FROM GRDG80, JOSE SIERRA, 6/2002
C
#include "gamess.par"
C
      LOGICAL OUT,DBG
      COMMON /xDSHLT / RTOL,DTOL,VTOL1,VTOL2,VTOLS,OUT,DBG
      COMMON /xIOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     2                CF(MXGTOT),CG(MXGTOT),
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
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
CJMS  LABELLED COMMON JMSG80 IS ONLY USED IN THIS MODULE GRD2B.
C
      LOGICAL IJSAME,KLSAME,IKSMJL
      COMMON /xJMSG80/PI4,PIF,P12(3,3),R12,P34(3,3),R34,FCS(4,3),E34MAX,
     2                TX21(MXG2),TH21(MXG2),TY01(MXG2),TY1Y(MXG2),
     3                TF12(MXG2),P0(3,MXG2),TE12(MXG2),CSMAB(MXG2),
     4                TX43(MXG2),TH43(MXG2),TY03(MXG2),TY3Y(MXG2),
     5                TF34(MXG2),Q0(3,MXG2),TE34(MXG2),CSMCD(MXG2),
     6                IGBEG,JGBEG,KGBEG,LGBEG,IGEND,JGEND,KGEND,LGEND,
     7                IJSAME,KLSAME,IKSMJL
C
      DIMENSION P1(3),FQ(0:5),F1(0:4),F2(0:4),WS(8),WP(8)
      DIMENSION VEA(12),VES(9),WES(6)
C
      DIMENSION T(0:20),U(0:20)
      DIMENSION P1D(3,3),P1P(3,3),P24(3,3),P3D(3,3),P3P(3,3)
      DIMENSION Q1(6,3,3),R1(6,3,3),R2(6,3,3),R3(0:9,3,3,3)
      DIMENSION T1R(6,3,3,3),T9A(0:13,3,3,3)
      DIMENSION T3R(6,3,3,3),T9B(0:20,3,3,3)
C
      DIMENSION F5X(0:12,3,3,3),F6X(0: 6,3,3,3,3)
C
      DIMENSION A31(0:3,0:3),C31(2,0:3),A32(0:3,0:3),C32(2,0:3)
      DIMENSION A41(0:3,0:3),C41(2,0:3),A42(0:3,0:3),C42(2,0:3)
      DIMENSION A33(16),FIJ00(2),A43(16),FI0K0(2)
      DIMENSION A51(0:3,3),C51(2,  3),A52(0:3,3),C52(2,  3)
      DIMENSION SEIKK0(0:3),A53(16,0:3),C53(2,0:3)
      DIMENSION SEJJK0(  3),A54(16,  3),C54(2,  3)
      DIMENSION SEKJK0(  3),A55(16,  3),C55(2,  3)
      DIMENSION A56(0:22,3,0:3),C56(2,0:3)
C
      DIMENSION A60(0:3,0:3),C60(2,0:3),A61(0:3,0:3),C61(2,0:3)
      DIMENSION A62(16),FI00L(2),A63(16),F0J0L(2)
      DIMENSION A64(0:3,0:3),C64(2,  3),A65(0:3,0:3),C65(2,  3)
      DIMENSION A66(18,  3),C66(2,  3)
      DIMENSION A67(18,  3),C67(2,  3)
      DIMENSION A68(18,  3),C68(2,  3)
      DIMENSION A69(0:3,  3),C69(2,0:3),A70(0:3,  3),C70(2,0:3)
      DIMENSION A71(18,  3),C71(2,  3)
      DIMENSION A72(18,  3),C72(2,  3)
      DIMENSION A73(18,0:3),C73(2,0:3)
C
      DIMENSION SE0KKL(  3),A74(16,3),C74(2,0:3)
      DIMENSION SE0LKL(  3),A75(16,3),C75(2,0:3)
      DIMENSION SE0JLL(  3),A76(16,3),C76(2,0:3)
      DIMENSION A77(0:25,3,0:3),C77(2,0:3),A78(0:31,3,0:3),C78(2,0:3)
      DIMENSION A79(0:31,3,0:3),C79(2,0:3)
      DIMENSION A80(0: 2,2),A81(0:24,3),A82(0:31,2),A83(0:22,2)
      DIMENSION A84(0:13,2),A85(0:13,2),A86(0: 6)
      DIMENSION S4(0:14),Q4(0:4),FIJKL(2)
C
      PARAMETER (ZER=0.0D+00, PT5=0.5D+00, ONE=1.0D+00, TWO=2.0D+00)
      PARAMETER (THR=3.0D+00, SIX=6.0D+00)
C
      IF(JTYPE.LT.3) THEN
         WRITE(*,9010)
cwc         CALL ABRT
cwc1
         call aces_exit(1)
cwc0
      ENDIF
 9010 FORMAT(' GRDPER ONLY VALID FOR JTYPE= 3, 4, 5 OR 6')
C
      DO 140 N=1,12
  140 VEA(N)= ZER
C
C LOOP OVER THE UNCONTRACTED GAUSSIANS WITHIN THE SHELLS
C
C ... GAUSSIANS AT CENTER A.
C
C     ----- I PRIM -----
C
      JI= 0
      I0= 1
      DO 940 IG=IGBEG,IGEND
         X01= EX(IG)
         CSA= CS(IG)
         CPA= CP(IG)
C
C ... GAUSSIANS AT CENTER B.
C
C     ----- J PRIM -----
C
         JGND=JGEND
         IF(IJSAME) JGND=IG
         IJ=I0
         DO 930 JG=JGBEG,JGND
            IF(TE12(IJ).EQ.ZER) GO TO 930
            X02= EX(JG)
            CSB= CS(JG)
            CPB= CP(JG)
            X12= X01+X02
            H21= TH21(IJ)
            Y01= TY01(IJ)
            Y02= ONE-Y01
            Y1Y= TY1Y(IJ)
            F12= TF12(IJ)
            X24= X12+X12
            CSB= CSB*TE12(IJ)
            CPB= CPB*TE12(IJ)
            CSS= CSA*CSB
            CPS= CPA*CSB
            CSP= CSA*CPB
            CPP= CPA*CPB
C
            VE0SUM= ZER
            DO 158 N=1,9
  158       VES(N)= ZER
C
C ... GAUSSIANS AT CENTER C.
C
C     ----- K PRIM -----
C
            KGND=KGEND
            IF(IKSMJL) KGND=IG
            LK= 0
            K0= 1
            DO 920 KG=KGBEG,KGND
               X03= EX(KG)
               CSC= CS(KG)
               CPC= CP(KG)
               CSSS= CSS*CSC
               CPSS= CPS*CSC
               CSPS= CSP*CSC
               CPPS= CPP*CSC
               CSSP= CSS*CPC
               CPSP= CPS*CPC
               CSPP= CSP*CPC
               CPPP= CPP*CPC
C
C ... GAUSSIANS AT CENTER D.
C
C     ----- L PRIM -----
C
               LGND=LGEND
               IF(KLSAME) LGND=KG
               IF(IKSMJL.AND.(IG.EQ.KG)) LGND=JG
               KL=K0
               DO 910 LG=LGBEG,LGND
                  IF(TE34(KL).EQ.ZER) GO TO 910
                  X04= EX(LG)
                  CSD= CS(LG)
                  CPD= CP(LG)
                  X34= X03+X04
                  X14= X12+X34
                  X41= ONE/X14
                  VTEST= CSMAB(IJ)*CSMCD(KL)*X41
                  IF(VTEST.LT.VTOLS) GO TO 910
                  Z02= X12*X41
                  Z04= ONE-Z02
C
                  X43= TX43(KL)
                  H43= TH43(KL)
                  Y03= TY03(KL)
                  Y04= ONE-Y03
                  Y3Y= TY3Y(KL)
C
CJMS              RHO= X12*X34*X41
                  RHO= X12*Z04
                  XVA= ZER
                  DO 160 M=1,3
                     T00 = P0(M,IJ)-Q0(M,KL)
                     P1(M)=T00*Z04
  160             XVA= XVA+T00*T00
                  XVA= XVA*RHO
                  IF(XVA.LT.CUG) THEN
                     TMP= XVA*RGI
                     N  = INT(TMP)
                     C02= TMP-N
                     TMP= C02*(ONE-C02)
                     C03=-TMP*(TWO-C02)
                     C04=-TMP*(ONE+C02)
                     DO 165 M=0,MFQ
                        FQ(M)= R(1,M,N)+C02*R(2,M,N)+C03*R(3,M,N)+
     2                         C04*R(3,M,N+1)
  165                CONTINUE
                  ELSE
                     IF(VTEST.LT.VTOLS*XVA) GO TO 910
                     XIN= ONE/XVA
                     FQ(0)= SQRT(PI4*XIN)
                     TMP= PT5*XIN
                     DO 170 M=1,MFQ
                        FQ(M)= FQ(M-1)*TMP
  170                TMP= TMP+XIN
                  ENDIF
C
C COEFFICIENTS FOR THE TWO PASSES OF THE INTEGRAL CALCULATIONS
C
                  TMP= ONE
                  DO 171 M=0,MFQ-1
                     F1(M)= FQ(M  )*TMP
                     F2(M)= FQ(M+1)*TMP
  171             TMP= TMP*Z04
C
                  DO 175 N=1,6
  175             WES(N)= ZER
C
                  F34= TF34(KL)
                  E34= TE34(KL)* SQRT(X41)
                  IF(IKSMJL.AND.(IJ.NE.KL)) E34=E34+E34
                  CSD= CSD*E34
                  CPD= CPD*E34
                  WS(1)= CSSS*CSD
                  WS(2)= CPSS*CSD
                  ES1= E(0,0,0,0)*WS(1)
C
                  XX1=-X12*X43
                  XX2= XX1*XX1
                  Y31= Y03*XX1
                  Y41= XX1-Y31
                  IF(JTYPE.NE.4) THEN
                     WS(3)= CSPS*CSD
                     WS(4)= CPPS*CSD
                  ENDIF
                  IF(JTYPE.GT.3) THEN
                     WS(5)= CSSP*CSD
                     WS(6)= CPSP*CSD
                     IF(JTYPE.GT.4) THEN
                        WS(7)= CSPP*CSD
                        WS(8)= CPPP*CSD
                        IF(JTYPE.GT.5) THEN
                           WP(1)= CSSS*CPD
                           WP(2)= CPSS*CPD
                           WP(3)= CSPS*CPD
                           WP(4)= CPPS*CPD
                           WP(5)= CSSP*CPD
                           WP(6)= CPSP*CPD
                           WP(7)= CSPP*CPD
                           WP(8)= CPPP*CPD
                        ENDIF
                     ENDIF
                  ENDIF
C
      JI= JI+1
      IF(JI.EQ.1) THEN
C
CJMS  THE BULK CALCULATION IS MADE IN THIS BLOCK.
CJMS  IT IS ONLY CALCULATED ONCE PER SUBROUTINE CALL
C
         DO 210 N=1,3
            D1N= P12(N,1)-P34(N,1)
            DO 210 M=1,3
               P1D(M,N)= P12(M,3)*D1N
               P1P(M,N)= P12(M,3)*P12(N,3)
               P24(M,N)= P12(M,3)*P34(N,3)
               P3D(M,N)= P34(M,3)*D1N
               P3P(M,N)= P34(M,3)*P34(N,3)
  210    CONTINUE
         DO 220 N=1,3
            DO 220 M=1,3
               Q1(1,M,N)= P1D(M,N)
               Q1(2,M,N)= P1P(M,N)
               Q1(3,M,N)= P24(M,N)
               Q1(4,M,N)= P3D(N,M)
               Q1(5,M,N)= P3P(M,N)
               R1(1,M,N)= P1D(M,N)+P1P(M,N)
               R1(2,M,N)= P1D(M,N)+P1D(N,M)+P1P(M,N)+P1P(M,N)+P1P(M,N)
               R1(3,M,N)= P1P(M,N)+P1P(M,N)
               R1(4,M,N)= P24(M,N)
               R1(5,M,N)= P24(M,N)+P24(N,M)
C P2(M,N)
               R2(1,M,N)=(P12(M,1)-P34(M,1))*(P12(N,1)-P34(N,1))
               R2(2,M,N)= P1D(M,N)+P1D(N,M)
               R2(3,M,N)= P1P(M,N)
               R2(4,M,N)= P3D(M,N)+P3D(N,M)
               R2(5,M,N)= P24(M,N)+P24(N,M)
               R2(6,M,N)= P3P(M,N)
  220    CONTINUE
C
         IF(JTYPE.NE.4) THEN
            DO 233 J=0,3
               DO 231 M=0,2
                  A31(M,J)= ZER
                  A32(M,J)= ZER
  231          CONTINUE
               DO 233 I=1,3
                  D1I= P12(I,1)-P34(I,1)
                  A31(0,J)= A31(0,J)+D1I     *E(J,I,0,0)
                  A31(1,J)= A31(1,J)+P12(I,3)*E(J,I,0,0)
                  A31(2,J)= A31(2,J)+P34(I,3)*E(J,I,0,0)
                  A32(0,J)= A32(0,J)+D1I     *E(I,J,0,0)
                  A32(1,J)= A32(1,J)+P12(I,3)*E(I,J,0,0)
                  A32(2,J)= A32(2,J)+P34(I,3)*E(I,J,0,0)
  233       CONTINUE
            AEIJ00= ZER
            DO 234 J=1,3
               DO 234 I=1,3
  234       AEIJ00= AEIJ00+P1P(I,J)*E(I,J,0,0)
            SEJJ00= E(1,1,0,0)+E(2,2,0,0)+E(3,3,0,0)
            DO 235 N=1,4
               A33(  N)= ZER
               A33(8+N)= ZER
  235       CONTINUE
            A33(13)= ZER
            A33(14)= ZER
            DO 239 J=1,3
               DO 239 I=1,3
C V(I,J,0,0)
                  GEIJ00= E(I,J,0,0)
                  DO 237 N=1,4
                     A33(  N)= A33(  N)+R1(N,I,J)*GEIJ00
                     A33(8+N)= A33(8+N)+R2(N,I,J)*GEIJ00
  237             CONTINUE
                  A33(13)= A33(13)+R2(5,I,J)*GEIJ00
                  A33(14)= A33(14)+R2(6,I,J)*GEIJ00
  239       CONTINUE
            A33( 5)= A33(13)
         ENDIF
         IF(JTYPE.LT.4) GO TO 490
C
         DO 243 K=0,3
            DO 241 M=0,2
               A41(M,K)= ZER
               A42(M,K)= ZER
  241       CONTINUE
            DO 243 I=1,3
               D1I= P12(I,1)-P34(I,1)
               A41(0,K)= A41(0,K)+D1I     *E(K,0,I,0)
               A41(1,K)= A41(1,K)+P12(I,3)*E(K,0,I,0)
               A41(2,K)= A41(2,K)+P34(I,3)*E(K,0,I,0)
               A42(0,K)= A42(0,K)+D1I     *E(I,0,K,0)
               A42(1,K)= A42(1,K)+P12(I,3)*E(I,0,K,0)
               A42(2,K)= A42(2,K)+P34(I,3)*E(I,0,K,0)
  243    CONTINUE
         DO 245 N=1,5
            A43(  N)= ZER
            A43(8+N)= ZER
  245    CONTINUE
         SEK0K0= E(1,0,1,0)+E(2,0,2,0)+E(3,0,3,0)
         DO 249 K=1,3
            DO 249 I=1,3
C V(I,0,K,0)
               GEI0K0= E(I,0,K,0)
               DO 247 N=1,5
                  A43(  N)= A43(  N)+Q1(N,I,K)*GEI0K0
                  A43(8+N)= A43(8+N)+R2(N,I,K)*GEI0K0
  247          CONTINUE
  249    CONTINUE
         A43(14)= A43( 5)
         IF(JTYPE.LT.5) GO TO 490
C
         DO 255 K=1,3
            DO 250 M=0,2
               A51(M,K)= ZER
               A52(M,K)= ZER
  250       CONTINUE
            DO 255 J=1,3
               D1J = P12(J,1)-P34(J,1)
               A51(0,K)= A51(0,K)+D1J     *E(0,K,J,0)
               A51(1,K)= A51(1,K)+P12(J,3)*E(0,K,J,0)
               A51(2,K)= A51(2,K)+P34(J,3)*E(0,K,J,0)
               A52(0,K)= A52(0,K)+D1J     *E(0,J,K,0)
               A52(1,K)= A52(1,K)+P12(J,3)*E(0,J,K,0)
               A52(2,K)= A52(2,K)+P34(J,3)*E(0,J,K,0)
  255    CONTINUE
         DO 268 K=1,3
            SEJJK0(K)= E(1,1,K,0)+E(2,2,K,0)+E(3,3,K,0)
            DO 262 N=1,4
               A54(  N,K)= ZER
               A54(8+N,K)= ZER
  262       CONTINUE
            A54(13,K)= ZER
            A54(14,K)= ZER
            DO 266 J=1,3
               DO 266 I=1,3
C V(I,J,0,0)
                  GEIJK0= E(I,J,K,0)
                  DO 264 N=1,4
                     A54(  N,K)= A54(  N,K)+R1(N,I,J)*GEIJK0
                     A54(8+N,K)= A54(8+N,K)+R2(N,I,J)*GEIJK0
  264             CONTINUE
                  A54(13,K)= A54(13,K)+R2(5,I,J)*GEIJK0
                  A54(14,K)= A54(14,K)+R2(6,I,J)*GEIJK0
  266       CONTINUE
            A54( 5,K)= A54(13,K)
  268    CONTINUE
         DO 278 J=1,3
            SEKJK0(J)= E(1,J,1,0)+E(2,J,2,0)+E(3,J,3,0)
            DO 272 N=1,5
               A55(  N,J)= ZER
               A55(8+N,J)= ZER
  272       CONTINUE
            DO 276 K=1,3
               DO 276 I=1,3
C V(I,0,K,0)
                  GEIJK0= E(I,J,K,0)
                  DO 274 N=1,5
                     A55(  N,J)= A55(  N,J)+Q1(N,I,K)*GEIJK0
                     A55(8+N,J)= A55(8+N,J)+R2(N,I,K)*GEIJK0
  274             CONTINUE
  276       CONTINUE
            A55(14,J)= A55( 5,J)
  278    CONTINUE
         DO 288 I=0,3
            SEIKK0(I)= E(I,1,1,0)+E(I,2,2,0)+E(I,3,3,0)
            DO 282 N=1,5
               A53(  N,I)= ZER
               A53(8+N,I)= ZER
  282       CONTINUE
            DO 286 K=1,3
               DO 286 J=1,3
C V(0,J,K,0)
                  GEIJK0= E(I,J,K,0)
                  DO 284 N=1,5
                     A53(  N,I)= A53(  N,I)+Q1(N,J,K)*GEIJK0
                     A53(8+N,I)= A53(8+N,I)+R2(N,J,K)*GEIJK0
  284             CONTINUE
  286       CONTINUE
            A53(14,I)= A53( 5,I)
  288    CONTINUE
C
         DO 293 K=1,3
            DO 293 J=1,3
               DO 293 I=1,3
                  DO 291 M=1,6
                     T1R(M,I,J,K)= P12(I,3)*R2(M,J,K)
                     T3R(M,I,J,K)= P34(I,3)*R2(M,J,K)
  291             CONTINUE
                  T1R(6,I,J,K)=-T1R(6,I,J,K)
                  T3R(6,I,J,K)=-T3R(6,I,J,K)
  293    CONTINUE
C
         DO 295 N=1,3
            DO 295 M=1,3
               DO 295 L=1,3
C P3(L,M,N)
                  D1L= P12(L,1)-P34(L,1)
                  R3(0,L,M,N)= D1L*R2(1,M,N)
                  R3(1,L,M,N)= D1L*R2(2,M,N)+T1R(1,L,M,N)
                  R3(2,L,M,N)= D1L*R2(3,M,N)+T1R(2,L,M,N)
                  R3(3,L,M,N)= T1R(3,L,M,N)
                  R3(4,L,M,N)=-D1L*R2(4,M,N)-T3R(1,L,M,N)
                  R3(5,L,M,N)=-D1L*R2(5,M,N)-T1R(4,L,M,N)-T3R(2,L,M,N)
                  R3(6,L,M,N)=-T1R(5,L,M,N) -T3R(3,L,M,N)
                  R3(7,L,M,N)= D1L*R2(6,M,N)+T3R(4,L,M,N)
                  R3(8,L,M,N)=-T1R(6,L,M,N) +T3R(5,L,M,N)
                  R3(9,L,M,N)= T3R(6,L,M,N)
  295    CONTINUE
C
         DO 299 I=1,3
            D1I= P12(I,1)-P34(I,1)
            F5X(0,I,I,I)=-D1I     *THR
            F5X(1,I,I,I)=-P12(I,3)*THR
            F5X(2,I,I,I)= P34(I,3)*THR
            DO 297 M=0,9
  297       F5X(3+M,I,I,I)=-R3(M,I,I,I)
            DO 299 J=1,3
               IF(I.NE.J) THEN
                  K=6-I-J
                  F5X(0,K,J,I)= ZER
                  F5X(1,K,J,I)= ZER
                  F5X(2,K,J,I)= ZER
                  D1J= P12(J,1)-P34(J,1)
                  F5X(0,J,I,I)=-D1J
                  F5X(1,J,I,I)=-P12(J,3)
                  F5X(2,J,I,I)= P34(J,3)
                  F5X(0,I,J,I)=-D1J
                  F5X(1,I,J,I)=-P12(J,3)
                  F5X(2,I,J,I)= P34(J,3)
                  F5X(0,I,I,J)=-D1J
                  F5X(1,I,I,J)=-P12(J,3)
                  F5X(2,I,I,J)= P34(J,3)
                  DO 298 M=0,9
                     F5X(3+M,K,J,I)=-R3(M,3,2,1)
                     F5X(3+M,J,I,I)=-R3(M,J,I,I)
                     F5X(3+M,I,J,I)=-R3(M,J,I,I)
                     F5X(3+M,I,I,J)=-R3(M,J,I,I)
  298             CONTINUE
               ENDIF
  299    CONTINUE
C
         M56=(JTYPE-5)*3
         DO 300 L=0,M56
            DO 300 M=1,3
               DO 300 N=0,17
  300    A56(N,M,L)= ZER
         DO 309 K=1,3
            DO 309 J=1,3
               DO 309 I=1,3
                  MM= 8
                  IF(I.EQ.J) MM=11
                  T01= P12(I,3)*P1D(J,K)
                  T02= P12(I,3)*P1P(J,K)
                  T03= P1D(J,I)*P34(K,3)
                  T04= P24(J,I)*P34(K,3)
                  T( 0)=-T01
                  T( 1)= T01-T02
                  T( 2)= T02
                  T( 3)= P12(I,3)*P24(J,K)
                  T( 4)=-T03
                  T( 5)=(P1D(I,J)-P1P(I,J))*P34(K,3)+T03
                  T( 6)=(P1P(I,J)+P1P(I,J))*P34(K,3)
                  T( 7)= T04
                  T( 8)=-P24(I,J)*P34(K,3)-T04
                  T( 9)=-P12(K,1)+P34(K,1)
                  T(10)=-P12(K,3)
                  T(11)= P34(K,3)
                  T(12)=-P12(I,3)
                  T(13)= P12(J,3)
                  U( 0)= F5X(0,I,J,K)
                  U( 1)= F5X(1,I,J,K)
                  U( 2)= F5X(2,I,J,K)
                  U( 3)=-T1R(1,J,I,K)
                  U( 4)= T1R(1,I,J,K)+T1R(1,J,I,K)-T1R(2,J,I,K)
                  U( 5)= T1R(2,I,J,K)+T1R(2,J,I,K)-T1R(3,J,I,K)
                  U( 6)= T1R(3,I,J,K)+T1R(3,J,I,K)
                  U( 7)= T1R(4,J,I,K)
                  U( 8)= T1R(4,I,J,K)+T1R(4,J,I,K)-T1R(5,J,I,K)
                  U( 9)= T1R(5,I,J,K)+T1R(5,J,I,K)
                  U(10)= T1R(6,J,I,K)
                  U(11)=-T1R(6,I,J,K)-T1R(6,J,I,K)
                  DO 301 M=1,6
  301             U(11+M)= T3R(M,K,I,J)
                  DO 309 L=0,M56
                     GEIJKL= E(I,J,K,L)
C V(I,J,K,0)
                     DO 303 M=0,MM
  303                A56( M,1,L)= A56( M,1,L)+T(M)*GEIJKL
                     IF(J.EQ.K) A56(12,1,L)= A56(12,1,L)+T(12)*GEIJKL
                     IF(I.EQ.K) A56(13,1,L)= A56(13,1,L)+T(13)*GEIJKL
                     DO 305 M=0,17
  305                A56( M,2,L)= A56( M,2,L)+U( M)*GEIJKL
                     DO 307 M=0,9
  307                A56( M,3,L)= A56( M,3,L)+F5X(3+M,I,J,K)*GEIJKL
  309    CONTINUE
         IF(JTYPE.LT.6) GO TO 490
C
         DO 312 L=0,3
            DO 310 M=0,2
               A60(M,L)= ZER
               A61(M,L)= ZER
  310       CONTINUE
            DO 312 I=1,3
               D1I= P12(I,1)-P34(I,1)
C P1(I)
               A60(0,L)= A60(0,L)+D1I     *E(I,0,0,L)
               A60(1,L)= A60(1,L)+P12(I,3)*E(I,0,0,L)
               A60(2,L)= A60(2,L)+P34(I,3)*E(I,0,0,L)
C V(0,0,0,L)
               A61(0,L)= A61(0,L)+D1I     *E(L,0,0,I)
               A61(1,L)= A61(1,L)+P12(I,3)*E(L,0,0,I)
               A61(2,L)= A61(2,L)+P34(I,3)*E(L,0,0,I)
  312    CONTINUE
         DO 314 N=1,6
            A62(  N)= ZER
            A62(8+N)= ZER
            A63(  N)= ZER
            A63(8+N)= ZER
  314    CONTINUE
         SEL00L= E(1,0,0,1)+E(2,0,0,2)+E(3,0,0,3)
         DO 318 L=1,3
            DO 318 I=1,3
               GEI00L= E(I,0,0,L)
C V(I,0,0,L)
               A62( 1)= A62( 1)+ P1D(I,L)*GEI00L
               A62( 2)= A62( 2)+ P1P(I,L)*GEI00L
               A62( 3)= A62( 3)+ P24(I,L)*GEI00L
               A62( 4)= A62( 4)+ P3D(L,I)*GEI00L
               A62( 5)= A62( 5)- P3P(I,L)*GEI00L
               DO 316 N=1,6
  316          A62(8+N)= A62(8+N)+R2(N,I,L)*GEI00L
  318    CONTINUE
         SE0L0L= E(0,1,0,1)+E(0,2,0,2)+E(0,3,0,3)
         DO 322 L=1,3
            DO 322 J=1,3
               GE0J0L= E(0,J,0,L)
C V(0,J,0,L)
               A63( 1)= A63( 1)+ P1D(J,L)*GE0J0L
               A63( 2)= A63( 2)+ P1P(J,L)*GE0J0L
               A63( 3)= A63( 3)+ P24(J,L)*GE0J0L
               A63( 4)= A63( 4)+ P3D(L,J)*GE0J0L
               A63( 5)= A63( 5)- P3P(J,L)*GE0J0L
               DO 320 N=1,6
  320          A63(8+N)= A63(8+N)+R2(N,J,L)*GE0J0L
  322    CONTINUE
         DO 326 L=1,3
            DO 324 M=0,2
               A64(M,L)= ZER
               A65(M,L)= ZER
  324       CONTINUE
            DO 326 J=1,3
               D1J= P12(J,1)-P34(J,1)
               A64(0,L)= A64(0,L)+D1J     *E(0,J,0,L)
               A64(1,L)= A64(1,L)+P12(J,3)*E(0,J,0,L)
               A64(2,L)= A64(2,L)+P34(J,3)*E(0,J,0,L)
               A65(0,L)= A65(0,L)+D1J     *E(0,L,0,J)
               A65(1,L)= A65(1,L)+P12(J,3)*E(0,L,0,J)
               A65(2,L)= A65(2,L)+P34(J,3)*E(0,L,0,J)
  326    CONTINUE
         DO 328 L=1,3
            DO 328 M=0,2
               A69(M,L)= ZER
               A70(M,L)= ZER
  328    CONTINUE
         DO 329 L=1,3
            DO 329 N=1,6
               A66(   N,L)= ZER
               A66(10+N,L)= ZER
               A67(   N,L)= ZER
               A67(10+N,L)= ZER
               A68(   N,L)= ZER
               A68(10+N,L)= ZER
               A71(   N,L)= ZER
               A71(10+N,L)= ZER
               A72(   N,L)= ZER
               A72(10+N,L)= ZER
  329    CONTINUE
         DO 330 L=0,3
            DO 330 N=1,6
               A73(   N,L)= ZER
               A73(10+N,L)= ZER
  330    CONTINUE
         DO 338 L=1,3
            A66( 1,L)= E(1,1,0,L)+E(2,2,0,L)+E(3,3,0,L)
            A67( 1,L)= E(1,L,0,1)+E(2,L,0,2)+E(3,L,0,3)
            A68( 1,L)= E(L,1,0,1)+E(L,2,0,2)+E(L,3,0,3)
            DO 338 J=1,3
               DO 338 I=1,3
                  GEIJ0L= E(I,J,0,L)
C V(I,J,0,0)
                  A66( 2,L)= A66( 2,L)+ P1D(J,I)*GEIJ0L
                  A66( 3,L)= A66( 3,L)+ P1P(I,J)*GEIJ0L
                  A66( 4,L)= A66( 4,L)+ P24(J,I)*GEIJ0L
                  A66( 5,L)= A66( 5,L)+ P1D(I,J)*GEIJ0L
                  A66( 6,L)= A66( 6,L)- P24(I,J)*GEIJ0L
                  DO 332 N=1,6
  332             A66(10+N,L)= A66(10+N,L)+R2(N,I,J)*GEIJ0L
C V(I,0,0,L)
                  A67( 2,J)= A67( 2,J)+ P1D(I,L)*GEIJ0L
                  A67( 3,J)= A67( 3,J)+ P1P(I,L)*GEIJ0L
                  A67( 4,J)= A67( 4,J)+ P24(I,L)*GEIJ0L
                  A67( 5,J)= A67( 5,J)+ P3D(L,I)*GEIJ0L
                  A67( 6,J)= A67( 6,J)- P3P(I,L)*GEIJ0L
                  DO 334 N=1,6
  334             A67(10+N,J)= A67(10+N,J)+R2(N,I,L)*GEIJ0L
C V(0,J,0,L)
                  A68( 2,I)= A68( 2,I)+ P1D(J,L)*GEIJ0L
                  A68( 3,I)= A68( 3,I)+ P1P(J,L)*GEIJ0L
                  A68( 4,I)= A68( 4,I)+ P24(J,L)*GEIJ0L
                  A68( 5,I)= A68( 5,I)+ P3D(L,J)*GEIJ0L
                  A68( 6,I)= A68( 6,I)- P3P(J,L)*GEIJ0L
                  DO 336 N=1,6
  336             A68(10+N,I)= A68(10+N,I)+R2(N,J,L)*GEIJ0L
  338    CONTINUE
         DO 342 L=1,3
            D1L= P12(L,1)-P34(L,1)
            DO 342 K=1,3
               D1K= P12(K,1)-P34(K,1)
C V(0,0,L,0)
               A69(0,L)= A69(0,L)+D1K     *E(0,0,K,L)
               A69(1,L)= A69(1,L)+P12(K,3)*E(0,0,K,L)
               A69(2,L)= A69(2,L)+P34(K,3)*E(0,0,K,L)
C V(0,0,0,L)
               A70(0,K)= A70(0,K)+D1L     *E(0,0,K,L)
               A70(1,K)= A70(1,K)+P12(L,3)*E(0,0,K,L)
               A70(2,K)= A70(2,K)+P34(L,3)*E(0,0,K,L)
  342    CONTINUE
         DO 348 L=1,3
            A71( 1,L)= E(1,0,1,L)+E(2,0,2,L)+E(3,0,3,L)
            A72( 1,L)= E(1,0,L,1)+E(2,0,L,2)+E(3,0,L,3)
            DO 348 K=1,3
               DO 348 I=1,3
                  GEI0KL= E(I,0,K,L)
C V(I,0,K,0)
                  A71( 2,L)= A71( 2,L)+ P1D(I,K)*GEI0KL
                  A71( 3,L)= A71( 3,L)+ P1P(I,K)*GEI0KL
                  A71( 4,L)= A71( 4,L)+ P24(I,K)*GEI0KL
                  A71( 5,L)= A71( 5,L)+ P3D(K,I)*GEI0KL
                  A71( 6,L)= A71( 6,L)- P3P(I,K)*GEI0KL
                  DO 344 N=1,6
  344             A71(10+N,L)= A71(10+N,L)+R2(N,I,K)*GEI0KL
C V(I,0,0,L)
                  A72( 2,K)= A72( 2,K)+ P1D(I,L)*GEI0KL
                  A72( 3,K)= A72( 3,K)+ P1P(I,L)*GEI0KL
                  A72( 4,K)= A72( 4,K)+ P24(I,L)*GEI0KL
                  A72( 5,K)= A72( 5,K)+ P3D(L,I)*GEI0KL
                  A72( 6,K)= A72( 6,K)- P3P(I,L)*GEI0KL
                  DO 346 N=1,6
  346             A72(10+N,K)= A72(10+N,K)+R2(N,I,L)*GEI0KL
  348    CONTINUE
            A73( 1,0)= E(0,0,1,1)+E(0,0,2,2)+E(0,0,3,3)
         DO 353 L=1,3
            A73( 1,L)= E(L,0,1,1)+E(L,0,2,2)+E(L,0,3,3)
            DO 353 K=1,3
               DO 353 I=0,3
                  GEI0KL= E(I,0,K,L)
                  A73( 2,I)= A73( 2,I)+ P3D(L,K)*GEI0KL
                  A73( 3,I)= A73( 3,I)+ P24(K,L)*GEI0KL
                  A73( 4,I)= A73( 4,I)- P3P(K,L)*GEI0KL
                  A73( 5,I)= A73( 5,I)+ P3D(K,L)*GEI0KL
                  A73( 6,I)= A73( 6,I)+ P24(L,K)*GEI0KL
                  DO 351 N=1,6
  351             A73(10+N,I)= A73(10+N,I)+R2(N,K,L)*GEI0KL
  353    CONTINUE
C
         DO 359 N=1,3
            DO 357 M=1,5
               A74(  M,N)= ZER
               A74(8+M,N)= ZER
               A75(  M,N)= ZER
               A75(8+M,N)= ZER
               A76(  M,N)= ZER
               A76(8+M,N)= ZER
  357       CONTINUE
            A75(14,N)= ZER
            A76(14,N)= ZER
  359    CONTINUE
         DO 369 L=1,3
            SE0KKL(L)= E(0,1,1,L)+E(0,2,2,L)+E(0,3,3,L)
            SE0LKL(L)= E(0,1,L,1)+E(0,2,L,2)+E(0,3,L,3)
            SE0JLL(L)= E(0,L,1,1)+E(0,L,2,2)+E(0,L,3,3)
            DO 367 K=1,3
               DO 367 J=1,3
                  GE0JKL= E(0,J,K,L)
C V(0,J,K,0)
                  DO 361 N=1,5
                     A74(  N,L)= A74(  N,L)+Q1(N,J,K)*GE0JKL
                     A74(8+N,L)= A74(8+N,L)+R2(N,J,K)*GE0JKL
  361             CONTINUE
C V(0,J,0,L)
                  A75( 1,K)= A75( 1,K)+ P1D(J,L)*GE0JKL
                  A75( 2,K)= A75( 2,K)+ P1P(J,L)*GE0JKL
                  A75( 3,K)= A75( 3,K)+ P24(J,L)*GE0JKL
                  A75( 4,K)= A75( 4,K)+ P3D(L,J)*GE0JKL
                  A75( 5,K)= A75( 5,K)- P3P(J,L)*GE0JKL
                  DO 363 N=1,6
  363             A75(8+N,K)= A75(8+N,K)+R2(N,J,L)*GE0JKL
C V(0,0,K,L)
                  A76( 1,J)= A76( 1,J)+ P3D(L,K)*GE0JKL
                  A76( 2,J)= A76( 2,J)+ P24(K,L)*GE0JKL
                  A76( 3,J)= A76( 3,J)- P3P(K,L)*GE0JKL
                  A76( 4,J)= A76( 4,J)+ P3D(K,L)*GE0JKL
                  A76( 5,J)= A76( 5,J)+ P24(L,K)*GE0JKL
                  DO 365 N=1,6
  365             A76(8+N,J)= A76(8+N,J)+R2(N,K,L)*GE0JKL
  367       CONTINUE
            A74(14,L)= A74( 5,L)
  369    CONTINUE
C
         DO 372 M=0,14
  372    S4(M)= ZER
         DO 378 L=1,3
            DO 378 K=1,3
               DO 374 N=1,6
  374          U(N)= R2(N,K,L)
               DO 378 J=1,3
                  DO 378 I=1,3
                     T( 0)= R2(1,I,J)*U(1)
                     T( 1)= R2(1,I,J)*U(2)+R2(2,I,J)*U(1)
                     T( 2)= R2(1,I,J)*U(3)+R2(2,I,J)*U(2)+R2(3,I,J)*U(1)
                     T( 3)= R2(2,I,J)*U(3)+R2(3,I,J)*U(2)
                     T( 4)= R2(3,I,J)*U(3)
                     T( 5)=-R2(1,I,J)*U(4)-R2(4,I,J)*U(1)
                     T( 6)=-R2(1,I,J)*U(5)-R2(2,I,J)*U(4)
     2                     -R2(4,I,J)*U(2)-R2(5,I,J)*U(1)
                     T( 7)=-R2(2,I,J)*U(5)-R2(3,I,J)*U(4)
     2                     -R2(4,I,J)*U(3)-R2(5,I,J)*U(2)
                     T( 8)=-R2(3,I,J)*U(5)-R2(5,I,J)*U(3)
                     T( 9)= R2(1,I,J)*U(6)+R2(4,I,J)*U(4)+R2(6,I,J)*U(1)
                     T(10)= R2(2,I,J)*U(6)+R2(4,I,J)*U(5)
     2                     +R2(5,I,J)*U(4)+R2(6,I,J)*U(2)
                     T(11)= R2(3,I,J)*U(6)+R2(5,I,J)*U(5)+R2(6,I,J)*U(3)
                     T(12)=-R2(4,I,J)*U(6)-R2(6,I,J)*U(4)
                     T(13)=-R2(5,I,J)*U(6)-R2(6,I,J)*U(5)
                     T(14)= R2(6,I,J)*U(6)
                     GEIJKL= E(I,J,K,L)
                     DO 376 M=0,14
  376                S4(M)= S4(M)+T(M)*GEIJKL
  378    CONTINUE
C
         DO 386 I=1,3
            F6X(0,I,I,I,I)= THR
            DO 382 M=1,6
  382       F6X(M,I,I,I,I)= R2(M,I,I)*SIX
            DO 386 J=1,3
               IF(I.NE.J) THEN
                  K=6-I-J
                  F6X(0,J,I,I,I)= ZER
                  F6X(0,I,J,I,I)= ZER
                  F6X(0,I,I,J,I)= ZER
                  F6X(0,I,I,I,J)= ZER
                  F6X(0,J,J,I,I)= ONE
                  F6X(0,J,I,J,I)= ONE
                  F6X(0,I,J,J,I)= ONE
                  F6X(0,J,I,I,J)= ONE
                  F6X(0,I,J,I,J)= ONE
                  F6X(0,I,I,J,J)= ONE
                  F6X(0,K,J,I,I)= ZER
                  F6X(0,K,I,J,I)= ZER
                  F6X(0,I,K,J,I)= ZER
                  F6X(0,K,I,I,J)= ZER
                  F6X(0,I,K,I,J)= ZER
                  F6X(0,I,I,K,J)= ZER
                  DO 384 M=1,6
                     EQ0= R2(M,J,I)*THR
                     F6X(M,J,I,I,I)= EQ0
                     F6X(M,I,J,I,I)= EQ0
                     F6X(M,I,I,J,I)= EQ0
                     F6X(M,I,I,I,J)= EQ0
                     EQ1= R2(M,I,I)+R2(M,J,J)
                     F6X(M,J,J,I,I)= EQ1
                     F6X(M,J,I,J,I)= EQ1
                     F6X(M,I,J,J,I)= EQ1
                     F6X(M,J,I,I,J)= EQ1
                     F6X(M,I,J,I,J)= EQ1
                     F6X(M,I,I,J,J)= EQ1
                     F6X(M,K,J,I,I)= R2(M,K,J)
                     F6X(M,K,I,J,I)= R2(M,K,J)
                     F6X(M,I,K,J,I)= R2(M,K,J)
                     F6X(M,K,I,I,J)= R2(M,K,J)
                     F6X(M,I,K,I,J)= R2(M,K,J)
                     F6X(M,I,I,K,J)= R2(M,K,J)
  384             CONTINUE
               ENDIF
  386    CONTINUE
C
         DO 395 N=1,2
            DO 391 M=0,20
  391       A82( M,N)= ZER
            DO 392 M=0,14
  392       A83( M,N)= ZER
            DO 393 M=0,9
               A84(M,N)= ZER
               A85(M,N)= ZER
  393       CONTINUE
  395    CONTINUE
         DO 397 M=0,6
  397    A86( M)= ZER
C
         DO 399 K=1,3
            DO 399 J=1,3
               DO 399 I=1,3
                  T9A( 0,I,J,K)= P12(I,3)*P3D(J,K)
                  T9A( 1,I,J,K)= P12(I,3)*P24(K,J)
                  T9A( 2,I,J,K)=-P12(I,3)*P3P(J,K)
                  T9A( 3,I,J,K)=-P12(I,3)*P3D(K,J)
                  T9A( 4,I,J,K)=-P12(I,3)*P24(J,K)
                  T9A( 5,I,J,K)= P12(I,3)*P3P(J,K)
                  T9A( 6,I,J,K)=-P34(J,3)*P3D(K,I)
                  T9A( 7,I,J,K)=-P34(J,3)*P24(I,K)
                  T9A( 8,I,J,K)= P34(J,3)*P3P(I,K)
                  T9A( 9,I,J,K)= P12(I,1)-P34(I,1)
                  T9A(10,I,J,K)=-P12(I,3)
                  T9A(11,I,J,K)=-P34(I,3)
                  T9A(12,I,J,K)= P34(K,3)
                  T9A(13,I,J,K)=-P34(J,3)
                  T9B( 0,I,J,K)= F5X(0,I,J,K)
                  T9B( 1,I,J,K)= F5X(1,I,J,K)
                  T9B( 2,I,J,K)= F5X(2,I,J,K)
                  DO 398 M=1,6
                     T9B(M+ 2,I,J,K)= T3R(M,K,I,J)
                     T9B(M+ 8,I,J,K)= T1R(M,I,J,K)
                     T9B(M+14,I,J,K)= T3R(M,J,I,K)
  398             CONTINUE
  399    CONTINUE
C
         A80(0,1)= ZER
         A80(1,1)= ZER
         A80(0,2)= ZER
         A80(1,2)= ZER
         DO 408 L=1,3
            DO 408 K=1,3
               A80(1,2)= A80(1,2)+E(K,K,L,L)
               DO 408 J=1,3
                  DO 408 I=1,3
                     GEIJKL= E(I,J,K,L)
                     A80(0,1)= A80(0,1)+P1P(I,J)*P3P(K,L)*GEIJKL
                     T02=-P12(J,3)*GEIJKL
                     T03= P12(I,3)*GEIJKL
                     DO 401 M=0,20
                        A82( M,1)= A82( M,1)+T9B( M,I,K,L)*T02
                        A82( M,2)= A82( M,2)+T9B( M,J,K,L)*T03
  401                CONTINUE
                     T04=-P34(L,3)*GEIJKL
                     T05= P34(K,3)*GEIJKL
                     DO 402 M=0,2
                        A83( M,1)= A83( M,1)+F5X( M,I,J,K)*T04
                        A83( M,2)= A83( M,2)+F5X( M,I,J,L)*T05
  402                CONTINUE
                     T00=-P1P(I,J)*GEIJKL
                     T01=-P3P(K,L)*GEIJKL
                     DO 403 M=1,6
                        A83(2+M,1)= A83(2+M,1)+R2(M,K,L)*T00
                        A83(2+M,2)= A83(2+M,2)+R2(M,I,J)*T01
  403                CONTINUE
                     IF(I.EQ.J) THEN
                        A80(1,1)= A80(1,1)+T01
                        DO 404 M=1,6
  404                   A83(8+M,1)= A83(8+M,1)+R2(M,K,L)*GEIJKL
                     ENDIF
                     IF(K.EQ.L) THEN
                        A80(0,2)= A80(0,2)+T00
                        DO 405 M=1,6
  405                   A83(8+M,2)= A83(8+M,2)+R2(M,I,J)*GEIJKL
                     ENDIF
                     DO 406 M=0,9
                        A84(M,1)= A84(M,1)+F5X(3+M,I,K,L)*T02
                        A84(M,2)= A84(M,2)+F5X(3+M,J,K,L)*T03
                        A85(M,1)= A85(M,1)+F5X(3+M,I,J,K)*T04
                        A85(M,2)= A85(M,2)+F5X(3+M,I,J,L)*T05
  406                CONTINUE
                     DO 407 M=0,6
  407                A86(M)= A86(M)+F6X(M,I,J,K,L)*GEIJKL
  408    CONTINUE
C
         DO 409 L=0,3
            DO 409 M=1,3
               DO 409 N=0,20
                  A77(N,M,L)= ZER
                  A78(N,M,L)= ZER
                  A79(N,M,L)= ZER
  409    CONTINUE
         DO 419 L=1,3
            DO 419 J=1,3
               DO 419 I=1,3
                  MM= 8
                  IF(I.EQ.J) MM=11
                  T( 0)=-P12(I,3)*P1D(J,L)
                  T( 1)=-P12(I,3)*P1P(J,L)
                  T( 2)= P12(I,3)*P24(J,L)
                  T( 3)=-P12(I,3)*P3D(L,J)
                  T( 4)=-P12(I,3)*P24(J,L)
                  T( 5)= P12(I,3)*P3P(J,L)
                  T( 6)= P12(J,3)*P3D(L,I)
                  T( 7)= P12(J,3)*P24(I,L)
                  T( 8)=-P12(J,3)*P3P(I,L)
                  T( 9)= P12(L,1)-P34(L,1)
                  T(10)= P12(L,3)
                  T(11)=-P34(L,3)
                  T(12)= P12(I,3)
                  T(13)=-P12(J,3)
                  U( 0)= F5X(0,I,J,L)
                  U( 1)= F5X(1,I,J,L)
                  U( 2)= F5X(2,I,J,L)
                  DO 411 M=1,6
                     U(M+ 2)= T1R(M,I,J,L)
                     U(M+ 8)=-T1R(M,J,I,L)
                     U(M+14)=-T3R(M,L,I,J)
  411             CONTINUE
                  DO 419 K=0,3
                     GEIJKL= E(I,J,K,L)
C V(I,J,0,L)
                     DO 413 M=0,MM
  413                A77( M,1,K)= A77( M,1,K)+T( M)*GEIJKL
                     IF(J.EQ.L) A77(12,1,K)= A77(12,1,K)+T(12)*GEIJKL
                     IF(I.EQ.L) A77(13,1,K)= A77(13,1,K)+T(13)*GEIJKL
                     DO 415 M=0,20
  415                A77( M,2,K)= A77( M,2,K)+U( M)*GEIJKL
                     DO 417 M=0,9
  417                A77( M,3,K)= A77( M,3,K)+F5X(3+M,I,J,L)*GEIJKL
  419    CONTINUE
         DO 429 L=1,3
            DO 429 K=1,3
               MM= 8
               IF(K.EQ.L) MM=11
                  A78(12,1,0)= A78(12,1,0)+T9A(12,K,K,L)*E(K,0,K,L)
                  A78(13,1,0)= A78(13,1,0)+T9A(13,L,K,L)*E(L,0,K,L)
               DO 429 I=1,3
                  A78(12,1,I)= A78(12,1,I)+T9A(12,K,K,L)*E(K,I,K,L)
                  A78(13,1,I)= A78(13,1,I)+T9A(13,L,K,L)*E(L,I,K,L)
                  DO 429 J=0,3
                     GEIJKL= E(I,J,K,L)
C V(I,0,K,L)
                     DO 425 M=0,MM
  425                A78( M,1,J)= A78( M,1,J)+T9A( M,I,K,L)*GEIJKL
                     DO 427 M=0,20
  427                A78( M,2,J)= A78( M,2,J)+T9B( M,I,K,L)*GEIJKL
                     DO 428 M=0,9
  428                A78( M,3,J)= A78( M,3,J)+F5X(3+M,I,K,L)*GEIJKL
  429    CONTINUE
         DO 439 L=1,3
            DO 439 K=1,3
               MM=11
               IF(K.EQ.L) MM=14
                  A79(15,1,0)= A79(15,1,0)+T9A(12,K,K,L)*E(0,K,K,L)
                  A79(16,1,0)= A79(16,1,0)+T9A(13,L,K,L)*E(0,L,K,L)
               DO 439 J=1,3
                  A79(15,1,J)= A79(15,1,J)+T9A(12,K,K,L)*E(J,K,K,L)
                  A79(16,1,J)= A79(16,1,J)+T9A(13,L,K,L)*E(J,L,K,L)
                  T( 0)= T9A( 0,J,L,K)
                  T( 1)= T9A( 1,J,L,K)
                  T( 2)= T9A( 2,J,L,K)
                  DO 431 M=0,13
  431             T(3+M)=T9A( M,J,K,L)
                  DO 439 I=0,3
                     GEIJKL= E(I,J,K,L)
C V(0,J,K,L)
                     IF(K.EQ.L) A79(17,1,I)=A79(17,1,I)-P12(J,3)*GEIJKL
                     DO 435 M=0,MM
  435                A79( M,1,I)= A79( M,1,I)+T( M)*GEIJKL
                     DO 437 M=0,20
  437                A79( M,2,I)= A79( M,2,I)+T9B( M,J,K,L)*GEIJKL
                     DO 438 M=0,9
  438                A79( M,3,I)= A79( M,3,I)+F5X(3+M,J,K,L)*GEIJKL
  439    CONTINUE
         DO 441 N=1,3
            DO 441 M=0,13
  441    A81( M,N)= ZER
         DO 449 L=1,3
            DO 449 K=1,3
               A81(12,3)= A81(12,3)+E(K,K,L,L)
               MM= 8
               IF(K.EQ.L) MM=11
               T( 0)= P3D(L,K)
               T( 1)= P24(K,L)
               T( 2)=-P3P(K,L)
               T( 3)= P3D(K,L)
               T( 4)= P24(L,K)
               T( 6)= P3P(K,L)
               T( 7)= T( 0)
               T( 8)= T( 1)
               T( 9)= T( 2)
               T(10)= T( 3)
               T(11)= T( 4)
               DO 449 J=1,3
                  DO 449 I=1,3
                     GEIJKL= E(I,J,K,L)
                     T01= P12(J,3)*GEIJKL
                     T02= P12(I,3)*GEIJKL
                     DO 443 M=0,MM
                        A81( M,1)= A81( M,1)-T9A(M,I,K,L)*T01
                        A81( M,2)= A81( M,2)+T9A(M,J,K,L)*T02
  443                CONTINUE
                     IF(I.EQ.K) A81(12,1)= A81(12,1)-T9A(12,I,K,L)*T01
                     IF(I.EQ.L) A81(13,1)= A81(13,1)-T9A(13,I,K,L)*T01
                     IF(J.EQ.K) A81(12,2)= A81(12,2)+T9A(12,J,K,L)*T02
                     IF(J.EQ.L) A81(13,2)= A81(13,2)+T9A(13,J,K,L)*T02
                     TMP= P1P(I,J)*GEIJKL
                     DO 445 M=0,4
  445                A81( M,3) = A81( M,3)+T( M)*TMP
                     IF(K.EQ.L)A81( 5,3)=A81( 5,3)+TMP
                     IF(I.EQ.J) THEN
                        DO 447 M=6,11
  447                   A81( M,3) = A81( M,3)+T( M)*GEIJKL
                     ENDIF
  449    CONTINUE
      ENDIF
  490 CONTINUE
C
      LK= LK+1
      IF(LK.EQ.1) THEN
C
CJMS  THE DEPENDENCY ON X01 AND X02 IS TREATED IN THIS BLOCK.
CJMS  IT IS ONLY CALCULATED ONCE PER OUTER LOOP (LOOPS 940 AND 930)
C
         IF(JTYPE.NE.4) THEN
            DO 510 J=0,3
               A31(3,J)= A31(0,J)+ A31(1,J)*Y02
               A32(3,J)= A32(0,J)+ A32(1,J)*Y02
  510       CONTINUE
            T01= SEJJ00*H21
            A33( 6)=-AEIJ00*Y1Y+T01
            A33( 7)= A33( 1)-(A33( 2)- A33( 3)*Y01)*Y01+T01
            A33( 8)= A33( 4)- A33( 5)*Y01
            A33(15)= A33( 9)+(A33(10)+ A33(11)*Y02)*Y02
            A33(16)= A33(12)+ A33(13)*Y02
         ENDIF
         IF(JTYPE.LT.4) GO TO 590
C
         DO 520 K=0,3
            A41(3,K)= A41(0,K)+ A41(1,K)*Y02
            A42(3,K)= A42(0,K)+ A42(1,K)*Y02
  520    CONTINUE
         A43( 6)= SEK0K0*H21
         A43( 7)=(A43( 1)+ A43( 2)*Y02)*Y02
         A43( 8)= A43( 3)*Y02
         A43(15)= A43( 9)+(A43(10)+ A43(11)*Y02)*Y02
         A43(16)= A43(12)+ A43(13)*Y02
         IF(JTYPE.LT.5) GO TO 590
C
         DO 530 N=1,3
            A51(3,N)= A51(0,N)+ A51(1,N)*Y02
            A52(3,N)= A52(0,N)+ A52(1,N)*Y02
C V(I,J,0,0)
            T01= SEJJK0(N)*H21
            A54( 6,N)=-A54(11,N)*Y1Y+T01
            A54( 7,N)= A54( 1,N)-(A54( 2,N)- A54( 3,N)*Y01)*Y01+T01
            A54( 8,N)= A54( 4,N)- A54( 5,N)*Y01
            A54(15,N)= A54( 9,N)+(A54(10,N)+ A54(11,N)*Y02)*Y02
            A54(16,N)= A54(12,N)+ A54(13,N)*Y02
C V(I,0,K,0)
            A55( 6,N)= SEKJK0(N)*H21
            A55( 7,N)=(A55( 1,N)+ A55( 2,N)*Y02)*Y02
            A55( 8,N)= A55( 3,N)*Y02
            A55(15,N)= A55( 9,N)+(A55(10,N)+ A55(11,N)*Y02)*Y02
            A55(16,N)= A55(12,N)+ A55(13,N)*Y02
  530    CONTINUE
         DO 535 M=0,3
C V(0,J,K,0)
            A53( 6,M)= SEIKK0(M)*H21
            A53( 7,M)=(A53( 1,M)+ A53( 2,M)*Y02)*Y01
            A53( 8,M)=-A53( 3,M)*Y01
            A53(15,M)= A53( 9,M)+(A53(10,M)+ A53(11,M)*Y02)*Y02
            A53(16,M)= A53(12,M)+ A53(13,M)*Y02
  535    CONTINUE
C
         DO 540 M=0,M56
            T01        = A56( 9,1,M)+ A56(13,1,M)
     2                 +(A56(10,1,M)- A56(13,1,M)+ A56(12,1,M))*Y02
            T02        = A56(11,1,M)*H21
            A56(18,1,M)=(A56( 0,1,M)+(A56( 1,1,M)+ A56( 2,1,M)*Y02)*Y02)
     2                  *Y02- T01*H21
            A56(19,1,M)= A56( 3,1,M)*Y1Y- T02
            A56(20,1,M)= A56( 4,1,M)+(A56( 5,1,M)+ A56( 6,1,M)*Y02)*Y02
     2                  +T02
            A56(21,1,M)= A56( 7,1,M)+ A56( 8,1,M)*Y02
            A56(18,2,M)=(A56( 0,2,M)+ A56( 1,2,M)*Y02)*H21
     2                 -(A56( 3,2,M)+(A56( 4,2,M)+(A56( 5,2,M)
     3                              + A56( 6,2,M)*Y02)*Y02)*Y02)
            A56(19,2,M)= A56( 2,2,M)*H21
     2                 -(A56( 7,2,M)-(A56( 8,2,M)+ A56( 9,2,M)*Y02)*Y02)
            A56(20,2,M)= A56(10,2,M)+ A56(11,2,M)*Y02
            A56(21,2,M)= A56(12,2,M)+(A56(13,2,M)+ A56(14,2,M)*Y02)*Y02
            A56(22,2,M)= A56(15,2,M)+ A56(16,2,M)*Y02
            A56(18,3,M)= A56( 0,3,M)+(A56( 1,3,M)
     2                 +(A56( 2,3,M)+ A56( 3,3,M)*Y02)*Y02)*Y02
            A56(19,3,M)= A56( 4,3,M)+(A56( 5,3,M)+ A56( 6,3,M)*Y02)*Y02
            A56(20,3,M)= A56( 7,3,M)+ A56( 8,3,M)*Y02
            A56(21,3,M)= A56( 9,3,M)
  540    CONTINUE
         IF(JTYPE.LT.6) GO TO 590
C
         A62( 6)= SEL00L*H21
         A62( 7)=(A62( 1)+ A62( 2)*Y02)*Y02
         A62( 8)= A62( 4)+ A62( 3)*Y02
         A62(15)= A62( 9)+(A62(10)+A62(11)*Y02)*Y02
         A62(16)= A62(12)+ A62(13)*Y02
         A63( 6)= SE0L0L*H21
         A63( 7)=(A63( 1)+ A63( 2)*Y02)*Y01
         A63( 8)= A63( 4)+ A63( 3)*Y02
         A63(15)= A63( 9)+(A63(10)+A63(11)*Y02)*Y02
         A63(16)= A63(12)+ A63(13)*Y02
         DO 550 N=1,3
            A64(3,N)= A64(0,N)+ A64(1,N)*Y02
            A65(3,N)= A65(0,N)+ A65(1,N)*Y02
            T02      = A66( 3,N)*Y02
            A66( 7,N)= A66( 1,N)*H21-(A66( 2,N)+T02)
     2               +(A66( 2,N)+ A66( 5,N)+T02+T02)*Y02
            A66( 8,N)= A66( 4,N)*Y01+ A66( 6,N)*Y02
            A66( 9,N)= A66( 1,N)*H21-T02*Y01
            A66(17,N)= A66(11,N)+(A66(12,N)+ A66(13,N)*Y02)*Y02
            A66(18,N)= A66(14,N)+ A66(15,N)*Y02
            T02      = A67( 4,N)*Y02
            A67( 7,N)= A67( 1,N)*H21
            A67( 8,N)=(A67( 2,N)+ A67( 3,N)*Y02)*Y02
            A67( 9,N)= A67( 5,N)+T02
            A67(10,N)= T02
            A67(17,N)= A67(11,N)+(A67(12,N)+ A67(13,N)*Y02)*Y02
            A67(18,N)= A67(14,N)+ A67(15,N)*Y02
            T02      = A68( 4,N)*Y02
            A68( 7,N)= A68( 1,N)*H21
            A68( 8,N)=(A68( 2,N)+ A68( 3,N)*Y02)*Y01
            A68( 9,N)= A68( 5,N)+T02
            A68(10,N)= A68( 4,N)-T02
            A68(17,N)= A68(11,N)+(A68(12,N)+ A68(13,N)*Y02)*Y02
            A68(18,N)= A68(14,N)+ A68(15,N)*Y02
            A69(3,N)= A69(0,N)+ A69(1,N)*Y02
            A70(3,N)= A70(0,N)+ A70(1,N)*Y02
            T02      = A71( 4,N)*Y02
            A71( 7,N)= A71( 1,N)*H21
            A71( 8,N)=(A71( 2,N)+ A71( 3,N)*Y02)*Y02
            A71( 9,N)= A71( 5,N)+T02
            A71(10,N)= T02
            A71(17,N)= A71(11,N)+(A71(12,N)+ A71(13,N)*Y02)*Y02
            A71(18,N)= A71(14,N)+ A71(15,N)*Y02
            T02      = A72( 4,N)*Y02
            A72( 7,N)= A72( 1,N)*H21
            A72( 8,N)=(A72( 2,N)+ A72( 3,N)*Y02)*Y02
            A72( 9,N)= A72( 5,N)+T02
            A72(10,N)= T02
            A72(17,N)= A72(11,N)+(A72(12,N)+ A72(13,N)*Y02)*Y02
            A72(18,N)= A72(14,N)+ A72(15,N)*Y02
C V(0,J,K,0)
            A74( 6,N)= SE0KKL(N)*H21
            A74( 7,N)=(A74( 1,N)+ A74( 2,N)*Y02)*Y01
            A74( 8,N)=-A74( 3,N)*Y01
            A74(15,N)= A74( 9,N)+(A74(10,N)+ A74(11,N)*Y02)*Y02
            A74(16,N)= A74(12,N)+ A74(13,N)*Y02
            A75( 6,N)= SE0LKL(N)*H21
            A75( 7,N)=(A75( 1,N)+ A75( 2,N)*Y02)*Y01
            A75( 8,N)= A75( 4,N)+ A75( 3,N)*Y02
            A75(15,N)= A75( 9,N)+(A75(10,N)+ A75(11,N)*Y02)*Y02
            A75(16,N)= A75(12,N)+ A75(13,N)*Y02
            A76( 6,N)= SE0JLL(N)*H21
            A76( 7,N)= A76( 1,N)+ A76( 2,N)*Y02
            A76( 8,N)= A76( 4,N)+ A76( 5,N)*Y02
            A76(15,N)= A76( 9,N)+(A76(10,N)+ A76(11,N)*Y02)*Y02
            A76(16,N)= A76(12,N)+ A76(13,N)*Y02
  550    CONTINUE
         DO 560 M=0,3
            A60(3,M)= A60(1,M)*Y02
            A61(3,M)= A61(0,M)+ A61(1,M)*Y02
            A73( 7,M)= A73( 1,M)*H21
            A73( 8,M)= A73( 2,M)+ A73( 3,M)*Y02
            A73( 9,M)= A73( 5,M)+ A73( 6,M)*Y02
            A73(17,M)= A73(11,M)+(A73(12,M)+ A73(13,M)*Y02)*Y02
            A73(18,M)= A73(14,M)+ A73(15,M)*Y02
            A77(21,1,M)=(A77( 0,1,M)+ A77( 1,1,M)*Y02)*Y1Y
     2                 +(A77( 9,1,M)+ A77(13,1,M)
     3                 +(A77(10,1,M)+ A77(12,1,M)- A77(13,1,M))*Y02)*H21
            T01= A77(11,1,M)*H21
            A77(22,1,M)= A77( 2,1,M)*Y1Y+T01
            A77(23,1,M)= A77( 6,1,M)
     2                 +(A77( 3,1,M)- A77( 6,1,M)+ A77( 7,1,M)
     3                 +(A77( 4,1,M)- A77( 7,1,M))*Y02)*Y02+T01
            A77(24,1,M)= A77( 8,1,M)+(A77( 5,1,M)- A77( 8,1,M))*Y02
            A77(21,2,M)=(A77( 0,2,M)+ A77( 1,2,M)*Y02)*H21
     2                 -(A77( 9,2,M)
     3                 +(A77( 3,2,M)- A77( 9,2,M)+ A77(10,2,M)
     4                 +(A77( 4,2,M)- A77(10,2,M)+ A77(11,2,M)
     5                 +(A77( 5,2,M)- A77(11,2,M))*Y02)*Y02)*Y02)
            A77(22,2,M)= A77( 2,2,M)*H21
     2                 + A77(12,2,M)
     3                 +(A77( 6,2,M)- A77(12,2,M)+ A77(13,2,M)
     4                 +(A77( 7,2,M)- A77(13,2,M))*Y02)*Y02
            A77(23,2,M)= A77(14,2,M)+(A77( 8,2,M)- A77(14,2,M))*Y02
            A77(24,2,M)= A77(15,2,M)+(A77(16,2,M)+ A77(17,2,M)*Y02)*Y02
            A77(25,2,M)= A77(18,2,M)+ A77(19,2,M)*Y02
            A77(21,3,M)= A77( 0,3,M)+(A77( 1,3,M)+(A77( 2,3,M)
     2                              + A77( 3,3,M)*Y02)*Y02)*Y02
            A77(22,3,M)= A77( 4,3,M)+(A77( 5,3,M)+ A77( 6,3,M)*Y02)*Y02
            A77(23,3,M)= A77( 7,3,M)+ A77( 8,3,M)*Y02
            A77(24,3,M)= A77( 9,3,M)
            T01=(A78( 3,1,M)+ A78( 4,1,M)*Y02)*Y02
            T02= A78( 5,1,M)*Y02
            A78(21,1,M)=(A78( 0,1,M)+ A78( 1,1,M)*Y02)*Y02-T01+T02
     2                 -(A78(11,1,M)- A78(12,1,M)+ A78(13,1,M))*H21
            A78(22,1,M)= A78( 2,1,M)*Y02-T02
            A78(23,1,M)= T01
     2                 -(A78( 9,1,M)- A78(10,1,M)*Y02+ A78(12,1,M))*H21
            A78(24,1,M)= A78( 5,1,M)*Y02
            A78(25,1,M)= A78( 6,1,M)+ A78( 7,1,M)*Y02
            A78(26,1,M)= A78( 8,1,M)
            A78(27,1,M)= A78(10,1,M)*Y02*H21
            T01=(A78( 9,2,M)+(A78(10,2,M)+ A78(11,2,M)*Y02)*Y02)*Y02
            T02=(A78(12,2,M)+ A78(13,2,M)*Y02)*Y02
            A78(21,2,M)=(A78( 0,2,M)+ A78( 1,2,M)*Y02)*H21-T01
            A78(22,2,M)= A78( 2,2,M)*H21+T02
            A78(23,2,M)= A78( 3,2,M)+(A78( 4,2,M)+ A78( 5,2,M)*Y02)*Y02
            A78(24,2,M)= A78( 6,2,M)+ A78( 7,2,M)*Y02
            A78(25,2,M)= A78( 8,2,M)
            A78(26,2,M)= T01
            A78(27,2,M)= T02
            A78(28,2,M)= A78(14,2,M)*Y02
            A78(29,2,M)= A78( 3,2,M)+ A78(15,2,M)
     2                 +(A78( 4,2,M)+ A78(16,2,M)
     3                 +(A78( 5,2,M)+ A78(17,2,M))*Y02)*Y02
            A78(30,2,M)= A78( 6,2,M)+ A78(18,2,M)
     2                 +(A78( 7,2,M)+ A78(19,2,M))*Y02
            A78(31,2,M)= A78( 8,2,M)+ A78(20,2,M)
            A78(21,3,M)= A78( 0,3,M)+(A78( 1,3,M)+(A78( 2,3,M)
     2                              + A78( 3,3,M)*Y02)*Y02)*Y02
            A78(22,3,M)= A78( 4,3,M)+(A78( 5,3,M)+ A78( 6,3,M)*Y02)*Y02
            A78(23,3,M)= A78( 7,3,M)+ A78( 8,3,M)*Y02
            A78(24,3,M)= A78( 9,3,M)
            T01= A79( 0,1,M)
     2         +(A79( 1,1,M)+ A79( 6,1,M)+ A79( 7,1,M)*Y02)*Y02
            T02= A79( 2,1,M)+ A79( 8,1,M)*Y02
            T03= A79( 5,1,M)*Y01
            T04= A79(12,1,M)- A79(13,1,M)*Y02+ A79(15,1,M)
            A79(21,1,M)= T01-T04*H21
            A79(22,1,M)= T02+T03
            A79(23,1,M)=(A79( 3,1,M)+ A79( 4,1,M)*Y02)*Y01+T01-T02
     2                 +(A79(14,1,M)- A79(15,1,M)+ A79(16,1,M))*H21
            A79(24,1,M)= T03
            A79(25,1,M)= A79( 9,1,M)+ A79(10,1,M)*Y02
            A79(26,1,M)= A79(11,1,M)
            A79(27,1,M)= T04
            A79(28,1,M)= A79(14,1,M)- A79(15,1,M)+ A79(16,1,M)
            A79(29,1,M)=(A79(13,1,M)*Y02- A79(13,1,M))*H21
            A79(30,1,M)= A79(17,1,M)*Y01
            T01=(A79( 9,2,M)+(A79(10,2,M)+ A79(11,2,M)*Y02)*Y02)*Y01
            T02=(A79(12,2,M)+ A79(13,2,M)*Y02)*Y01
            A79(21,2,M)=(A79( 0,2,M)+ A79( 1,2,M)*Y02)*H21+T01
            A79(22,2,M)= A79( 2,2,M)*H21-T02
            A79(23,2,M)= A79( 3,2,M)+(A79( 4,2,M)+ A79( 5,2,M)*Y02)*Y02
            A79(24,2,M)= A79( 6,2,M)+ A79( 7,2,M)*Y02
            A79(25,2,M)= A79( 8,2,M)
            A79(26,2,M)= T01
            A79(27,2,M)= T02
            A79(28,2,M)= A79(14,2,M)*Y01
            A79(29,2,M)= A79( 3,2,M)+ A79(15,2,M)
     2                 +(A79( 4,2,M)+ A79(16,2,M)
     3                 +(A79( 5,2,M)+ A79(17,2,M))*Y02)*Y02
            A79(30,2,M)= A79( 6,2,M)+ A79(18,2,M)
     2                 +(A79( 7,2,M)+ A79(19,2,M))*Y02
            A79(31,2,M)= A79( 8,2,M)+ A79(20,2,M)
            A79(21,3,M)= A79( 0,3,M)+(A79( 1,3,M)+(A79( 2,3,M)
     2                              + A79( 3,3,M)*Y02)*Y02)*Y02
            A79(22,3,M)= A79( 4,3,M)+(A79( 5,3,M)+ A79( 6,3,M)*Y02)*Y02
            A79(23,3,M)= A79( 7,3,M)+ A79( 8,3,M)*Y02
            A79(24,3,M)= A79( 9,3,M)
  560    CONTINUE
         A80(2,1)= A80(0,1)*Y1Y+A80(1,1)*H21
         A80(2,2)= A80(0,2)*Y1Y+A80(1,2)*H21
         DO 562 N=1,2
            A81(14,N)=(A81( 0,N)+A81( 1,N)*Y02)*Y02
            A81(15,N)= A81( 2,N)*Y02
            A81(16,N)=(A81( 3,N)+A81( 4,N)*Y02)*Y02
            A81(17,N)= A81( 5,N)*Y02
            A81(18,N)= A81( 6,N)+A81( 7,N)*Y02
            A81(19,N)= A81( 8,N)
            A81(20,N)= A81( 9,N)-A81(10,N)*Y02
            A81(21,N)= A81(10,N)*Y02
            A81(22,N)= A81(11,N)
            A81(23,N)= A81(12,N)
            A81(24,N)= A81(13,N)
  562    CONTINUE
         DO 564 M=14,24
  564    A81( M,1)= A81( M,1)+(A81( M,2)- A81( M,1))*Y02
         Y22= Y02*Y02
         DO 566 M=0,5
            T02   = A81(7+M,3)*H21
  566    A81(14+M,3)=A81( M,3)*Y22-T02
         A81(14,1)= A81(14,1)-(A81(17,3)+ A81(18,3)*Y02)
         A81(15,1)= A81(15,1)- A81(16,3)
         A81(16,1)= A81(16,1)+(A81(14,3)+ A81(15,3)*Y02)
         A81(17,1)= A81(17,1)+ A81(16,3)
         A81(18,1)= A81(18,1)- A81( 6,3)*H21
         T02 =(T02+ A81(20,1)+ A81(23,1))*H21
         A81(20,1)=(A81(21,1)+ A81(19,3))*H21
         A81(14,1)= A81(14,1)- A81(16,1)+ A81(17,1)
     2            -(A81(22,1)- A81(23,1)+ A81(24,1))*H21
         A81(15,1)= A81(15,1)- A81(17,1)
         A81(16,1)= A81(16,1)- T02
         DO 570 N=1,2
            A82(21,N)=(A82( 0,N)+ A82( 1,N)*Y02)
            A82(22,N)= A82( 2,N)
            A82(23,N)=(A82( 3,N)+(A82( 4,N)+ A82( 5,N)*Y02)*Y02)
            A82(24,N)=(A82( 6,N)+ A82( 7,N)*Y02)
            A82(25,N)= A82( 8,N)
            A82(26,N)=(A82( 9,N)+(A82(10,N)+ A82(11,N)*Y02)*Y02)*Y02
            A82(27,N)=(A82(12,N)+ A82(13,N)*Y02)*Y02
            A82(28,N)= A82(14,N)*Y02
            A82(29,N)= A82( 3,N)+ A82(15,N)+(A82( 4,N)+ A82(16,N)
     2               +(A82( 5,N)+ A82(17,N))*Y02)*Y02
            A82(30,N)= A82( 6,N)+ A82(18,N)+(A82( 7,N)+ A82(19,N))*Y02
            A82(31,N)= A82( 8,N)+ A82(20,N)
  570    CONTINUE
         DO 572 M=21,31
  572    A82( M,1)= A82( M,1)+(A82( M,2)-A82( M,1))*Y02
         DO 574 N=1,2
            A83(15,N)= A83( 0,N)+ A83( 1,N)*Y02
            A83(16,N)= A83( 2,N)
            A83(17,N)= A83( 3,N)+(A83( 4,N)+ A83( 5,N)*Y02)*Y02
            A83(18,N)= A83( 6,N)+ A83( 7,N)*Y02
            A83(19,N)= A83( 8,N)
            A83(20,N)=(A83( 9,N)+(A83(10,N)+ A83(11,N)*Y02)*Y02)*H21
            A83(21,N)=(A83(12,N)+ A83(13,N)*Y02)*H21
            A83(22,N)= A83(14,N)*H21
  574    CONTINUE
         T01= A83(15,1)*H21
         T02= A83(16,1)*H21
         A82(21,1)=(A82(21,1)- A86( 0)*H21)*H21
     2                           - A82(26,1)-(A83(17,1)*Y22+ A83(20,1))
         A82(22,1)= A82(22,1)*H21+ A82(27,1)+(A83(18,1)*Y22+ A83(21,1))
         A82(23,1)= A82(23,1)+ A83(20,2)+T01
         A82(24,1)= A82(24,1)+ A83(21,2)+T01-T02
         A82(25,1)= A82(25,1)- A83(22,2)+T02
         A82(28,1)= A82(28,1)-(A83(19,1)*Y22+ A83(22,1))
         A83(15,2)= A83(15,2)*H21- A82(29,1)
         A83(16,2)= A83(16,2)*H21+ A82(30,1)
         DO 576 N=1,2
            A84(10,N)= A84( 0,N)+(A84( 1,N)+(A84( 2,N)
     2                                     + A84( 3,N)*Y02)*Y02)*Y02
            A84(11,N)= A84( 4,N)+(A84( 5,N)+ A84( 6,N)*Y02)*Y02
            A84(12,N)= A84( 7,N)+ A84( 8,N)*Y02
            A84(13,N)= A84( 9,N)
  576    CONTINUE
         DO 578 M=10,13
  578    A84( M,1)= A84( M,1)+(A84( M,2)- A84( M,1))*Y02
         A84(10,1)= A84(10,1)-(A86(1)+(A86(2)+ A86(3)*Y02)*Y02)*H21
         A84(11,1)= A84(11,1)+(A86(4)+ A86(5)*Y02)*H21
         A84(12,1)= A84(12,1)- A86(6)*H21
         DO 580 N=1,2
            A85(10,N)= A85( 0,N)+(A85( 1,N)+(A85( 2,N)
     2                                     + A85( 3,N)*Y02)*Y02)*Y02
            A85(11,N)= A85( 4,N)+(A85( 5,N)+ A85( 6,N)*Y02)*Y02
            A85(12,N)= A85( 7,N)+ A85( 8,N)*Y02
            A85(13,N)= A85( 9,N)
  580    CONTINUE
         A85(10,2)= A85(10,2)- A85(10,1)+ A85(11,1)
         A85(11,2)= A85(11,2)- A85(11,1)+ A85(12,1)
         A85(12,2)= A85(12,2)- A85(12,1)+ A85(13,1)
         A85(13,2)= A85(13,2)- A85(13,1)
         Q4(0)= S4( 0)+(S4( 1)+(S4( 2)+(S4( 3)+S4( 4)*Y02)*Y02)*Y02)*Y02
         Q4(1)= S4( 5)+(S4( 6)+(S4( 7)+ S4( 8)*Y02)*Y02)*Y02
         Q4(2)= S4( 9)+(S4(10)+ S4(11)*Y02)*Y02
         Q4(3)= S4(12)+ S4(13)*Y02
         Q4(4)= S4(14)
      ENDIF
  590 CONTINUE
C
      IF(JTYPE.NE.4) THEN
         DO 610 M=0,3
            T(0)= A31(3,M)- A31(0,M)- A31(1,M)
            T(1)= A31(3,M)- A31(2,M)*Y04
            C31(1,M)= T(0)*F1(0)-T(1)*F1(1)
            C31(2,M)= T(0)*F2(0)-T(1)*F2(1)
            T(0)= A32(3,M)- A32(0,M)
            T(1)= A32(3,M)- A32(2,M)*Y04
            C32(1,M)= T(0)*F1(0)-T(1)*F1(1)
            C32(2,M)= T(0)*F2(0)-T(1)*F2(1)
  610    CONTINUE
         T(0)= A33( 6)
         T(1)= A33( 7)- A33( 8)*Y04
         T(2)= A33(15)-(A33(16)- A33(14)*Y04)*Y04
         FIJ00(1)= T(0)*F1(0)-T(1)*F1(1)+T(2)*F1(2)
         FIJ00(2)= T(0)*F2(0)-T(1)*F2(1)+T(2)*F2(2)
      ENDIF
      IF(JTYPE.LT.4) GO TO 690
C
      DO 620 M=0,3
         T(0)= A41(2,M)*Y04
         T(1)=(A41(3,M)-T(0))*XX1
         C41(1,M)= T(0)*F1(0)-T(1)*F1(1)
         C41(2,M)= T(0)*F2(0)-T(1)*F2(1)
         T(0)= A42(3,M)- A42(0,M)
         T(1)= A42(3,M)- A42(2,M)*Y04
         C42(1,M)= T(0)*F1(0)-T(1)*F1(1)
         C42(2,M)= T(0)*F2(0)-T(1)*F2(1)
  620 CONTINUE
      T01 = A43( 8)*Y04
      T(0)= T01
      T(1)= T01+(A43( 4)- A43( 5)*Y04)*Y04
     2         +(A43( 6)+ A43( 7)-T01)*XX1
      T(2)=(A43(15)-(A43(16)- A43(14)*Y04)*Y04)*XX1
      FI0K0(1)= T(0)*F1(0)-T(1)*F1(1)+T(2)*F1(2)
      FI0K0(2)= T(0)*F2(0)-T(1)*F2(1)+T(2)*F2(2)
      IF(JTYPE.LT.5) GO TO 690
C
      DO 630 N=1,3
         T(0)= A51(2,N)*Y04
         T(1)=(A51(3,N)-T(0))*XX1
         C51(1,N)= T(0)*F1(0)-T(1)*F1(1)
         C51(2,N)= T(0)*F2(0)-T(1)*F2(1)
         T(0)= A52(3,N)- A52(0,N)- A52(1,N)
         T(1)= A52(3,N)- A52(2,N)*Y04
         C52(1,N)= T(0)*F1(0)-T(1)*F1(1)
         C52(2,N)= T(0)*F2(0)-T(1)*F2(1)
C V(I,J,0,0)
         T(0)= A54( 6,N)
         T(1)= A54( 7,N)- A54( 8,N)*Y04
         T(2)= A54(15,N)-(A54(16,N)- A54(14,N)*Y04)*Y04
         C54(1,N)= T(0)*F1(0)-T(1)*F1(1)+T(2)*F1(2)
         C54(2,N)= T(0)*F2(0)-T(1)*F2(1)+T(2)*F2(2)
C V(I,0,K,0)
         T01 = A55( 8,N)*Y04
         T(0)= T01
         T(1)= T01+(A55( 4,N)- A55( 5,N)*Y04)*Y04
     2            +(A55( 6,N)+ A55( 7,N)-T01)*XX1
         T(2)=(A55(15,N)-(A55(16,N)- A55(14,N)*Y04)*Y04)*XX1
         C55(1,N)= T(0)*F1(0)-T(1)*F1(1)+T(2)*F1(2)
         C55(2,N)= T(0)*F2(0)-T(1)*F2(1)+T(2)*F2(2)
  630 CONTINUE
      DO 640 M=0,3
C V(0,J,K,0)
         T01 = A53( 8,M)*Y04
         T(0)= T01
         T(1)= T01+(A53( 3,M)+ A53( 4,M)- A53( 5,M)*Y04)*Y04
     2            +(A53( 6,M)- A53( 7,M)- T01)*XX1
         T(2)=(A53(15,M)-(A53(16,M)- A53(14,M)*Y04)*Y04)*XX1
         C53(1,M)= T(0)*F1(0)-T(1)*F1(1)+T(2)*F1(2)
         C53(2,M)= T(0)*F2(0)-T(1)*F2(1)+T(2)*F2(2)
  640 CONTINUE
      DO 650 M=0,M56
         T(0)=-A56(19,1,M)*Y04
         T(1)=(A56(18,1,M)-T(0))*XX1
     2       +(A56(20,1,M)+ A56(21,1,M)*Y04)*Y04
         T(2)=(A56(18,2,M)+(A56(19,2,M)- A56(20,2,M)*Y04)*Y04)*XX1
     2       -(A56(21,2,M)-(A56(22,2,M)+ A56(17,2,M)*Y04)*Y04)*Y04
         T(3)=(A56(18,3,M)+(A56(19,3,M)+(A56(20,3,M)
     2                                 + A56(21,3,M)*Y04)*Y04)*Y04)*XX1
         C56(1,M)= T(0)*F1(0)-T(1)*F1(1)-T(2)*F1(2)+T(3)*F1(3)
         C56(2,M)= T(0)*F2(0)-T(1)*F2(1)-T(2)*F2(2)+T(3)*F2(3)
  650 CONTINUE
      IF(JTYPE.LT.6) GO TO 690
C
      T02 = A62( 4)-A62( 8)
      T(0)= T02*Y03
      T(1)=(A62( 6)+ A62( 7)-T(0)+T02)*XX1
     2    -(A62( 8)+ A62( 5)*Y04)*Y03
      T(2)=(A62(15)-(A62(16)-A62(14)*Y04)*Y04)*XX1
      FI00L(1)= T(0)*F1(0)-T(1)*F1(1)+T(2)*F1(2)
      FI00L(2)= T(0)*F2(0)-T(1)*F2(1)+T(2)*F2(2)
      T02 = A63( 3)+A63( 4)-A63( 8)
      T(0)= T02*Y03
      T(1)=(A63( 6)- A63( 7)-T(0)+T02)*XX1
     2    -(A63( 8)+ A63( 5)*Y04)*Y03
      T(2)=(A63(15)-(A63(16)-A63(14)*Y04)*Y04)*XX1
      F0J0L(1)= T(0)*F1(0)-T(1)*F1(1)+T(2)*F1(2)
      F0J0L(2)= T(0)*F2(0)-T(1)*F2(1)+T(2)*F2(2)
      DO 660 N=1,3
         T(0)= A64(3,N)- A64(0,N)- A64(1,N)
         T(1)= A64(3,N)- A64(2,N)*Y04
         C64(1,N)= T(0)*F1(0)-T(1)*F1(1)
         C64(2,N)= T(0)*F2(0)-T(1)*F2(1)
         T(0)=-A65(2,N)*Y03
         T(1)=(A65(3,N)- A65(2,N)-T(0))*XX1
         C65(1,N)= T(0)*F1(0)-T(1)*F1(1)
         C65(2,N)= T(0)*F2(0)-T(1)*F2(1)
         T(0)= A66( 9,N)
         T(1)= A66( 7,N)+ A66( 8,N)*Y04
         T(2)= A66(17,N)-(A66(18,N)- A66(16,N)*Y04)*Y04
         C66(1,N)= T(0)*F1(0)-T(1)*F1(1)+T(2)*F1(2)
         C66(2,N)= T(0)*F2(0)-T(1)*F2(1)+T(2)*F2(2)
         T(0)=-A67(10,N)*Y03
         T(1)=(A67( 7,N)+ A67( 8,N)- A67(10,N)-T(0))*XX1
     2       -(A67( 9,N)+ A67( 6,N)*Y04)*Y03
         T(2)=(A67(17,N)-(A67(18,N)- A67(16,N)*Y04)*Y04)*XX1
         C67(1,N)= T(0)*F1(0)-T(1)*F1(1)+T(2)*F1(2)
         C67(2,N)= T(0)*F2(0)-T(1)*F2(1)+T(2)*F2(2)
         T(0)= A68(10,N)*Y03
         T(1)=(A68( 7,N)- A68( 8,N)+ A68(10,N)-T(0))*XX1
     2       -(A68( 9,N)+ A68( 6,N)*Y04)*Y03
         T(2)=(A68(17,N)-(A68(18,N)- A68(16,N)*Y04)*Y04)*XX1
         C68(1,N)= T(0)*F1(0)-T(1)*F1(1)+T(2)*F1(2)
         C68(2,N)= T(0)*F2(0)-T(1)*F2(1)+T(2)*F2(2)
C V(0,0,L,0)
         T(0)= A69(2,N)*Y04
         T(1)=(A69(3,N)-T(0))*XX1
         C69(1,N)= T(0)*F1(0)-T(1)*F1(1)
         C69(2,N)= T(0)*F2(0)-T(1)*F2(1)
C V(0,0,0,L)
         T02 = A70(2,N)*Y04
         T(0)=-A70(2,N)+T02
         T(1)=(A70(3,N)-T02)*XX1
         C70(1,N)= T(0)*F1(0)-T(1)*F1(1)
         C70(2,N)= T(0)*F2(0)-T(1)*F2(1)
         T(0)= A71(10,N)*Y04
         T(1)=(A71( 7,N)+ A71( 8,N)-T(0))*XX1
     2       +(A71( 9,N)+ A71( 6,N)*Y04)*Y04
         T(2)=(A71(17,N)-(A71(18,N)- A71(16,N)*Y04)*Y04)*XX1
         C71(1,N)= T(0)*F1(0)-T(1)*F1(1)+T(2)*F1(2)
         C71(2,N)= T(0)*F2(0)-T(1)*F2(1)+T(2)*F2(2)
         T(0)=-A72(10,N)*Y03
         T(1)=(A72( 7,N)+ A72( 8,N)- A72(10,N)-T(0))*XX1
     2       -(A72( 9,N)+ A72( 6,N)*Y04)*Y03
         T(2)=(A72(17,N)-(A72(18,N)- A72(16,N)*Y04)*Y04)*XX1
         C72(1,N)= T(0)*F1(0)-T(1)*F1(1)+T(2)*F1(2)
         C72(2,N)= T(0)*F2(0)-T(1)*F2(1)+T(2)*F2(2)
C V(0,J,K,0)
         T01 = A74( 8,N)*Y04
         T(0)= T01
         T(1)= T01+(A74( 3,N)+ A74( 4,N)- A74( 5,N)*Y04)*Y04
     2            +(A74( 6,N)- A74( 7,N)-T01)*XX1
         T(2)=(A74(15,N)-(A74(16,N)- A74(14,N)*Y04)*Y04)*XX1
         C74(1,N)= T(0)*F1(0)-T(1)*F1(1)+T(2)*F1(2)
         C74(2,N)= T(0)*F2(0)-T(1)*F2(1)+T(2)*F2(2)
         T02 = A75( 3,N)+ A75( 4,N)- A75( 8,N)
         T(0)= T02*Y03
         T(1)=(A75( 6,N)- A75( 7,N)-T(0)+T02)*XX1
     2       -(A75( 8,N)+ A75( 5,N)*Y04)*Y03
         T(2)=(A75(15,N)-(A75(16,N)- A75(14,N)*Y04)*Y04)*XX1
         C75(1,N)= T(0)*F1(0)-T(1)*F1(1)+T(2)*F1(2)
         C75(2,N)= T(0)*F2(0)-T(1)*F2(1)+T(2)*F2(2)
         T01 = A76( 6,N)*XX1
         T(0)= A76( 3,N)*Y3Y-T01
         T(1)=(T01- A76( 7,N)-(A76( 3,N)- A76( 7,N)- A76( 8,N)
     2                       -(A76( 3,N)+ A76( 3,N))*Y04)*Y04)*XX1
         T(2)=(A76(15,N)-(A76(16,N)- A76(14,N)*Y04)*Y04)*XX2
         C76(1,N)= T(0)*F1(0)-T(1)*F1(1)+T(2)*F1(2)
         C76(2,N)= T(0)*F2(0)-T(1)*F2(1)+T(2)*F2(2)
  660 CONTINUE
      DO 670 M=0,3
         T(0)= A60(3,M)
         T(1)= A60(3,M)+ A60(0,M)- A60(2,M)*Y04
         C60(1,M)= T(0)*F1(0)-T(1)*F1(1)
         C60(2,M)= T(0)*F2(0)-T(1)*F2(1)
         T01 =-A61(2,M)*Y03
         T(0)= T01
         T(1)=(A61(3,M)- A61(2,M)-T01)*XX1
         C61(1,M)= T(0)*F1(0)-T(1)*F1(1)
         C61(2,M)= T(0)*F2(0)-T(1)*F2(1)
         T(0)= A73( 4,M)*Y3Y- A73( 7,M)*XX1
         T(1)=(A73( 8,M)+T(0)
     2       -(A73( 8,M)+ A73( 9,M)+ A73( 4,M)*Y04)*Y04)*XX1
         T(2)=(A73(17,M)-(A73(18,M)- A73(16,M)*Y04)*Y04)*XX2
         C73(1,M)= T(0)*F1(0)+T(1)*F1(1)+T(2)*F1(2)
         C73(2,M)= T(0)*F2(0)+T(1)*F2(1)+T(2)*F2(2)
         T(0)= A77(22,1,M)*Y03
         T(1)=(A77(21,1,M)+ A77(22,1,M)*Y04)*XX1
     2       +(A77(23,1,M)+ A77(24,1,M)*Y04)*Y03
         T(2)=(A77(21,2,M)+(A77(22,2,M)+ A77(23,2,M)*Y04)*Y04)*XX1
     2       -(A77(24,2,M)-(A77(25,2,M)+ A77(20,2,M)*Y04)*Y04)*Y03
         T(3)=(A77(21,3,M)+(A77(22,3,M)+(A77(23,3,M)
     2                                 + A77(24,3,M)*Y04)*Y04)*Y04)*XX1
         C77(1,M)= T(0)*F1(0)-T(1)*F1(1)-T(2)*F1(2)+T(3)*F1(3)
         C77(2,M)= T(0)*F2(0)-T(1)*F2(1)-T(2)*F2(2)+T(3)*F2(3)
         T(0)=(A78(22,1,M)+ A78(24,1,M))*Y3Y+ A78(27,1,M)*XX1
         T(1)=(A78(23,1,M)+(A78(21,1,M)+A78(22,1,M)*Y04)*Y04)*XX1
     2       +(A78(25,1,M)+ A78(26,1,M)*Y04)*Y3Y- A78(27,1,M)*XX2
         T(2)=(A78(21,2,M)+(A78(22,2,M)+ A78(28,2,M)*Y04)*Y04)*XX2
     2       +(A78(23,2,M)-(A78(24,2,M)+ A78(25,2,M)*Y04)*Y04)*XX1
     3       -(A78(29,2,M)-(A78(30,2,M)+ A78(31,2,M)*Y04)*Y04)*Y41
         T(3)=(A78(21,3,M)+(A78(22,3,M)+(A78(23,3,M)
     2                                 + A78(24,3,M)*Y04)*Y04)*Y04)*XX2
         C78(1,M)= T(0)*F1(0)-T(1)*F1(1)-T(2)*F1(2)+T(3)*F1(3)
         C78(2,M)= T(0)*F2(0)-T(1)*F2(1)-T(2)*F2(2)+T(3)*F2(3)
         T(0)=-A79(24,1,M)*Y3Y+ A79(30,1,M)*H43
         T(1)=(A79(21,1,M)-(A79(23,1,M)+ A79(22,1,M)*Y04)*Y04)*XX1
     2       +(A79(25,1,M)+ A79(26,1,M)*Y04)*Y3Y- A79(29,1,M)*XX2
         T(2)=(A79(21,2,M)+(A79(22,2,M)- A79(28,2,M)*Y04)*Y04)*XX2
     2       +(A79(23,2,M)-(A79(24,2,M)+ A79(25,2,M)*Y04)*Y04)*XX1
     3       -(A79(29,2,M)-(A79(30,2,M)+ A79(31,2,M)*Y04)*Y04)*Y41
         T(3)=(A79(21,3,M)+(A79(22,3,M)+(A79(23,3,M)
     2                                 + A79(24,3,M)*Y04)*Y04)*Y04)*XX2
         C79(1,M)= T(0)*F1(0)-T(1)*F1(1)-T(2)*F1(2)+T(3)*F1(3)
         C79(2,M)= T(0)*F2(0)-T(1)*F2(1)-T(2)*F2(2)+T(3)*F2(3)
  670 CONTINUE
      T(0)= A80( 2,1)*Y3Y+ A80( 2,2)*H43
      T(1)=(A81(16,1)+(A81(14,1)+A81(15,1)*Y04)*Y04)*XX1
     2    +(A81(18,1)+ A81(19,1)*Y04)*Y3Y- A81(20,1)*XX2
      T(2)=(A82(21,1)+(A82(22,1)+ A82(28,1)*Y04)*Y04)*XX2
     2    +(A82(23,1)-(A82(24,1)+ A82(25,1)*Y04)*Y04)*XX1
     3    +(A83(15,2)+(A83(16,2)+ A82(31,1)*Y04)*Y04)*Y41
     4    -(A83(17,2)-(A83(18,2)- A83(19,2)*Y04)*Y04)*Y3Y
      T(3)=(A84(10,1)+(A84(11,1)+(A84(12,1)+A84(13,1)*Y04)*Y04)*Y04)*XX2
     2    +(A85(10,1)+(A85(10,2)+(A85(11,2)+(A85(12,2)
     3                          + A85(13,2)*Y04)*Y04)*Y04)*Y04)*XX1
      T(4)=(Q4(0)+(Q4(1)+(Q4(2)+(Q4(3)+Q4(4)*Y04)*Y04)*Y04)*Y04)*XX2
      FIJKL(1)= T(0)*F1(0)-T(1)*F1(1)-T(2)*F1(2)+T(3)*F1(3)+T(4)*F1(4)
      FIJKL(2)= T(0)*F2(0)-T(1)*F2(1)-T(2)*F2(2)+T(3)*F2(3)+T(4)*F2(4)
  690 CONTINUE
C
C BEGINING OF THE TWO PASSES THROUGH THE INTEGRAL CALCULATIONS
C
                  DO 800 ICP=1,2
                     T20= FQ(ICP-1)*WS(2)
                     VE0= FQ(ICP-1)*ES1
                     IF(JTYPE.EQ.3) THEN
C
C (SS,SS), (PS,SS), (PP,SS) AND (SP,SS) SECTION
C
                        T30= FQ(ICP-1)*WS(3)
                        DO 693 N=1,3
                           VEA(  N)= E(N,0,0,0)*T20  +C31(ICP,N)*WS(4)
                           VEA(3+N)= E(0,N,0,0)*T30  +C32(ICP,N)*WS(4)
  693                   CONTINUE
                        VE0= VE0+C32(ICP,0)*WS(2)
     2                          +C31(ICP,0)*WS(3)+FIJ00(ICP)*WS(4)
                        GO TO 700
                     ENDIF
                     T50= FQ(ICP-1)*WS(5)
                     IF(JTYPE.EQ.4) THEN
C
C (SS,SS), (PS,SS), (PS,PS) AND (SS,PS) SECTION
C
                        DO 694 N=1,3
                           VEA(  N)= E(N,0,0,0)*T20  +C41(ICP,N)*WS(6)
                           VEA(6+N)= E(0,0,N,0)*T50  +C42(ICP,N)*WS(6)
  694                   CONTINUE
                        VE0= VE0+C42(ICP,0)*WS(2)
     2                          +C41(ICP,0)*WS(5)+FI0K0(ICP)*WS(6)
                        GO TO 700
                     ENDIF
C
C (SS,SS), (PS,SS), (PP,SS), (SP,SS), (PS,PS), (SS,PS), (PP,PS) AND
C (SP,PS) SECTION
C
                     T30= FQ(ICP-1)*WS(3)
                     DO 695 N=1,3
                        VEA(  N)= E(N,0,0,0)*T20  +C31(ICP,N)*WS(4)
     2                           +C41(ICP,N)*WS(6)+C53(ICP,N)*WS(8)
                        VEA(3+N)= E(0,N,0,0)*T30  +C32(ICP,N)*WS(4)
     2                           +C51(ICP,N)*WS(7)+C55(ICP,N)*WS(8)
                        VEA(6+N)= E(0,0,N,0)*T50  +C42(ICP,N)*WS(6)
     2                           +C52(ICP,N)*WS(7)+C54(ICP,N)*WS(8)
  695                CONTINUE
                     VE0= VE0+C32(ICP,0)*WS(2)
     2                       +C31(ICP,0)*WS(3)+FIJ00(ICP)*WS(4)
     3                       +C41(ICP,0)*WS(5)+FI0K0(ICP)*WS(6)
     4                       +C53(ICP,0)*WS(7)+C56(ICP,0)*WS(8)
                     IF(JTYPE.LT.6) GO TO 700
C
C (PP,PP), (PP,SP), (PS,PP), (SP,PP), (SS,PP), (SP,SP), (PS,SP) AND
C (SS,SP) SECTION
C
                     T10= FQ(ICP-1)*WP(1)
                     DO 696 N=1,3
                        VEA(  N)= VEA(  N)
     2                           +C61(ICP,N)*WP(2)+C68(ICP,N)*WP(4)
     3                           +C73(ICP,N)*WP(6)+C79(ICP,N)*WP(8)
                        VEA(3+N)= VEA(3+N)
     2                           +C65(ICP,N)*WP(3)+C67(ICP,N)*WP(4)
     3                           +C76(ICP,N)*WP(7)+C78(ICP,N)*WP(8)
                        VEA(6+N)= VEA(6+N)
     2                           +C70(ICP,N)*WP(5)+C72(ICP,N)*WP(6)
     3                           +C75(ICP,N)*WP(7)+C77(ICP,N)*WP(8)
                        VEA(9+N)= E(0,0,0,N)*T10  +C60(ICP,N)*WP(2)
     2                           +C64(ICP,N)*WP(3)+C66(ICP,N)*WP(4)
     3                           +C69(ICP,N)*WP(5)+C71(ICP,N)*WP(6)
     4                           +C74(ICP,N)*WP(7)+C56(ICP,N)*WP(8)
  696                CONTINUE
                     VE0= VE0+C61(ICP,0)*WP(1)+FI00L(ICP)*WP(2)
     2                       +F0J0L(ICP)*WP(3)+C77(ICP,0)*WP(4)
     3                       +C73(ICP,0)*WP(5)+C78(ICP,0)*WP(6)
     4                       +C79(ICP,0)*WP(7)+FIJKL(ICP)*WP(8)
  700                CONTINUE
                     IF(ICP.EQ.1) THEN
                        VE0SUM= VE0SUM+VE0
                        FE0= F34*VE0
                        DO 720 M=1,3
                           VES(  M)= VES(  M)+VEA(  M)
                           VES(3+M)= VES(3+M)+VEA(3+M)
                           T01 = VEA(6+M)
                           T02 = VEA(9+M)
                           WES(  M)=-T01+(T01+T02)*Y03
  720                   CONTINUE
                     ELSE
                        FE1= X24*VE0
                        DO 730 M=1,3
                           T01 = VEA(  M)+VEA(3+M)
                           T02 = VEA(6+M)+VEA(9+M)
                           WES(3+M)=-T01+(T01+T02)*Z02
  730                   CONTINUE
                     ENDIF
  800             CONTINUE
C
C SUMMATION OF CONTRIBUTIONS FROM THE UNCONTRACTED GAUSSIANS
C
                  DO 905 M=1,3
                     WES(  M)= WES(  M)+P34(M,3)*FE0
                     WES(3+M)= WES(3+M)-P1(M)*FE1
                     VES(6+M)= VES(6+M)+WES(3+M)
                     FCS(3,M)= FCS(3,M)-WES(3+M)*Y03+WES(  M)
  905             CONTINUE
C
  910          KL=KL+1
  920       K0=K0+MXGSH
C
            TMP= F12*VE0SUM
            DO 925 M=1,3
               T01= VES(  M)
               T02= VES(3+M)+VES(6+M)
               T01=-T01+(T01+T02)*Y01+P12(M,3)*TMP
               FCS(1,M)= FCS(1,M)+T01
               FCS(2,M)= FCS(2,M)-T01+VES(6+M)
  925       CONTINUE
C
  930    IJ=IJ+1
  940 I0=I0+MXGSH
C
      RETURN
      END
