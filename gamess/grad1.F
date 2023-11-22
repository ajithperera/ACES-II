C*MODULE GRD1    *DECK STVDER
      SUBROUTINE xSTVDER
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL SOME,OUT,GOPARR,DSKWRK,MASWRK,EFLDL,SVDSKW,ISEPS,USEPS,
     *        UNCON,NXT
C
#include "gamess.par"
      PARAMETER (MXPT=100, MXFRG=50, MXFGPT=MXPT*MXFRG)
      PARAMETER (MXPSH=5*MXPT, MXPG=5*MXPSH)
      PARAMETER (MAXDEN=25*MXATM, LENABC=2000)
C
cwc      COMMON /xEFLDC / EVEC(3),EFLDL
cwc      COMMON /xEFPBUF/ POLCHG(10),NBUFMO,LBUFF(MXAO),LBF,NAPOL,IAPOL(10)
cwc      COMMON /xEFPOTD/ METHOF,ISCCHG,MOVE
cwc      COMMON /xFGRAD / DEF(3,MXFGPT),DEFT(3,MXFRG),TORQ(3,MXFRG),
cwc     *                EFCENT(3,MXFRG),FRGMAS(MXFRG),FRGMI(6,MXFRG),
cwc     *                ATORQ(3,MXFRG)
      COMMON /xFMCOM / X(1)
cwc      COMMON /xFRGINF/ NMPTS(MXFRG),NMTTPT,IEFC,IEFD,IEFQ,IEFO,
cwc     *                NPPTS(MXFRG),NPTTPT,IEFP,
cwc     *                NRPTS(MXFRG),NRTTPT,IREP,ICHGP,NFRG
      COMMON /xGRAD  / DE(3,MXATM)
cwc1
c buffer for mpi_allreduce
      dimension buf(3*MXATM)
cwc0
      COMMON /xINFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM)
      COMMON /xINTOPT/ ISCHWZ,IECP,NECP,IEFLD
cwc      COMMON /xISEPS / ISEPS, USEPS
      COMMON /xIOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /xOUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /xPAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /xPAULIN/ EX(MXPG,MXFRG),CS(MXPG,MXFRG),CP(MXPG,MXFRG),
     *                CD(MXPG,MXFRG),CF(MXPG,MXFRG),CG(MXPG,MXFRG),
     *                PRNAME(MXFGPT),PRCORD(3,MXFGPT),EFZNUC(MXFGPT),
     *                KSTART(MXPSH,MXFRG),KATOM(MXPSH,MXFRG),
     *                KTYPE(MXPSH,MXFRG),KNG(MXPSH,MXFRG),
     *                KLOC(MXPSH,MXFRG), KMIN(MXPSH,MXFRG),
     *                KMAX(MXPSH,MXFRG),NSHELL(MXFRG),NGAUSS(MXFRG),
     *                NATEF(MXFRG),NUMEF(MXFRG),NTPATM
C
cwc      COMMON /xPAULMO/ CENTNM(MXFGPT),CENTCD(3,MXFGPT),NORB(MXFRG),
cwc     *                NPBF(MXFRG),NTMO
      COMMON /xRESTAR/ TIMLIM,IREST,NREC,INTLOC,IST,JST,KST,LST
      COMMON /xRUNOPT/ RUNTYP,EXETYP,NEVALS
      COMMON /xSOLV  / FEPSI,RDS,DISEX2,COSURF(3,LENABC),SRAD(MXATM),
     *                QDEN(MAXDEN),AR(LENABC),
     *                NSPA,NPS,NPS2,NDEN,NPSPHER
      COMMON /xSIMDAT/ NACC,NREJ,IGOMIN,NRPA,IBWM,NACCT,NREJT,NRPAT,
     *                NPRTGO,IDPUNC,IGOFLG
cwc      COMMON /xZRFPAR/ GZRF,FIND(3),GNUCF,EBORN,DIELEC,IZRF,ICALC
cwc      COMMON /xRELWFN/ RMETHOD,QRQMT,CLIG,CLIG2,QRTOL,IQRORD,MODQR,NESOC,
cwc     *                NRATOM,NUMU,NQMTR,NQRDAF,MORDA,NDARELB
      COMMON /xMACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
C
      DATA CHECK/8HCHECK   /, ZRF/8HZRF     /,GRD1/8HGRD1    /
      DATA DEBUG/8HDEBUG   /,RESC/8HRESC    /,ANONE/8HNONE    /,
     *     ANESC/8HNESC    /
C
C     ----- DRIVER FOR THE GRADIENT OF THE ONE ELECTRON TERMS -----
C
cwc      IPOT=IEFC+IEFD+IEFQ+IEFO+IEFP+IREP+NTMO
C
      SOME = NPRINT.NE.-5 .AND. MASWRK
      OUT = (EXETYP.EQ.DEBUG  .OR.  EXETYP.EQ.GRD1   .OR.  NPRINT.EQ.-3)
     *         .AND.  MASWRK
cwc      UNCON=RMETHOD.EQ.RESC.AND.MOD(MODQR,2).EQ.1
cwc      IF(SOME) THEN
cwc         T0 = 0.0D+00
cwc         IF(IECP.GT.0  .OR.  IZRF.GT.0  .OR.  IPOT.GT.0 .OR. EFLDL)
cwc     *     CALL TSECND(T0)
cwc           WRITE(IW,9000)
cwc      END IF
cwc1
c      write (96,9000)
cwc0
C
C     INITIALIZE GRADIENT TO ZERO
C
cwc      NFRPTS = NMTTPT + NPTTPT + NRTTPT + NTPATM + NTMO
cwc1
      nmttpt=0
      nfrpts=0
cwc0
      CALL xVCLR(DE,1,3*NAT)
cwc      CALL xVCLR(DEF,1,3*MXFGPT)
cwc      IF(NFRG.GT.0) CALL VCLR(ATORQ,1,3*NFRG)
cwc      IF(NTMO.GT.0) CALL DPAULA
C
C     ----- GET WORKING STORAGE -----
C     NOTE THAT THESE THREE SECTIONS OVERLAP, AND THAT THEREFORE
C     EIJDEN MUST BE CALLED BEFORE DENDD1.  AFTER THOSE TWO CALLS,
C     ONLY -EPS- AND TOTAL DENSITY -DA- FROM THE 1ST TWO SECTIONS
C     WILL NEEDED BY ANY OF THE FOLLOWING CALLS.
C
      L1 = NUM
      L2 = (NUM*NUM+NUM)/2
      L3 = NUM*NUM
C
      CALL yVALFM(LOADFM)
      LEPS   = LOADFM + 1
      LVEC   = LEPS   + L2
      LEIG   = LVEC   + L3
      LIA    = LEIG   + L3
      LWRK   = LIA    + L1
      LAST1  = LWRK   + L1
C
      LDA    = LVEC
      LDB    = LDA+L2
      LAST2  = LDB+L2
C
      NDIMCD = 1
cwc      IF(ISEPS) NDIMCD=NAT+NPS
      LDRG   = LDB
      LEF3   = LDRG   + NAT*NAT
      LCOSDG = LEF3   + 3*NFRPTS
      LWORK  = LCOSDG + NDIMCD*NDIMCD
      LAST3  = LWORK  + 10*NMTTPT
C
      LAST = MAX(LAST1,LAST2,LAST3)
      LDERIV=LAST
      LDERIV1=LAST
C
cwc      IF(RMETHOD.EQ.RESC) THEN
cwc         L2U=L2
cwc         IF(UNCON) L2U=(NUMU*NUMU+NUMU)/2
cwc         NQRDAF=30
cwc         NTOTREC=4*3*NAT
cwc         LDERIV1=LDERIV+L2U*3*NAT
cwc         MORDA=LDERIV1+L2U*3*NAT
cwc         LAST = MORDA+(NTOTREC-1)/NWDVAR+1
cwc         LWRKL2=LAST
cwc         LUU=LAST
cwc         LWRKU=LAST
cwc         IF(UNCON) THEN
cwc            LUU=LWRKL2+L2
cwc            LWRKU=LUU+NUMU*L1
cwc            LAST=LWRKU+NUMU
cwc         ENDIF
cwc      ENDIF
      NEED = LAST - LOADFM - 1
      CALL yGETFM(NEED)
      IF(EXETYP.EQ.CHECK) GO TO 100
C
C     ----- GET LAGRANGIAN -EPS- AND TOTAL DENSITY -DA- MATRICES -----
C
      CALL xEIJDEN(X(LEPS),X(LVEC),X(LEIG),
     *            X(LIA),X(LWRK),L1,L2,L3,NQMT,X(LEIG))
cwc      IF (NTMO.GT.0) CALL QMADD(X(LEPS),X(LDA),L2,78)
      IF(OUT) THEN
cwc         WRITE(IW,*) 'THE LAGRANGIAN MATRIX IS'
cwc         CALL PRTRI(X(LEPS),L1)
      END IF
cwc1
c         WRITE(96,*) 'THE LAGRANGIAN MATRIX IS'
c         CALL yPRTRI(X(LEPS),L1)
cwc0
C
C     ----- GET TOTAL DENSITY -DA- MATRIX -----
C
      CALL xDENDD1(X(LDA),X(LDB),L2)
      IF(OUT) THEN
cwc         WRITE(IW,*) 'THE DENSITY MATRIX IS'
cwc         CALL PRTRI(X(LDA),L1)
      END IF
cwc1
c         WRITE(96,*) 'THE DENSITY MATRIX IS'
c         CALL yPRTRI(X(LDA),L1)
cwc0
C
cwc      IF(RMETHOD.EQ.RESC) THEN
cwc         SVDSKW=DSKWRK
cwc         DSKWRK=.TRUE.
cwc         CALL RAOPEN(NQRDAF,X(MORDA),0,NTOTREC,L3,NPRINT)
cwc         DSKWRK=SVDSKW
cwc      ENDIF
C
cwc      IF(RMETHOD.EQ.ANESC) CALL FLIPBASIS(0)
cwc      IF(UNCON) THEN
cwc         CALL FLIPBASIS(17)
cwc         L1=NUMU
cwc         L2 = (L1*L1+L1)/2
cwc         L3=L1*L1
cwc      ENDIF
C
C     ----- NUCLEAR REPULSION FORCE -----
C     ----- DENSITY FORCE -----
C
      CALL xVNNDER(X(LDRG),NAT,X(LEF3),NFRPTS,X(LCOSDG),NDIMCD,X(LWORK))
cwc1
      uncon=.false.
cwc0
      CALL xSDER(X(LEPS),X(LDERIV),L1,L2,UNCON)
cwc      IF (NTMO.GT.0) CALL QMADD(X(LDA),X(LEPS),L2,79)
C
C       QMADD CAN FAIL TO WORK WITH UNCONTRACTED BASIS (?)
C       IN THIS CASE SDER COULD NOT CONTRACT LAGRANGIAN SINCE
C       IT IS COMPUTED IN THE UNCONTRACTED SPACE AND THE DENSITY
C       IN CONTRACTED.
C       READ IN THE TRANSFORMATION MATRIX (CONTRACTION COEFFICIENTS)
C
cwc      IF(UNCON) THEN
cwc        CALL DAREAD(IDAF,IODA,X(LUU),L1*NUM,NDARELB+19,0)
cwc        LI=LDERIV
cwc        NXT = IBTYP.EQ.1
cwc        NEXT = -1
cwc        MINE = -1
cwc        DO IATM=1,NAT
cwc          DO IXYZ=1,3
cwc            MINE = MINE + 1
cwc            IF(GOPARR.AND.NXT) THEN
cwc              IF(MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
cwc              IF(NEXT.NE.MINE) GO TO 50
cwc            END IF
cwc            IF(GOPARR.AND..NOT.NXT.AND.MOD(MINE,NPROC).NE.ME) GOTO 50
cwc            CALL GRADEN(DE(IXYZ,IATM),X(LI),X(LWRKL2),X(LUU),X(LEPS),
cwc     *                  X(LWRKU),UNCON,0)
cwc   50       CONTINUE
cwc            LI=LI+L2
cwc          ENDDO
cwc        ENDDO
cwc        IF(GOPARR.AND.NXT) CALL DDI_DLBRESET
cwc      ENDIF
C
C     ----- HELLMANN-FEYNMAN FORCE -----
C     ----- INTEGRAL FORCE (AO DERIVATIVE CONTRIBUTION) -----
C
      CALL xHELFEY(X(LDA),X(LDERIV),L2)
      CALL xTVDER(X(LDA),X(LDERIV),X(LDERIV1),L2,0)
C
C     DERIVATIVES OF PVP INTEGRALS (HELLMANN-FEYNMAN AND AO DERIVATIVE)
C     CALLS MUST BE IN THIS ORDER SINCE THE ARRAY IS ACCUMULATED
C
cwc      IF(RMETHOD.NE.ANONE) THEN
cwc         IF(RMETHOD.EQ.ANESC) CALL FLIPBASIS(1)
cwc         CALL HFPVP(X(LDA),X(LDERIV),L2)
cwc         CALL DPVPINT(X(LDA),X(LDERIV),L2)
cwcC
cwcC        NESC ADDS CORRECTIONS TO S AND T
cwcC
cwc         IF(RMETHOD.EQ.ANESC) THEN
cwcC                GET CORRECTION TO S (USE LAGRANGIAN!)
cwc            CALL TVDER(X(LEPS),X(LDERIV),X(LDERIV1),L2,1)
cwc            CALL FLIPBASIS(2)
cwcC                GET CORRECTION TO T (USE DENSITY!)
cwc            CALL TVDER(X(LDA),X(LDERIV),X(LDERIV1),L2,2)
cwc            CALL FLIPBASIS(3)
cwc         ENDIF
cwcC
cwc         IF(RMETHOD.EQ.RESC) CALL RESCGX(X(LDA),L1,L2,L3)
cwcC
cwcC        THERE IS NO GRADIENT DRIVER FOR NESC - IT IS NOT NEEDED!
cwcC
cwc         IF(UNCON) THEN
cwc            CALL FLIPBASIS(15)
cwc            L1=NUM
cwc            L2 = (L1*L1+L1)/2
cwc            L3=L1*L1
cwc         ENDIF
cwc         IF(SOME.AND.(IECP.GT.0.OR.IZRF.GT.0.OR.IPOT.GT.0.OR.EFLDL))THEN
cwc           CALL TSECND(T1)
cwc           TR = T1-T0
cwc           WRITE(IW,9015) TR
cwc           T0 = T1
cwc         END IF
cwc      ENDIF
      CALL xDENDD1(X(LDA),X(LDB),L2)
C
cwc      IF(SOME  .AND.
cwc     *  (IECP.GT.0 .OR. IZRF.GT.0 .OR. IPOT.GT.0 .OR. EFLDL)) THEN
cwc         CALL TSECND(T1)
cwc         TG = T1-T0
cwc         WRITE(IW,9020) TG
cwc         T0 = T1
cwc      END IF
C
C     ----- EFFECTIVE CORE POTENTIAL MODIFICATIONS -----
C
  100 CONTINUE
cwc      IF(IECP.GT.0) THEN
cwc        CALL ECPGRD(DE, X(LDA), EXETYP, SOME)
cwc        CALL TSECND(T0)
cwc      END IF
C
C     ----- GET SELF-CONSISTENT REACTION FIELD GRADIENT -----
C
cwc      IF(IZRF.EQ.1) THEN
cwc         DUMMY = 0.0D+00
cwc         CALL DIPDER(X(LDA),DUMMY,DUMMY,DUMMY,L2,ZRF)
cwc         IF(SOME) THEN
cwc            CALL TSECND(T1)
cwc            TG = T1-T0
cwc            WRITE(IW,9040) TG
cwc            T0 = T1
cwc         END IF
cwc      END IF
C
C     ----- EFFECTIVE FRAGMENT GRADIENT -----
C     EFFT CALCULATES TOTAL FORCE AND TORQUE ON EACH FRAGMENT
C
cwc      IF(IPOT.GT.0) THEN
cwcC
cwcC   IF NBUFMO.GT.0 THEN WE HAVE TO USE THE BUFFERLESS DENSITY,
cwcC   WHICH REQUIRES MORE MEMORY
cwcC
cwc         IF(NBUFMO.GT.0)THEN
cwc            LVF = LAST
cwc            LDBF  = LVF + L3
cwc            LOCCF   = LDBF  +  L2
cwc            LWRKF  = LOCCF  +  L1
cwc            LASTF   = LWRKF  +  L2
cwc            NEEDF = LASTF - LAST
cwc            CALL GETFM(NEEDF)
cwc            CALL XCOPY(L2,X(LDA),1,X(LWRKF),1)
cwc            CALL DENDBF(X(LDA),X(LDB),X(LDBF),X(LVF),X(LOCCF),NBUFMO,
cwc     *                L1,L2,L3)
cwc         END IF
cwc         CALL FRGDER(X(LDA),L2)
cwc         IF(NBUFMO.GT.0) THEN
cwc             CALL XCOPY(L2,X(LWRKF),1,X(LDA),1)
cwc             CALL RETFM(NEEDF)
cwc         END IF
cwc         IF(MOVE.NE.0) CALL EFFT
cwc         IF(SOME) THEN
cwc            CALL TSECND(T1)
cwc            TG = T1-T0
cwc            WRITE(IW,9050) TG
cwc            T0 = T1
cwc         END IF
cwc      END IF
C
C     ----- GET ELECTRIC FIELD CONTRIBUTION -----
C
cwc      IF (EFLDL) THEN
cwc         CALL EFLD2(X(LDA),L2)
cwc         IF(SOME) THEN
cwc            CALL TSECND(T1)
cwc            TG = T1-T0
cwc            WRITE(IW,9060) TG
cwc            T0 = T1
cwc         END IF
cwc      END IF
C
cwc      IF(RMETHOD.EQ.RESC) CALL RACLOS(NQRDAF,'KEEP')
      CALL yRETFM(NEED)
C
C     ----- GATHER ALL GRADIENT CONTRIBUTIONS -----
C     ----- SAVE COMPLETED 1E-GRADIENT -----
C
cwc      IF(GOPARR) CALL DDI_GSUMF(1500,DE,3*NAT)
cwc1
      call MPIAR (DE,buf,3*NAT)
cwc0
      CALL yDAWRIT(IDAF,IODA,DE,3*NAT,3,0)
C
C     ----- CHECK TIME -----
C
      IREST=4
      IST=1
      JST=1
      KST=1
      LST=1
cwc      IF(MASWRK .AND. NPRTGO.NE.2) WRITE(IW,9010)
cwc      CALL TEXIT(1,IREST)
cwc1
c      write (96,9010)
cwc0
      RETURN
C
 9000 FORMAT(/1X,'BEGINNING ONE ELECTRON GRADIENT...')
 9010 FORMAT(1X,'..... END OF 1-ELECTRON GRADIENT ......')
 9015 FORMAT(1X,'TIME TO DO RELATIVISTIC GRADIENT INTEGRALS=',F10.2)
 9020 FORMAT(1X,'TIME TO DO ORDINARY GRADIENT INTEGRALS=',F10.2)
 9040 FORMAT(1X,'TIME TO DO     SCRF GRADIENT INTEGRALS=',F10.2)
 9050 FORMAT(1X,'TIME TO DO FRAGMENT GRADIENT INTEGRALS=',F10.2)
 9060 FORMAT(1X,'TIME TO DO ELECTRIC FIELD CONTRIBUTION=',F10.2)
      END
c--------------------------------------------------------------------------
C*MODULE GRD1    *DECK TVDER
      SUBROUTINE xTVDER(DENAB,DTINT,DVINT,L2,MODUS)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
cwc1
#undef _PARALLEL_ACES_COM_
#include "parallel_aces.com"
cwc0
      LOGICAL OUT,DBG,NORM
      LOGICAL GOPARR,DSKWRK,MASWRK,NXT,SVDSKW
C
#include "gamess.par"
      PARAMETER (MXCHRM=1)
C
      DIMENSION DENAB(L2),DTINT(L2,3,NAT),DVINT(L2,3,NAT)
      DIMENSION DIJ(225), IJX(35), IJY(35), IJZ(35),
     *          XS(6,7), YS(6,7), ZS(6,7), XT(6,5), YT(6,5), ZT(6,5),
     *          DXS(5,5), DYS(5,5), DZS(5,5),
     *          DXT(5,5), DYT(5,5), DZT(5,5),
     *          XV(6,5,5), YV(6,5,5), ZV(6,5,5),
     *          DXV(5,5,5), DYV(5,5,5), DZV(5,5,5)
C
cwc      COMMON /xCHMGMS/ XCHM(MXCHRM),YCHM(MXCHRM),ZCHM(MXCHRM),
cwc     *                DXELMM(MXCHRM),DYELMM(MXCHRM),DZELMM(MXCHRM),
cwc     *                QCHM(MXCHRM),NCHMAT,KCHRMM
      COMMON /xDSTV  / XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,
     *                XJ,YJ,ZJ,NI,NJ,CX,CY,CZ
      COMMON /xFMCOM / X(1)
      COMMON /xGRAD  / DE(3,MXATM)
      COMMON /xIJPAIR/ IA(MXAO)
      COMMON /xINFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM)
      COMMON /xIOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /xOUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /xPAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /xRUNOPT/ RUNTYP,EXETYP,NEVALS
      COMMON /xROOT  / XX,U(13),W(13),NROOTS
cwc      COMMON /xRELWFN/ RMETHOD,QRQMT,CLIG,CLIG2,QRTOL,IQRORD,MODQR,NESOC,
cwc     *                NRATOM,NUMU,NQMTR,NQRDAF,MORDA,NDARELB
C
C     COSMO INFORMATION
C
      PARAMETER (MAXDEN=25*MXATM, LENABC=2000, NPPA=1082)
      LOGICAL ISEPS, USEPS
cwc      COMMON /xISEPS / ISEPS, USEPS
      COMMON /xSOLV  / FEPSI,RDS,DISEX2,COSURF(3,LENABC),SRAD(MXATM),
     *                QDEN(MAXDEN),AR(LENABC),
     *                NSPA,NPS,NPS2,NDEN,NPSPHER
      COMMON /xSOLVI / IATSP(LENABC+1),N0(2),NP1,NP2,ISKIP
      COMMON /xCOSMO3/ COSZAN(NPPA),CORZAN(3,NPPA)
C
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00,TWO=2.0D+00)
      PARAMETER (RLN10=2.30258D+00)
      PARAMETER (PI212=1.1283791670955D+00)
      PARAMETER (SQRT3=1.73205080756888D+00)
      PARAMETER (SQRT5=2.23606797749979D+00,SQRT7=2.64575131106459D+00)
C
      DATA DEBUG/8HDEBUG   /, DBUGME/8HTVDER   /, GRD1/8HGRD1    /,
     *     RESC/8HRESC    /,ANESC/8HNESC    /
C
      DATA IJX / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     *           4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     *           5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     *           3, 1, 3, 2, 2/
      DATA IJY / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     *           1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     *           1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     *           1, 3, 2, 3, 2/
      DATA IJZ / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     *           1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     *           1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     *           3, 3, 2, 2, 3/
C
C     ----- BASIS FUNCTION DERIVATIVE CONTRIBUTIONS TO GRADIENT -----
C     INTEGRALS ARE OF TYPE <II'/H/JJ> = <II'/T+V/JJ>
C     RESC RUNS STORE INTEGRALS WITHOUT CONTRACTING THEM WITH DENSITY
C
C     MODUS - USED ONLY WITH NESC METHOD AND ONLY FOR KIN. ENERGY T
C           = 0 BULK T GRADIENT
C           = 1 CORRECTION TO S DERIVATIVE (PROPORTIONAL TO T)
C           = 2 CORRECTION TO T DERIVATIVE
C
      DBG = (EXETYP.EQ.DEBUG  .OR.  EXETYP.EQ.DBUGME)  .AND.  MASWRK
      OUT = (EXETYP.EQ.DEBUG  .OR.  EXETYP.EQ.DBUGME  .OR.
     *       EXETYP.EQ.GRD1   .OR.  NPRINT.EQ.-3)
cwc      IF(DBG) WRITE(IW,9000)
cwc1
c      write (96,9000)
cwc0
      IAZ=0
C
      TOL = RLN10*ITOL
      NORM = NORMF .NE. 1 .OR. NORMP .NE. 1
cwc      IF(RMETHOD.EQ.ANESC) THEN
cwc         FSC=ONE/CLIG
cwc         FSC2=FSC*FSC/TWO
cwc      ENDIF
cwc      IF(RMETHOD.EQ.RESC) THEN
cwc         CALL VCLR(DTINT,1,L2*3*NAT)
cwc         CALL VCLR(DVINT,1,L2*3*NAT)
cwc      ENDIF
C
C INITIALIZE PARALLEL
C
      NXT = IBTYP.EQ.1
      IPCOUNT = ME - 1
      NEXT = -1
      MINE = -1
cwc1
      nloop=-1
cwc0
C
C     -D*ELMM- CONTAINS CONTRIBUTION TO THE FORCES ON THE
C     EXTERNAL CHARGES WHICH ARE DUE TO QM MOLECULE.
C
cwc      IF(NCHMAT.NE.0) THEN
cwc         DO 100 ICHARM=1,NCHMAT
cwc            DXELMM(ICHARM)=ZERO
cwc            DYELMM(ICHARM)=ZERO
cwc            DZELMM(ICHARM)=ZERO
cwc  100    CONTINUE
cwc      END IF
C
C     FIRST, LOOP OVER THE GENUINE QM MOLECULE (ICHARM=0 PASS),
C     AND THEN ANY EXTERNAL CHARGES WHICH CHARMM MAY BE USING.
C     THE LATTER MAKE A CONTRIBUTION TO THE GRADIENT OF THE ATOMS
C     IN THE QM MOLECULE, AS WELL AS TO THEMSELVES.
C
cwc1
      nchmat=0
cwc0
      DO 2000 ICHARM=0,NCHMAT
C
C     ----- I SHELL
C
      DO 1400 II = 1,NSHELL
C
C           GO PARALLEL!
C
        IF(NXT .AND. GOPARR) THEN
           MINE = MINE + 1
cwc           IF(MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
cwc           IF(NEXT.NE.MINE) GO TO 1400
        END IF
C
        IAT = KATOM(II)
        XI = C(1,IAT)
        YI = C(2,IAT)
        ZI = C(3,IAT)
        I1 = KSTART(II)
        I2 = I1+KNG(II)-1
        LIT = KTYPE(II)
        MINI = KMIN(II)
        MAXI = KMAX(II)
        LOCI = KLOC(II)-MINI
        LITDER = LIT + 1
C
C     ----- J SHELL
C
        DO 1300 JJ = 1,NSHELL
C
C           GO PARALLEL!
C
cwc1
cwc          IF((.NOT.NXT) .AND. GOPARR) THEN
cwc             IPCOUNT = IPCOUNT + 1
cwc             IF(MOD(IPCOUNT,NPROC).NE.0) GO TO 1300
cwc          END IF
              nloop=nloop+1
              nloop=mod(nloop,nprocs)
              if (nloop.ne.irank) goto 1300
cwc0
C
          JAT = KATOM(JJ)
          XJ = C(1,JAT)
          YJ = C(2,JAT)
          ZJ = C(3,JAT)
          J1 = KSTART(JJ)
          J2 = J1+KNG(JJ)-1
          LJT = KTYPE(JJ)
          MINJ = KMIN(JJ)
          MAXJ = KMAX(JJ)
          LOCJ = KLOC(JJ)-MINJ
          LJTMOD = LJT+2
          NROOTS = (LIT+LJT-1)/2 + 1
          RR = (XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
C
C     ----- I PRIMITIVE
C
        DO 1100 IG = I1,I2
          AI = EX(IG)
          ARRI = AI*RR
          AXI = AI*XI
          AYI = AI*YI
          AZI = AI*ZI
          CSI=CS(IG)
          CPI=CP(IG)
          CDI=CD(IG)
          CFI=CF(IG)
          CGI=CG(IG)
C
C     ----- J PRIMITIVE
C
          DO 1000 JG = J1,J2
            AJ = EX(JG)
            AA = AI+AJ
            AA1 = ONE/AA
            DUM = AJ*ARRI*AA1
            IF(DUM .GT. TOL) GO TO 1000
            FAC = EXP(-DUM)
            CSJ = CS(JG)
            CPJ = CP(JG)
            CDJ = CD(JG)
            CFJ = CF(JG)
            CGJ = CG(JG)
            AX = (AXI+AJ*XJ)*AA1
            AY = (AYI+AJ*YJ)*AA1
            AZ = (AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR
C
            IJ = 0
            DUM1 = ZERO
            DUM2 = DUM1
            DO 600 I=MINI,MAXI
              IF(I.EQ.1) DUM1=CSI*FAC
              IF(I.EQ.2) DUM1=CPI*FAC
              IF(I.EQ.5) DUM1=CDI*FAC
              IF(I.EQ.8.AND.NORM) DUM1=DUM1*SQRT3
              IF(I.EQ.11) DUM1=CFI*FAC
              IF(I.EQ.14.AND.NORM) DUM1=DUM1*SQRT5
              IF(I.EQ.20.AND.NORM) DUM1=DUM1*SQRT3
              IF(I.EQ.21) DUM1=CGI*FAC
              IF(I.EQ.24.AND.NORM) DUM1=DUM1*SQRT7
              IF(I.EQ.30.AND.NORM) DUM1=DUM1*SQRT5/SQRT3
              IF(I.EQ.33.AND.NORM) DUM1=DUM1*SQRT3
C
              DO 600 J = MINJ,MAXJ
                IF(J.EQ.1) DUM2=DUM1*CSJ
                IF(J.EQ.2) DUM2=DUM1*CPJ
                IF(J.EQ.5) DUM2=DUM1*CDJ
                IF(J.EQ.8.AND.NORM) DUM2=DUM2*SQRT3
                IF(J.EQ.11) DUM2=DUM1*CFJ
                IF(J.EQ.14.AND.NORM) DUM2=DUM2*SQRT5
                IF(J.EQ.20.AND.NORM) DUM2=DUM2*SQRT3
                IF(J.EQ.21) DUM2=DUM1*CGJ
                IF(J.EQ.24.AND.NORM) DUM2=DUM2*SQRT7
                IF(J.EQ.30.AND.NORM) DUM2=DUM2*SQRT5/SQRT3
                IF(J.EQ.33.AND.NORM) DUM2=DUM2*SQRT3
C
                IJ=IJ+1
cwc                IF(RMETHOD.EQ.RESC) THEN
cwc                  DEN=ONE
cwc                  IF(LOCI+I.EQ.LOCJ+J) DEN=TWO
cwc                ELSE
                  NN = IA(MAX0(LOCI+I,LOCJ+J))+MIN0(LOCI+I,LOCJ+J)
                  DEN = DENAB(NN)
                  DEN = DEN + DEN
cwc                  IF(RMETHOD.EQ.ANESC) THEN
cwc                    IF(MODUS.EQ.1) DEN=DEN*FSC2
cwc                    IF(MODUS.EQ.2) DEN=-DEN
cwc                  ENDIF
cwc                ENDIF
                DIJ(IJ)=DUM2*DEN
  600         CONTINUE
C
C     -----  KINETIC ENERGY
C
      IF(ICHARM.GT.0) GO TO 700
C
      T = SQRT(AA1)
      X0 = AX
      Y0 = AY
      Z0 = AZ
      DO 660 J = 1,LJTMOD
        NJ =J
        DO 660 I = 1,LITDER
          NI = I
          CALL xVINT
          XS(I,J)=XINT*T
          YS(I,J)=YINT*T
          ZS(I,J)=ZINT*T
  660 CONTINUE
      CALL xDTXYZ(XT,YT,ZT,XS,YS,ZS,LITDER,LJT,AJ)
      CALL xDERI(DXS,DYS,DZS,XS,YS,ZS,LIT,LJT,AI)
      CALL xDERI(DXT,DYT,DZT,XT,YT,ZT,LIT,LJT,AI)
      IJ=0
      DO 680 I=MINI,MAXI
        IX=IJX(I)
        IY=IJY(I)
        IZ=IJZ(I)
        DO 670 J=MINJ,MAXJ
          JX=IJX(J)
          JY=IJY(J)
          JZ=IJZ(J)
          DUMX=DXT(IX,JX)* YS(IY,JY)* ZS(IZ,JZ)
     1        +DXS(IX,JX)* YT(IY,JY)* ZS(IZ,JZ)
     2        +DXS(IX,JX)* YS(IY,JY)* ZT(IZ,JZ)
          DUMY= XT(IX,JX)*DYS(IY,JY)* ZS(IZ,JZ)
     1        + XS(IX,JX)*DYT(IY,JY)* ZS(IZ,JZ)
     2        + XS(IX,JX)*DYS(IY,JY)* ZT(IZ,JZ)
          DUMZ= XT(IX,JX)* YS(IY,JY)*DZS(IZ,JZ)
     1        + XS(IX,JX)* YT(IY,JY)*DZS(IZ,JZ)
     2        + XS(IX,JX)* YS(IY,JY)*DZT(IZ,JZ)
          IJ=IJ+1
cwc          IF(RMETHOD.EQ.RESC) THEN
cwc             NN=IA(MAX0(LOCI+I,LOCJ+J))+MIN0(LOCI+I,LOCJ+J)
cwc             DUMINT=DIJ(IJ)
cwc             DTINT(NN,1,IAT)=DTINT(NN,1,IAT)+DUMX*DUMINT
cwc             DTINT(NN,2,IAT)=DTINT(NN,2,IAT)+DUMY*DUMINT
cwc             DTINT(NN,3,IAT)=DTINT(NN,3,IAT)+DUMZ*DUMINT
cwc          ELSE
             DE(1,IAT)=DE(1,IAT)+ DUMX*DIJ(IJ)
             DE(2,IAT)=DE(2,IAT)+ DUMY*DIJ(IJ)
             DE(3,IAT)=DE(3,IAT)+ DUMZ*DIJ(IJ)
cwc          ENDIF
  670     CONTINUE
  680   CONTINUE
cwc        IF(RMETHOD.EQ.ANESC.AND.MODUS.NE.0) GOTO 1000
C
C     ..... NUCLEAR ATTRACTION
C
  700 CONTINUE
      AAX = AA*AX
      AAY = AA*AY
      AAZ = AA*AZ
      IF(ICHARM.EQ.0) THEN
         MAXIC=NAT
cwc         IF(ISEPS) MAXIC=NAT+NPS
      ELSE
         MAXIC=1
      END IF
C
      DO 790 IC = 1,MAXIC
         IF(ICHARM.EQ.0) THEN
            IF(IC.LE.NAT) THEN
               ZNUC = -ZAN(IC)
               CX = C(1,IC)
               CY = C(2,IC)
               CZ = C(3,IC)
            ELSE
               IAZ=IATSP(IC-NAT)
               ZNUC = -COSZAN(IC-NAT)
               CX = CORZAN(1,IC-NAT)
               CY = CORZAN(2,IC-NAT)
               CZ = CORZAN(3,IC-NAT)
            END IF
         ELSE
cwc            ZNUC = -QCHM(ICHARM)
cwc            CX = XCHM(ICHARM)
cwc            CY = YCHM(ICHARM)
cwc            CZ = ZCHM(ICHARM)
         END IF
         XX = AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
         IF(NROOTS.LE.3) CALL xRT123
         IF(NROOTS.EQ.4) CALL xROOT4
         IF(NROOTS.EQ.5) CALL xROOT5
         DO 740 K = 1,NROOTS
            UU = AA*U(K)
            WW = W(K)*ZNUC
            TT = ONE/(AA+UU)
            T = SQRT(TT)
            X0 = (AAX+UU*CX)*TT
            Y0 = (AAY+UU*CY)*TT
            Z0 = (AAZ+UU*CZ)*TT
            DO 730 J = 1,LJT
               NJ = J
               DO 720 I = 1,LITDER
                  NI = I
                  CALL xVINT
                  XV(I,J,K) = XINT
                  YV(I,J,K) = YINT
                  ZV(I,J,K) = ZINT*WW
  720          CONTINUE
  730       CONTINUE
            CALL xDERI(DXV(1,1,K),DYV(1,1,K),DZV(1,1,K),
     *                 XV(1,1,K), YV(1,1,K), ZV(1,1,K),LIT,LJT,AI)
  740    CONTINUE
         IJ=0
         DO 780 I=MINI,MAXI
           IX=IJX(I)
           IY=IJY(I)
           IZ=IJZ(I)
           DO 780 J=MINJ,MAXJ
             JX=IJX(J)
             JY=IJY(J)
             JZ=IJZ(J)
             DUMX=ZERO
             DUMY=ZERO
             DUMZ=ZERO
             DO 770 K=1,NROOTS
               DUMX=DUMX+DXV(IX,JX,K)* YV(IY,JY,K)* ZV(IZ,JZ,K)
               DUMY=DUMY+ XV(IX,JX,K)*DYV(IY,JY,K)* ZV(IZ,JZ,K)
               DUMZ=DUMZ+ XV(IX,JX,K)* YV(IY,JY,K)*DZV(IZ,JZ,K)
  770        CONTINUE
             IJ=IJ+1
             IF((IC.GT.NAT).AND.(IAT.EQ.IAZ)) GOTO 780
             DUMINT=DIJ(IJ)*AA1*PI212
cwc             IF(RMETHOD.EQ.RESC) THEN
cwc                NN=IA(MAX0(LOCI+I,LOCJ+J))+MIN0(LOCI+I,LOCJ+J)
cwc                DVINT(NN,1,IAT)=DVINT(NN,1,IAT)+DUMX*DUMINT
cwc                DVINT(NN,2,IAT)=DVINT(NN,2,IAT)+DUMY*DUMINT
cwc                DVINT(NN,3,IAT)=DVINT(NN,3,IAT)+DUMZ*DUMINT
cwcC               IE. NO CHARMS NOW WITH RESC
cwc             ELSE
                DE(1,IAT)=DE(1,IAT)+DUMX*DUMINT
                DE(2,IAT)=DE(2,IAT)+DUMY*DUMINT
                DE(3,IAT)=DE(3,IAT)+DUMZ*DUMINT
cwc                IF(ISEPS.AND.(IC.GT.NAT)) THEN
cwc                   DE(1,IAZ)=DE(1,IAZ)-DUMX*(DIJ(IJ)*AA1*PI212)
cwc                   DE(2,IAZ)=DE(2,IAZ)-DUMY*(DIJ(IJ)*AA1*PI212)
cwc                   DE(3,IAZ)=DE(3,IAZ)-DUMZ*(DIJ(IJ)*AA1*PI212)
cwc                END IF
cwc                IF(ICHARM.GT.0) THEN
cwc                  DXELMM(ICHARM)=DXELMM(ICHARM)+DUMX*DUMINT
cwc                  DYELMM(ICHARM)=DYELMM(ICHARM)+DUMY*DUMINT
cwc                  DZELMM(ICHARM)=DZELMM(ICHARM)+DUMZ*DUMINT
cwc                END IF
cwc             ENDIF
  780    CONTINUE
  790 CONTINUE
C
 1000 CONTINUE
 1100 CONTINUE
C
C     ----- END OF PRIMITIVE LOOPS -----
C
      IF(DBG) THEN
cwc         WRITE(IW,9100) II,JJ
cwc         CALL EGOUT(DE,NAT)
cwc1
c         WRITE(96,9100) II,JJ
c         CALL yEGOUT(DE,NAT)
cwc0
      END IF
 1300 CONTINUE
 1400 CONTINUE
C
C     ----- END OF SHELL LOOPS -----
C
 2000 CONTINUE
      IF(OUT) THEN
cwc         IF(GOPARR) CALL DDI_GSUMF(1503,DE,3*NAT)
         IF(MASWRK) THEN
cwc            WRITE(IW,9000)
cwc            CALL EGOUT(DE,NAT)
         END IF
cwc         IF(GOPARR) CALL XSCAL(3*NAT,ONE/NPROC,DE,1)
      END IF
cwc1
c            WRITE(96,9000)
c            CALL yEGOUT(DE,NAT)
cwc0
C
cwc      IF(GOPARR .AND. NXT) CALL DDI_DLBRESET
C
cwc      IF(RMETHOD.EQ.RESC) THEN
cwc        IF(GOPARR) CALL DDI_GSUMF(1503,DTINT,L2*3*NAT)
cwc        IF(GOPARR) CALL DDI_GSUMF(1504,DVINT,L2*3*NAT)
cwc        SVDSKW=DSKWRK
cwc        DSKWRK=.TRUE.
cwc         LREC=3*NAT
cwcC        STARTING RECORD FOR THE T DERIVATIVE INTEGRALS
cwc        DO I=1,NAT
cwc          DO J=1,3
cwc            LREC=LREC+1
cwc            CALL RAWRIT(NQRDAF,X(MORDA),DTINT(1,J,I),L2,LREC,0)
cwcC
cwcC           ADD TWO CONTRIBUTIONS, HELLMANN-FEYMAN AND DV GOT HERE
cwcC
cwc            CALL RAREAD(NQRDAF,X(MORDA),DTINT(1,J,I),L2,LREC+6*NAT,0)
cwc            CALL XAXPY(L2,ONE,DTINT(1,J,I),1,DVINT(1,J,I),1)
cwc            CALL RAWRIT(NQRDAF,X(MORDA),DVINT(1,J,I),L2,LREC+6*NAT,0)
cwc          ENDDO
cwc        ENDDO
cwc        DSKWRK=SVDSKW
cwc      ENDIF
C
      RETURN
C
 9000 FORMAT(/10X,38(1H-)/10X,'GRADIENT INCLUDING AO DERIVATIVE TERMS'/
     *        10X,38(1H-))
 9100 FORMAT(1X,'TVDER: SHELLS II,JJ=',2I5)
      END
c--------------------------------------------------------------------------
C*MODULE GRD1    *DECK VINT
      SUBROUTINE xVINT
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION H(28),W(28),MIN(7),MAX(7)
C
      COMMON /xDSTV  / XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,
     *                XJ,YJ,ZJ,NI,NJ,CX,CY,CZ
      COMMON /xHERMIT/ H1,H2(2),H3(3),H4(4),H5(5),H6(6),H7(7)
      COMMON /xWERMIT/ W1,W2(2),W3(3),W4(4),W5(5),W6(6),W7(7)
C
      PARAMETER (ZERO=0.0D+00)
C
      EQUIVALENCE (H(1),H1),(W(1),W1)
C
      DATA MIN /1,2,4,7,11,16,22/
      DATA MAX /1,3,6,10,15,21,28/
C
C     ----- GAUSS-HERMITE QUADRATURE USING MINIMUM POINT FORMULA -----
C
      XINT = ZERO
      YINT = ZERO
      ZINT = ZERO
      NPTS = (NI+NJ-2)/2+1
      IMIN = MIN(NPTS)
      IMAX = MAX(NPTS)
      DO 340 I = IMIN,IMAX
         DUM = W(I)
         PX = DUM
         PY = DUM
         PZ = DUM
         DUM = H(I)*T
         PTX = DUM+X0
         PTY = DUM+Y0
         PTZ = DUM+Z0
         AX = PTX-XI
         AY = PTY-YI
         AZ = PTZ-ZI
         BX = PTX-XJ
         BY = PTY-YJ
         BZ = PTZ-ZJ
         IF(PX+AX .EQ. PX) AX = ZERO
         IF(PY+AY .EQ. PY) AY = ZERO
         IF(PZ+AZ .EQ. PZ) AZ = ZERO
         IF(PX+BX .EQ. PX) BX = ZERO
         IF(PY+BY .EQ. PY) BY = ZERO
         IF(PZ+BZ .EQ. PZ) BZ = ZERO
         GO TO (180,164, 162, 160,140,120,100),NI
C
  100    PX = PX*AX
         PY = PY*AY
         PZ = PZ*AZ
  120    PX = PX*AX
         PY = PY*AY
         PZ = PZ*AZ
  140    PX = PX*AX
         PY = PY*AY
         PZ = PZ*AZ
  160    PX = PX*AX
         PY = PY*AY
         PZ = PZ*AZ
  162    PX = PX*AX
         PY = PY*AY
         PZ = PZ*AZ
  164    PX = PX*AX
         PY = PY*AY
         PZ = PZ*AZ
  180    GO TO (320,310,300,280,260,240,220,200),NJ
C
  200    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  220    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  240    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  260    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  280    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  300    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  310    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  320    CONTINUE
         XINT = XINT+PX
         YINT = YINT+PY
         ZINT = ZINT+PZ
  340 CONTINUE
      RETURN
      END
c--------------------------------------------------------------------------
C*MODULE GRD1    *DECK DERI
      SUBROUTINE xDERI(DXDI,DYDI,DZDI,X,Y,Z,LIT,LJT,AI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  X  (6,1), Y  (6,1), Z  (6,1)
      DIMENSION DXDI(5,1),DYDI(5,1),DZDI(5,1)
C
C     ----- DXDI ... -----
C
      DO 110 J=1,LJT
        DXDI(1,J)=X(2,J)*(AI+AI)
        DYDI(1,J)=Y(2,J)*(AI+AI)
        DZDI(1,J)=Z(2,J)*(AI+AI)
  110 CONTINUE
C
      IF(LIT.EQ.1) RETURN
C
      DO 120 I=2,LIT
        DO 120 J=1,LJT
          DXDI(I,J)=X(I+1,J)*(AI+AI) - X(I-1,J)*(I-1)
          DYDI(I,J)=Y(I+1,J)*(AI+AI) - Y(I-1,J)*(I-1)
          DZDI(I,J)=Z(I+1,J)*(AI+AI) - Z(I-1,J)*(I-1)
  120 CONTINUE
C
      RETURN
      END
c--------------------------------------------------------------------------
C*MODULE GRD1    *DECK DTXYZ
      SUBROUTINE xDTXYZ(XT,YT,ZT,XS,YS,ZS,NI,NJ,AJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION XT(6,*),YT(6,*),ZT(6,*)
      DIMENSION XS(6,*),YS(6,*),ZS(6,*)
      PARAMETER (THREE=3.0D+00)
C
      DO 10 I=1,NI
        XT(I,1)=(XS(I,1  )       - XS(I,3  )*(AJ+AJ))*AJ
        YT(I,1)=(YS(I,1  )       - YS(I,3  )*(AJ+AJ))*AJ
        ZT(I,1)=(ZS(I,1  )       - ZS(I,3  )*(AJ+AJ))*AJ
   10 CONTINUE
C
      IF(NJ.EQ.1) RETURN
C
      DO 20 I=1,NI
        XT(I,2)=(XS(I,2  )*THREE - XS(I,4  )*(AJ+AJ))*AJ
        YT(I,2)=(YS(I,2  )*THREE - YS(I,4  )*(AJ+AJ))*AJ
        ZT(I,2)=(ZS(I,2  )*THREE - ZS(I,4  )*(AJ+AJ))*AJ
   20 CONTINUE
C
      IF(NJ.EQ.2) RETURN
C
      DO 30 J=3,NJ
        FACT1 = J+J-1
        FACT2 = (J-1)*(J-2)
        FACT2 = FACT2/2
        DO 30 I=1,NI
          XT(I,J)=(XS(I,J  )*FACT1 - XS(I,J+2)*(AJ+AJ))*AJ
     *                             - XS(I,J-2)*FACT2
          YT(I,J)=(YS(I,J  )*FACT1 - YS(I,J+2)*(AJ+AJ))*AJ
     *                             - YS(I,J-2)*FACT2
          ZT(I,J)=(ZS(I,J  )*FACT1 - ZS(I,J+2)*(AJ+AJ))*AJ
     *                             - ZS(I,J-2)*FACT2
   30 CONTINUE
C
      RETURN
      END
c--------------------------------------------------------------------------
C*MODULE GRD1    *DECK HELFEY
      SUBROUTINE xHELFEY(DENAB,HFINT,L2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
cwc1
#undef _PARALLEL_ACES_COM_
#include "parallel_aces.com"
cwc0
      LOGICAL IANDJ,OUT,DBG,NORM,DOUBLE,GOPARR,DSKWRK,MASWRK,NXT,SVDSKW
C
C  MAGIC NUMBERS: S(1)+P(3)+D(6)+F(10)+G(15)=35
C                 BASIS FUNCTION TYPES=5
C                 MAXIMUM BASIS TYPE SQUARED = (15)**2 = 225
C
#include "gamess.par"
C
      DIMENSION DENAB(L2),HFINT(L2,3,NAT)
      DIMENSION IJX(35),IJY(35),IJZ(35),DIJ(225),
     *          XIN(5,5,5,2),YIN(5,5,5,2),ZIN(5,5,5,2)
C
      COMMON /xDSTV  / XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,
     *                XJ,YJ,ZJ,NI,NJ,CX,CY,CZ
      COMMON /xGRAD  / DE(3,MXATM)
      COMMON /xIJPAIR/ IA(MXAO)
      COMMON /xINFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM)
      COMMON /xIOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /xOUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /xPAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /xROOT  / XX,U(13),W(13),NROOTS
      COMMON /xRUNOPT/ RUNTYP,EXETYP,NEVALS
cwc      COMMON /xRELWFN/ RMETHOD,QRQMT,CLIG,CLIG2,QRTOL,IQRORD,MODQR,NESOC,
cwc     *                NRATOM,NUMU,NQMTR,NQRDAF,MORDA,NDARELB
      COMMON /xFMCOM / X(1)
C
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00)
      PARAMETER (RLN10=2.30258D+00)
      PARAMETER (PI212=1.1283791670955D+00)
      PARAMETER (SQRT3=1.73205080756888D+00)
      PARAMETER (SQRT5=2.23606797749979D+00)
      PARAMETER (SQRT7=2.64575131106459D+00)
C
      DATA DEBUG/8HDEBUG   /, DBUGME/8HHELFEY  /, GRD1/8HGRD1    /,
     *     RESC/8HRESC    /
C
C  THE IJX, IJY, AND IJZ ARRAYS CONTAIN THE POWERS OF THE CARTESIAN
C  GAUSSIANS PLUS 1 IN EVERY PLACE.  THE REASON FOR THE PLUS ONE IS
C  THAT THEY ARE USED TO INDEX AN ARRAY RUNNING FROM 1-5 CORRESPONDING
C  TO 0-4 IN THE POWER.
C
      DATA IJX/ 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     *          4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     *          5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     *          3, 1, 3, 2, 2/
      DATA IJY/ 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     *          1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     *          1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     *          1, 3, 2, 3, 2/
      DATA IJZ/ 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     *          1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     *          1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     *          3, 3, 2, 2, 3/
C
C     ----- HELMANN-FEYNMAN GRADIENT TERM -----
C     INTEGRAL TYPE IS <II/H'/JJ> = <II/V'/JJ>
C     RESC RUNS DO NOT CONTRACT THE INTEGRALS WITH DENSITY
C
      DBG = (EXETYP.EQ.DEBUG  .OR.  EXETYP.EQ.DBUGME)  .AND.  MASWRK
      OUT = (EXETYP.EQ.DEBUG  .OR.  EXETYP.EQ.DBUGME  .OR.
     *       EXETYP.EQ.GRD1   .OR.  NPRINT.EQ.-3)
cwc      IF(DBG) WRITE(IW,9000)
cwc1
c      write (96,9000)
cwc0
C
      TOL = RLN10*ITOL
      NORM = NORMF.NE.1 .OR. NORMP.NE.1
cwc      IF(RMETHOD.EQ.RESC) CALL VCLR(HFINT,1,L2*3*NAT)
C
C     INITIALIZATION FOR PARALLEL
C
      NXT = IBTYP.EQ.1
      IPCOUNT = ME - 1
      NEXT = -1
      MINE = -1
cwc1
      nloop=-1
cwc0
C
C     ----- I SHELL
C
      DO 600 II = 1,NSHELL
C
C           GO PARALLEL!
C
      IF(NXT .AND. GOPARR) THEN
         MINE = MINE + 1
cwc         IF(MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
cwc         IF(NEXT.NE.MINE) GO TO 600
      END IF
C
      I = KATOM(II)
      XI = C(1,I)
      YI = C(2,I)
      ZI = C(3,I)
      I1 = KSTART(II)
      I2 = I1+KNG(II)-1
      LIT = KTYPE(II)
      MINI = KMIN(II)
      MAXI = KMAX(II)
      LOCI = KLOC(II)-MINI
C
C     ----- J SHELL
C
      DO 580 JJ = 1,II
C
C           GO PARALLEL!
C
cwc        IF((.NOT.NXT) .AND. GOPARR) THEN
cwc           IPCOUNT = IPCOUNT + 1
cwc           IF(MOD(IPCOUNT,NPROC).NE.0) GO TO 580
cwc        END IF
              nloop=nloop+1
              nloop=mod(nloop,nprocs)
              if (nloop.ne.irank) goto 580
cwc0
C
        J = KATOM(JJ)
        XJ = C(1,J)
        YJ = C(2,J)
        ZJ = C(3,J)
        J1 = KSTART(JJ)
        J2 = J1+KNG(JJ)-1
        LJT = KTYPE(JJ)
        MINJ = KMIN(JJ)
        MAXJ = KMAX(JJ)
        LOCJ = KLOC(JJ)-MINJ
        NROOTS = (LIT+LJT+1-2)/2 + 1
        RR = (XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
        IANDJ = II .EQ. JJ
C
C     ----- I PRIMITIVE
C
        DO 520 IG = I1,I2
        AI = EX(IG)
        ARRI = AI*RR
        AXI = AI*XI
        AYI = AI*YI
        AZI = AI*ZI
        CSI = CS(IG)
        CPI = CP(IG)
        CDI = CD(IG)
        CFI = CF(IG)
        CGI = CG(IG)
C
C     ----- J PRIMITIVE
C
        JGMAX = J2
        IF(IANDJ) JGMAX = IG
        DO 500 JG = J1,JGMAX
          AJ = EX(JG)
          AA = AI+AJ
          AA1 = ONE/AA
          DUM = AJ*ARRI*AA1
          IF(DUM .GT. TOL) GO TO 500
          FAC = EXP(-DUM)
          CSJ = CS(JG)
          CPJ = CP(JG)
          CDJ = CD(JG)
          CFJ = CF(JG)
          CGJ = CG(JG)
          AX = (AXI+AJ*XJ)*AA1
          AY = (AYI+AJ*YJ)*AA1
          AZ = (AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR
C
          DOUBLE=IANDJ.AND.IG.NE.JG
          JMAX = MAXJ
          NN = 0
          DUM1 = ZERO
          DUM2 = DUM1
          DO 360 I = MINI,MAXI
            IF(I.EQ.1) DUM1=CSI*FAC
            IF(I.EQ.2) DUM1=CPI*FAC
            IF(I.EQ.5) DUM1=CDI*FAC
            IF(I.EQ.8.AND.NORM) DUM1=DUM1*SQRT3
            IF(I.EQ.11) DUM1=CFI*FAC
            IF(I.EQ.14.AND.NORM) DUM1=DUM1*SQRT5
            IF(I.EQ.20.AND.NORM) DUM1=DUM1*SQRT3
            IF(I.EQ.21) DUM1=CGI*FAC
            IF(I.EQ.24.AND.NORM) DUM1=DUM1*SQRT7
            IF(I.EQ.30.AND.NORM) DUM1=DUM1*SQRT5/SQRT3
            IF(I.EQ.33.AND.NORM) DUM1=DUM1*SQRT3
C
            IF(IANDJ) JMAX = I
            DO 360 J = MINJ,JMAX
              IF(J.EQ.1) THEN
                DUM2=DUM1*CSJ
                IF( .NOT. DOUBLE) GO TO 350
                IF(I .GT. 1) THEN
                 DUM2 = DUM2+CSI*CPJ*FAC
                ELSE
                 DUM2 = DUM2+DUM2
                END IF
              ELSE IF(J.EQ.2) THEN
                DUM2=DUM1*CPJ
                IF(DOUBLE) DUM2 = DUM2+DUM2
              ELSE IF(J.EQ.5) THEN
                DUM2=DUM1*CDJ
                IF(DOUBLE) DUM2 = DUM2+DUM2
              ELSE IF(J.EQ.8.AND.NORM) THEN
                DUM2 = DUM2*SQRT3
              ELSE IF(J.EQ.11) THEN
                DUM2=DUM1*CFJ
                IF(DOUBLE) DUM2=DUM2+DUM2
              ELSE IF(J.EQ.14.AND.NORM) THEN
                DUM2=DUM2*SQRT5
              ELSE IF(J.EQ.20.AND.NORM) THEN
                DUM2=DUM2*SQRT3
              ELSE IF(J.EQ.21) THEN
                DUM2=DUM1*CGJ
                IF(DOUBLE) DUM2=DUM2+DUM2
              ELSE IF(J.EQ.24.AND.NORM) THEN
                DUM2=DUM2*SQRT7
              ELSE IF(J.EQ.30.AND.NORM) THEN
                DUM2=DUM2*SQRT5/SQRT3
              ELSE IF(J.EQ.33.AND.NORM) THEN
                DUM2=DUM2*SQRT3
              END IF
C
  350         NN = NN+1
cwc              IF(RMETHOD.EQ.RESC) THEN
cwc                DIJ(NN)=DUM2*PI212*AA1
cwc              ELSE
                NDUM = IA(LOCI+I)+(LOCJ+J)
                DEN = DENAB(NDUM)
                IF(.NOT.IANDJ.OR.I.NE.J) DEN=DEN+DEN
                DIJ(NN)=DUM2*DEN*PI212*AA1
cwc              ENDIF
  360     CONTINUE
C
C     ..... HELLMANN-FEYNMAN TERM .....
C
          AAX = AA*AX
          AAY = AA*AY
          AAZ = AA*AZ
          DO 480 IC = 1,NAT
            ZNUC = -ZAN(IC)
            CX = C(1,IC)
            CY = C(2,IC)
            CZ = C(3,IC)
            XX = AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
            IF(NROOTS.LE.3) CALL xRT123
            IF(NROOTS.EQ.4) CALL xROOT4
            IF(NROOTS.EQ.5) CALL xROOT5
            DO 420 K = 1,NROOTS
              UU = AA*U(K)
              WW = W(K)*ZNUC
              WW=WW*(UU+UU)
              TT = ONE/(AA+UU)
              T = SQRT(TT)
              X0 = (AAX+UU*CX)*TT
              Y0 = (AAY+UU*CY)*TT
              Z0 = (AAZ+UU*CZ)*TT
              DO 400 J = 1,LJT
                NJ = J
                DO 400 I = 1,LIT
                  NI = I
                  CALL xVINT
                  XIN(I,J,K,1) = XINT
                  YIN(I,J,K,1) = YINT
                  ZIN(I,J,K,1) = ZINT*WW
                  CALL xDVINT
                  XIN(I,J,K,2) = XINT
                  YIN(I,J,K,2) = YINT
                  ZIN(I,J,K,2) = ZINT*WW
  400         CONTINUE
  420       CONTINUE
           IJ=0
           DO 460 I=MINI,MAXI
            IX=IJX(I)
            IY=IJY(I)
            IZ=IJZ(I)
            JMAX=MAXJ
            IF(IANDJ) JMAX=I
            DO 460 J=MINJ,JMAX
             JX=IJX(J)
             JY=IJY(J)
             JZ=IJZ(J)
             DUMX = ZERO
             DUMY = ZERO
             DUMZ = ZERO
             DO 440 K = 1,NROOTS
              DUMX = DUMX+XIN(IX,JX,K,2)*YIN(IY,JY,K,1)*ZIN(IZ,JZ,K,1)
              DUMY = DUMY+XIN(IX,JX,K,1)*YIN(IY,JY,K,2)*ZIN(IZ,JZ,K,1)
  440         DUMZ = DUMZ+XIN(IX,JX,K,1)*YIN(IY,JY,K,1)*ZIN(IZ,JZ,K,2)
             IJ=IJ+1
             DUM=DIJ(IJ)
cwc             IF(RMETHOD.EQ.RESC) THEN
cwc               NN=IA(LOCI+I)+(LOCJ+J)
cwc               HFINT(NN,1,IC)=HFINT(NN,1,IC)+DUMX*DUM
cwc               HFINT(NN,2,IC)=HFINT(NN,2,IC)+DUMY*DUM
cwc               HFINT(NN,3,IC)=HFINT(NN,3,IC)+DUMZ*DUM
cwc             ELSE
               DE(1,IC) = DE(1,IC)+DUM*DUMX
               DE(2,IC) = DE(2,IC)+DUM*DUMY
               DE(3,IC) = DE(3,IC)+DUM*DUMZ
cwc             ENDIF
cwc             IF(DBG) WRITE(IW,9100) IC,II,JJ,IJ,DIJ(IJ),DUMX,DUMY,
cwc     +                         DUMZ,DE(1,IC),DE(2,IC),DE(3,IC)
cwc1
c             WRITE(96,9100) IC,II,JJ,IJ,DIJ(IJ),DUMX,DUMY,
c     +                         DUMZ,DE(1,IC),DE(2,IC),DE(3,IC)
cwc0
  460       CONTINUE
  480     CONTINUE
  500 CONTINUE
  520 CONTINUE
C
C     ----- END OF *PRIMITIVE* LOOPS -----
C
      IF(DBG) THEN
cwc         WRITE(IW,9010) II,JJ
cwc         CALL EGOUT(DE,NAT)
      END IF
cwc1
c         WRITE(96,9010) II,JJ
c         CALL yEGOUT(DE,NAT)
cwc0
  580 CONTINUE
  600 CONTINUE
C
C     ----- END OF *SHELL* LOOPS -----
C
      IF(OUT) THEN
cwc         IF(GOPARR) CALL DDI_GSUMF(1502,DE,3*NAT)
cwc         IF(MASWRK) THEN
cwc            WRITE(IW,9000)
cwc            CALL EGOUT(DE,NAT)
cwc         END IF
cwc         IF(GOPARR) CALL XSCAL(3*NAT,ONE/NPROC,DE,1)
cwc1
c            WRITE(96,9000)
c            CALL yEGOUT(DE,NAT)
cwc0
      END IF
cwc      IF(GOPARR .AND. NXT) CALL DDI_DLBRESET
cwc      IF(RMETHOD.EQ.RESC) THEN
cwc        IF(GOPARR) CALL DDI_GSUMF(1502,HFINT,L2*3*NAT)
cwc        SVDSKW=DSKWRK
cwc        DSKWRK=.TRUE.
cwc        LREC=3*3*NAT
cwcC       STARTING RECORD FOR THE HELLMANN-FEYNMAN INTEGRALS
cwc        DO I=1,NAT
cwc          DO J=1,3
cwc            LREC=LREC+1
cwc            CALL RAWRIT(NQRDAF,X(MORDA),HFINT(1,J,I),L2,LREC,0)
cwc          ENDDO
cwc        ENDDO
cwc        DSKWRK=SVDSKW
cwc      ENDIF
C
      RETURN
C
 9000 FORMAT(/10X,22(1H-)/10X,'HELLMANN-FEYNMAN FORCE'/10X,22(1H-))
 9010 FORMAT(1X,'HELFEY: SHELLS II,JJ=',2I5)
 9100 FORMAT(1X,'IC',I3,' II, JJ, I',3I5,1P,4E15.7/
     *       25X,'DE(1,2,3 AND IC)',1P,3E15.7)
      END
c--------------------------------------------------------------------------
C*MODULE GRD1    *DECK DVINT
      SUBROUTINE xDVINT
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D+00)
C
      COMMON /xDSTV  / XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,
     *                XJ,YJ,ZJ,NI,NJ,CX,CY,CZ
      COMMON /xHERMIT/ H1,H2(2),H3(3),H4(4),H5(5),H6(6),H7(7)
      COMMON /xWERMIT/ W1,W2(2),W3(3),W4(4),W5(5),W6(6),W7(7)
C
      DIMENSION H(28),W(28),MIN(7),MAX(7)
C
      EQUIVALENCE (H(1),H1),(W(1),W1)
C
      DATA MIN /1,2,4,7,11,16,22/
      DATA MAX /1,3,6,10,15,21,28/
C
C     ----- GAUSS-HERMITE QUADRATURE USING MINIMUM POINT FORMULA -----
C
      XINT = ZERO
      YINT = ZERO
      ZINT = ZERO
      NPTS = (NI+NJ+1-2)/2+1
      IMIN = MIN(NPTS)
      IMAX = MAX(NPTS)
      DO 340 I = IMIN,IMAX
         DUM = H(I)*T
         PTX = DUM+X0
         PTY = DUM+Y0
         PTZ = DUM+Z0
         PX = PTX-CX
         PY = PTY-CY
         PZ = PTZ-CZ
         AX = PTX-XI
         AY = PTY-YI
         AZ = PTZ-ZI
         BX = PTX-XJ
         BY = PTY-YJ
         BZ = PTZ-ZJ
         GO TO (180,170,160,140,120,100),NI
C
  100    PX = PX*AX
         PY = PY*AY
         PZ = PZ*AZ
  120    PX = PX*AX
         PY = PY*AY
         PZ = PZ*AZ
  140    PX = PX*AX
         PY = PY*AY
         PZ = PZ*AZ
  160    PX = PX*AX
         PY = PY*AY
         PZ = PZ*AZ
  170    PX = PX*AX
         PY = PY*AY
         PZ = PZ*AZ
  180    GO TO (320,300,280,260,240,220,200),NJ
C
  200    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  220    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  240    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  260    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  280    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  300    PX = PX*BX
         PY = PY*BY
         PZ = PZ*BZ
  320    DUM = W(I)
         XINT = XINT+DUM*PX
         YINT = YINT+DUM*PY
         ZINT = ZINT+DUM*PZ
  340 CONTINUE
      RETURN
      END
c--------------------------------------------------------------------------
C*MODULE GRD1    *DECK DENDD1
      SUBROUTINE xDENDD1(DA,DB,L2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /xIOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /xPSILVL/ IPSI
      COMMON /xWFNOPT/ SCFTYP,CITYP,DFTYPE,CCTYP,MPLEVL,MPCTYP
C
      DIMENSION DA(L2),DB(L2)
C
      DATA RHF,RMC/8HRHF     ,8HMCSCF   /
      DATA RNONE/8HNONE    /
C
C        READ TOTAL DENSITY MATRIX
C        -DENDD1- DIFFERS FROM -WFNDEN- IN THAT IT ADDS DA+DB.
C
      CALL yDAREAD(IDAF,IODA,DA,L2,16,0)
      IF(SCFTYP.EQ.RHF) RETURN
      IF(SCFTYP.EQ.RMC) RETURN
      IF(CITYP.NE.RNONE  .AND.  IPSI.EQ.1) RETURN
C
      CALL yDAREAD(IDAF,IODA,DB,L2,20,0)
C
C        RETURN TOTAL DENSITY IN -DA-
C
      CALL xVADD(DA,1,DB,1,DA,1,L2)
      RETURN
      END
c--------------------------------------------------------------------------
C*MODULE GRD1    *DECK EIJDEN
      SUBROUTINE xEIJDEN(EPS,V,E,IA,WRK,L1,L2,L3,L0,ECI)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION METHOD
C
      LOGICAL CANONC,FCORE,FORS,NOCI,EKT,LINSER,MFRZ
C
#include "gamess.par"
      PARAMETER (MXRT=100, MXNORO=250)
C
      DIMENSION V(L1,L0),EPS(L2),E(*),IA(L1),WRK(L1),ECI(L0,L0)
C
C        E(L1) FOR RHF, E(L3) FOR GVB/MCSCF/CI
C        NOTE FOR CI RUNS, -ECI- AND -E- MUST BE EQUIVALENT BY CALL
C
      COMMON /xDETWFN/ WSTATE(MXRT),SPINS(MXRT),CRIT,PRTTOL,SDET,SZDET,
     *                GRPDET,STSYM,GLIST,
     *                NFLGDM(MXRT),IWTS(MXRT),NCORSV,NCOR,NACT,NORBDT,
     *                NADET,NBDET,KDET,KSTDET,IROOT,IPURES,MAXW1,NITDET,
     *                MAXP,NCIDET,IGPDET,KSTSYM,NFTGCI
      COMMON /xGUGWFN/ NFZC,NMCC,NDOC,NAOS,NBOS,NALP,NVAL,NEXT,NFZV,
     *                IFORS,IEXCIT,ICICI,NOIRR
      COMMON /xGVBWFN/ CICOEF(2,12),F(25),ALPHA(325),BETA(325),NO(10),
     *                NCO,NSETO,NOPEN,NPAIR,NMOGVB,NCONF(MXAO),NHAM
      COMMON /xIOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /xINFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM)
      COMMON /xMCINP / METHOD,CISTEP,ACURCY,ENGTOL,DAMP,MICIT,NWORD,
     *                NORBMC,NOROT(2,MXNORO),MOFRZ(15),NPFLG(10),
     *                NOFO,CANONC,FCORE,FORS,NOCI,EKT,LINSER
      COMMON /xMFRPAR/ MFRZ,NUMFRZ,IFRZ(MXAO)
      COMMON /xWFNOPT/ SCFTYP,CITYP,DFTYPE,CCTYP,MPLEVL,MPCTYP
C
      PARAMETER (ZERO=0.0D+00, HALF=0.5D+00)
C
      DATA RHF,UHF,ROHF/8HRHF     ,8HUHF     ,8HROHF    /
      DATA  GVB,RMC/8HGVB     ,8HMCSCF   /
      DATA ALDET,GUGA,GENCI,ORMAS
     *      /8HALDET   ,8HGUGA    ,8HGENCI   ,8HORMAS   /
      DATA CIS/8HCIS     /
C
C     ----- COMPUTE "ENERGY WEIGHTED DENSITY MATRIX" -----
C     THIS QUANTITY IS ACTUALLY THE LAGRANGIAN MATRIX,
C     BACKTRANSFORMED INTO THE AO BASIS.
C
C           RHF CASE
C
      IF(SCFTYP.EQ.RHF .AND. MPLEVL.EQ.0 .AND.
     *   CITYP.NE.GUGA .AND. CITYP.NE.CIS) THEN
         IF(MFRZ) THEN
           L2 = (L1*L1+L1)/2
           CALL yDAREAD(IDAF,IODA,EPS,L2,36,0)
           RETURN
         END IF
         CALL yDAREAD(IDAF,IODA,V,L3,15,0)
         CALL yDAREAD(IDAF,IODA,E,L1,17,0)
         IJ = 0
         DO 130 I = 1,L1
            DO 120 J = 1,I
               IJ = IJ+1
               DUM = ZERO
               DO 100 K = 1,NA
                  DUM = DUM-E(K)*V(I,K)*V(J,K)
c                  write (96,*) 'i,j,k,v(i,k),v(j,k),e(k),dum'
c                  write (96,'(3i3,4g17.8)') i,j,k,v(i,k),v(j,k),e(k),dum
  100          CONTINUE
               EPS(IJ) = DUM+DUM
  120       CONTINUE
  130    CONTINUE
         RETURN
      END IF
C
C           UHF AND ROHF CASE
C
      IF((SCFTYP.EQ.UHF .OR. SCFTYP.EQ.ROHF)  .AND.  MPLEVL.EQ.0) THEN
         L2 = (L1*L1+L1)/2
         CALL yDAREAD(IDAF,IODA,EPS,L2,36,0)
         IJ = 0
         DO 160 I = 1,L1
            IJ = IJ+I
            EPS(IJ) = HALF*EPS(IJ)
  160    CONTINUE
         RETURN
      END IF
C
C     CLOSED SHELL MP2 CASE: JUST READ IT FROM DAF
C
      IF(SCFTYP.EQ.RHF .AND. MPLEVL.EQ.2) THEN
         L2 = (L1*L1+L1)/2
         CALL yDAREAD(IDAF,IODA,EPS,L2,309,0)
         RETURN
      END IF
C
C     CLOSED SHELL CIS CASE: JUST READ IT FROM DAF
C
cwc      IF(SCFTYP.EQ.RHF .AND. CITYP.EQ.CIS) THEN
cwc         L2 = (L1*L1+L1)/2
cwc         CALL DAREAD(IDAF,IODA,EPS,L2,309,0)
cwc         RETURN
cwc      END IF
C
C     UMP2 CASE: READ ALPHA AND BETA ENERGY WEIGHTED DENSITY, AND ADD
C
      IF(SCFTYP.EQ.UHF .AND. MPLEVL.EQ.2) THEN
         L2 = (L1*L1+L1)/2
         CALL yDAREAD(IDAF,IODA,EPS,L2,309,0)
         CALL yDAREAD(IDAF,IODA,E  ,L2,369,0)
         DO N=1,L2
            EPS(N)=EPS(N)+E(N)
         ENDDO
         RETURN
      END IF
C
C           ONLY GVB/MCSCF/CI GET THIS FAR...
C
cwc      IF(SCFTYP.EQ.GVB) NORB = NMOGVB
cwc      IF(SCFTYP.EQ.RMC) THEN
cwc         IF(CISTEP.EQ.ALDET) NORB = NCORSV+NACT
cwc         IF(CISTEP.EQ.GENCI) NORB = NCORSV+NACT
cwc         IF(CISTEP.EQ.ORMAS) NORB = NCORSV+NACT
cwc         IF(CISTEP.EQ.GUGA)  NORB = NMCC+NDOC+NAOS+NBOS+NALP+NVAL
cwc      END IF
cwc      IF(CITYP.EQ.GUGA) NORB = NFZC+NDOC+NAOS+NBOS+NALP+NVAL
cwc      DO 210 I = 1,L1
cwc         IA(I) = (I*I-I)/2
cwc  210 CONTINUE
cwcC
cwcC     ----- READ MOS AND LAGRANGIAN FROM DISK -----
cwcC
cwc      L3ORB = NORB*NORB
cwc      IF(SCFTYP.EQ.RMC  .OR.  CITYP.EQ.GUGA) L3ORB = L3
cwc      CALL DAREAD(IDAF,IODA,V,L3,15,0)
cwc      CALL DAREAD(IDAF,IODA,E,L3ORB,36,0)
cwcC
cwcC     ----- MAP MCSCF/CI LAGRANGIAN INTO NORB BY NORB FORMAT -----
cwcC
cwc      IF(SCFTYP.EQ.GVB .OR. CITYP.EQ.GUGA) GO TO 280
cwc      IF(NORB.GT.1) THEN
cwc         IADD = 1
cwc         IADDE = 1
cwc         DO 260 I = 1,NORB
cwc            IADD = IADD + NORB
cwc            IADDE = IADDE + L1
cwc            CALL XCOPY(NORB,E(IADDE),1,E(IADD),1)
cwc  260    CONTINUE
cwc      END IF
cwcC
cwcC     ----- CALCULATE C-DAGGER*E*C -----
cwcC
cwc  280 CONTINUE
cwc      CALL VCLR(EPS,1,(L1*L1+L1)/2)
cwcC
cwc      IF (CITYP.NE.GUGA) THEN
cwcC
cwcC     ---CALCULATE THE HALF TRANSFORM FIRST -----
cwcC     NOTE THAT E(KL) IS USED EXACTLY TWICE.  DIVIDE BY TWO TO GET
cwcC     THE VALUES APPROPRIATE FOR THE GENERALIZED LAGRANGE MULTIPLIERS.
cwcC
cwc      DO 400 I = 1,L1
cwc         KL = 0
cwc         DO 330 L = 1,NORB
cwc            WRK(L) = ZERO
cwc            DO 320 K = 1,NORB
cwc               KL = KL + 1
cwc               WRK(L) = WRK(L) - V(I,K)*E(KL)
cwc  320       CONTINUE
cwc  330    CONTINUE
cwc         DO 340 L = 1,NORB
cwc            WRK(L) = WRK(L)*HALF
cwc  340    CONTINUE
cwcC
cwc         DO 380 J = 1,L1
cwc            IJ = IA(I) + J
cwc            IF(J .GT. I) IJ = IA(J) + I
cwc            DO 360 L = 1,NORB
cwc               EPS(IJ) = EPS(IJ) + WRK(L)*V(J,L)
cwc  360       CONTINUE
cwc  380    CONTINUE
cwc  400 CONTINUE
cwcC
cwc      ELSE
cwcC
cwcC     --- FOR CI THE LAGRANGIAN -ECI- IS NOT SYMMETRIC ---
cwcC
cwc      DO 500 I = 1,L1
cwc         DO 430 L = 1,NORB
cwc            WRK(L) = ZERO
cwc            DO 420 K = 1,NORB
cwc               IF(K.LE.L) THEN
cwc                  WRK(L) = WRK(L) - V(I,K)*ECI(K,L)
cwc               ELSE
cwc                  WRK(L) = WRK(L) - V(I,K)*ECI(L,K)
cwc               END IF
cwc  420       CONTINUE
cwc  430    CONTINUE
cwc         DO 440 L = 1,NORB
cwc            WRK(L) = WRK(L)*HALF
cwc  440    CONTINUE
cwcC
cwc         DO 480 J = 1,L1
cwc            IJ = IA(I) + J
cwc            IF(J .GT. I) IJ = IA(J) + I
cwc            DO 460 L = 1,NORB
cwc               EPS(IJ) = EPS(IJ) + WRK(L)*V(J,L)
cwc  460       CONTINUE
cwc  480    CONTINUE
cwc  500 CONTINUE
cwcC
cwc      END IF
cwcC
cwc      RETURN
      END
c--------------------------------------------------------------------------
C*MODULE GRD1    *DECK VNNDER
      SUBROUTINE xVNNDER(DRG,NATM,EF3,NFRPTS,COSDRG,NDIMCD,WORK)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
cwc1
#undef _PARALLEL_ACES_COM_
#include "parallel_aces.com"
cwc0
C
      LOGICAL OUT,GOPARR,DSKWRK,MASWRK,EFLDL
C
#include "gamess.par"
      PARAMETER (MXDFG=5, MXPT=100, MXFRG=50,
     *           MXFGPT=MXPT*MXFRG)
      PARAMETER (MXPSH=5*MXPT, MXPG=5*MXPSH)
C
      DIMENSION DRG(NATM,NATM),EF3(3,*),COSDRG(NDIMCD,NDIMCD),
     *          WORK(10,*)
      DIMENSION EF3T(3)
C
cwc      COMMON /xEFLDC / EVEC(3),EFLDL
cwc      COMMON /xEFPOTD/ METHOD,INABIO,MOVE
cwc      COMMON /xFGRAD / DEF(3,MXFGPT),DEFT(3,MXFRG),TORQ(3,MXFRG),
cwc     *                EFCENT(3,MXFRG),FRGMAS(MXFRG),FRGMI(6,MXFRG),
cwc     *                ATORQ(3,MXFRG)
      COMMON /xFMCOM / XX(1)
cwc      COMMON /xFRGINF/ NMPTS(MXFRG),NMTTPT,IEFC,IEFD,IEFQ,IEFO,
cwc     *                NPPTS(MXFRG),NPTTPT,IEFP,
cwc     *                NRPTS(MXFRG),NRTTPT,IREP,ICHGP,NFRG
cwc      COMMON /xFRGTYP/ NDFRG,ISET(MXFRG),NAMIDX(MXDFG)
      COMMON /xGRAD  / DE(3,MXATM)
      COMMON /xINFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM)
      COMMON /xIOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /xOUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /xPAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /xPAULIN/ EX(MXPG,MXFRG),CS(MXPG,MXFRG),CP(MXPG,MXFRG),
     *                CD(MXPG,MXFRG),CF(MXPG,MXFRG),CG(MXPG,MXFRG),
     *                PRNAME(MXFGPT),PRCORD(3,MXFGPT),EFZNUC(MXFGPT),
     *                KSTART(MXPSH,MXFRG),KATOM(MXPSH,MXFRG),
     *                KTYPE(MXPSH,MXFRG),KNG(MXPSH,MXFRG),
     *                KLOC(MXPSH,MXFRG), KMIN(MXPSH,MXFRG),
     *                KMAX(MXPSH,MXFRG),NSHELL(MXFRG),NGAUSS(MXFRG),
     *                NATEF(MXFRG),NUMEF(MXFRG),NTPATM
C
cwc      COMMON /xPAULMO/ CENTNM(MXFGPT),CENTCD(3,MXFGPT),NORB(MXFRG),
cwc     *                NPBF(MXFRG),NTMO
cwc      COMMON /xPCMPAR/ IPCM,NFT26,NFT27,IKREP,IEF,IP_F
      COMMON /xRUNOPT/ RUNTYP,EXETYP,NEVALS
cwc      COMMON /xZRFPAR/ GZRF,FIND(3),GNUCF,EBORN,DIELEC,IZRF,ICALC
C
C     COSMO CHANGES
C
      PARAMETER (MAXDEN=25*MXATM, LENABC=2000, NPPA=1082)
      LOGICAL ISEPS,USEPS
cwc      COMMON /xISEPS / ISEPS, USEPS
      COMMON /xSOLV  / FEPSI,RDS,DISEX2,COSURF(3,LENABC),SRAD(MXATM),
     *                QDEN(MAXDEN),AR(LENABC),
     *                NSPA,NPS,NPS2,NDEN,NPSPHER
      COMMON /xSOLVI / IATSP(LENABC+1),N0(2),NP1,NP2,ISKIP
      COMMON /xCOSMO3/ COSZAN(NPPA),CORZAN(3,NPPA)
C
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00)
C
      DATA DEBUG/8HDEBUG   /, DBUGME/8HVNNDER  /, GRD1/8HGRD1    /
C
C     ----- GRADIENT OF NUCLEAR REPULSION ENERGY -----
C     ----- PLUS NUCLEAR CONTRIBUTIONS FROM EFFECTIVE FRAGMENTS -----
C
      OUT = (EXETYP.EQ.DEBUG  .OR.  EXETYP.EQ.DBUGME  .OR.
     *       EXETYP.EQ.GRD1   .OR.  NPRINT.EQ.-3)  .AND.  MASWRK
C
C     ----- FORM DISTANCE MATRIX -----
C
      DRG(1,1) = ZERO
      DO 130 K = 2,NAT
         DRG(K,K) = ZERO
         K1 = K-1
         DO 120 L = 1,K1
            RKL = ZERO
            DO 100 I = 1,3
               RKL = RKL+(C(I,K)-C(I,L))**2
  100       CONTINUE
            DRG(K,L) = -ONE/RKL
            DRG(L,K) = SQRT(RKL)
  120    CONTINUE
  130 CONTINUE
C
C     ----- FORM CONTRIBUTION TO GRADIENT -----
C
      DO 290 KK = 1,3
         DO 240 K = 2,NAT
            ZAK = ZAN(K)
            KM1 = K-1
            DO 230 L = 1,KM1
               ZAL = ZAN(L)
               PKL = (C(KK,K)-C(KK,L))/DRG(L,K)
               DE(KK,K) = DE(KK,K)+PKL*DRG(K,L)*ZAK*ZAL
  230       CONTINUE
  240    CONTINUE
C
         NAT1 = NAT-1
         DO 280 K = 1,NAT1
            ZAK = ZAN(K)
            KP1 = K+1
            DO 270 L = KP1,NAT
               ZAL = ZAN(L)
               PKL = (C(KK,K)-C(KK,L))/DRG(K,L)
               DE(KK,K) = DE(KK,K)+PKL*DRG(L,K)*ZAK*ZAL
  270      CONTINUE
  280    CONTINUE
  290 CONTINUE
C
C     COSMO CONTRIBUTION - KIM, PUT IN BY VOLKER OCT.1998
C
cwc      IF(ISEPS) THEN
cwc         WRITE(IW,*)"NOW CALCULATING SAS-NUC GRADIENT CONTRIBUTION:"
cwc         DO 350 IZ=1,NDIMCD
cwc            DO 355 JZ=1,NDIMCD
cwc               COSDRG(IZ,JZ)=ZERO
cwc  355       CONTINUE
cwc  350    CONTINUE
cwcC
cwc         DO 292 K=1,NAT
cwc            DO 294 J=1+NAT,NPS+NAT
cwcC              IAJ=IATSP(J-NAT)
cwc               RKJ=ZERO
cwc               DO 296  I=1,3
cwc                  RKJ=RKJ+(C(I,K)-CORZAN(I,J-NAT))**2
cwc  296          CONTINUE
cwc               COSDRG(K,J)=-1/RKJ
cwc               COSDRG(J,K)=SQRT(RKJ)
cwc  294       CONTINUE
cwc  292    CONTINUE
cwcC
cwc      DO 297 KK=1,3
cwc         DO  298 K=1,NAT
cwc            ZAK=ZAN(K)
cwc            DO 299 L=1+NAT,NPS+NAT
cwc               IAL=IATSP(L-NAT)
cwc               ZAL=COSZAN(L-NAT)
cwc               PKL=(C(KK,K)-CORZAN(KK,L-NAT))/COSDRG(L,K)
cwc               DE(KK,K)=DE(KK,K)+PKL*COSDRG(K,L)*ZAK*ZAL
cwc               DE(KK,IAL)=DE(KK,IAL)-PKL*COSDRG(K,L)*ZAK*ZAL
cwc  299       CONTINUE
cwc  298    CONTINUE
cwc  297 CONTINUE
cwcC
cwc      DO 360 I=1,NPS
cwcC        IAI=IATSP(I)
cwc         DO 362 J=1,I-1
cwcC           IAJ=IATSP(J)
cwcC            IF(IAI.EQ.IAJ) GOTO 362
cwc            RIJ=ZERO
cwc            DO 364 K=1,3
cwc               RIJ=RIJ+(CORZAN(K,I)-CORZAN(K,J))**2
cwc  364       CONTINUE
cwc            COSDRG(I+NAT,J+NAT)=-1/RIJ
cwc            COSDRG(J+NAT,I+NAT)=SQRT(RIJ)
cwc  362    CONTINUE
cwc  360  CONTINUE
cwcC
cwc      DO 366 KK=1,3
cwc         DO 368 K=1,NPS
cwc            IAK=IATSP(K)
cwc            ZAK=COSZAN(K)
cwc            DO 370 L=1,K-1
cwc               IAL=IATSP(L)
cwc               ZAL=COSZAN(L)
cwc               IF(IAK.EQ.IAL) GOTO 370
cwc               PKL=(CORZAN(KK,K)-CORZAN(KK,L))/COSDRG(L+NAT,K+NAT)
cwc              DE(KK,IAK)=DE(KK,IAK)+PKL*COSDRG(K+NAT,L+NAT)*ZAK*ZAL
cwc              DE(KK,IAL)=DE(KK,IAL)-PKL*COSDRG(K+NAT,L+NAT)*ZAK*ZAL
cwc  370      CONTINUE
cwc  368    CONTINUE
cwc  366  CONTINUE
cwc      ENDIF
C
C     ----- NUCLEAR CONTRIBUTIONS FROM EXTERNAL ELECTRIC FIELD -----
C
cwc      IF (EFLDL) CALL DNUCEF(DE)
C
C     ----- NUCLEAR CONTRIBUTIONS FROM ZERNER REACTION FIELD -----
C
cwc      IF(IZRF.EQ.1) CALL DNUCZ(DE)
C
C     ----- NUCLEAR CONTRIBUTIONS FROM EFFECTIVE FRAGMENTS -----
C
cwc      CALL VCLR(EF3,1,3*NFRPTS)
cwc      IF(IEFC.EQ.1) CALL DNUCC(DE,EF3)
cwc      IF(IEFD.EQ.1) CALL DNUCD(DE,EF3)
cwc      IF(IEFQ.EQ.1) CALL DNUCQ(DE,EF3)
cwc      IF(IEFO.EQ.1) CALL DNUCO(DE,EF3)
cwc      IF(IEFP.EQ.1) CALL DNUCP(DE,EF3)
cwc      IF(IP_F.EQ.1) CALL IEFFR(EF3)
cwcC
cwc      IF(NFRG.GT.1) THEN
cwc         IF(IEFC.EQ.1) CALL DCHCH(EF3)
cwc         IF(IEFD.EQ.1) CALL DCHDIP(EF3)
cwc         IF(IEFD.EQ.1) CALL DDPDP(EF3)
cwc         IF(IEFQ.EQ.1) CALL DCHQUA(EF3,WORK)
cwc         IF(IEFQ.EQ.1) CALL DDPQUD(EF3,WORK)
cwc         IF(IEFQ.EQ.1) CALL DQDQD(EF3,WORK)
cwc         IF(IEFO.EQ.1) CALL DCHOCT(EF3,WORK)
cwc         IF(IEFP.EQ.1.AND.IEFC.EQ.1) CALL DCHIND(EF3)
cwc         IF(IEFP.EQ.1) CALL DININ(EF3)
cwc         IF(IEFP.EQ.1.AND.IEFD.EQ.1) CALL DDPIND(EF3)
cwc         IF(IEFP.EQ.1.AND.IEFQ.EQ.1) CALL DQDIND(EF3,WORK)
cwcC
cwc         IF (METHOD .GE. 4) GO TO 400
cwcC
cwc         MXRPTS=0
cwc         DO 300 IFRG=1,NFRG
cwc            MXRPTS = MAX(NRPTS(IFRG),MXRPTS)
cwc  300    CONTINUE
cwc         MXRPTS=MXRPTS+1
cwc         NDFRG2 = (NDFRG*NDFRG+NDFRG)/2
cwc         LEN = MXRPTS*MXRPTS*NDFRG2
cwc         CALL VALFM(LOADFM)
cwc         LAREL = LOADFM + 1
cwc         LCREL = LAREL  + LEN
cwc         LAST  = LCREL  + LEN
cwc         NEED = LAST - LOADFM - 1
cwc         CALL GETFM(NEED)
cwc         IF(NTMO.EQ.0) CALL DAREAD(IDAF,IODA,XX(LAREL),2*LEN,75,0)
cwc         IF(IREP.EQ.1) CALL DREPUL(EF3,XX(LAREL),XX(LCREL),
cwc     *                             MXRPTS,NDFRG2)
cwc         CALL RETFM(NEED)
cwcC
cwc         IF(NTMO.EQ.0) GO TO 390
cwc         MXMO=0
cwc         MXBF=0
cwc         DO I = 1, NFRG
cwc            MXMO=MAX(MXMO,NORB(I))
cwc            MXBF=MAX(MXBF,NPBF(I))
cwc         END DO
cwc         MXMO2=(MXMO*MXMO+MXMO)/2
cwcC
cwc         CALL VALFM(LOADFM)
cwc         LPROVEC = LOADFM  + 1
cwc         LFOCKMA = LPROVEC + MXBF*NTMO
cwc         LSMAT   = LFOCKMA + MXMO2*NFRG
cwc         LTMAT   = LSMAT   + MXBF*MXBF
cwc         LEPS    = LTMAT   + MXBF*MXBF
cwc         LEPT    = LEPS    + MXBF*MXBF
cwc         LWRK    = LEPT    + MXBF*MXBF
cwc         LSIJ    = LWRK    + MXBF
cwc         LTIJ    = LSIJ    + MXMO*MXMO
cwc         LFASQ   = LTIJ    + MXMO*MXMO
cwc         LFBSQ   = LFASQ   + MXMO*MXMO
cwc         LAST    = LFBSQ   + MXMO*MXMO
cwc         NEED=LAST-LOADFM-1
cwc         CALL GETFM(NEED)
cwcC
cwc         LENPV=MXBF*NTMO
cwc         LENFM=MXMO2*NFRG
cwc         CALL DAREAD(IDAF,IODA,XX(LPROVEC),LENPV,267,0)
cwc         CALL DAREAD(IDAF,IODA,XX(LFOCKMA),LENFM,268,0)
cwc         CALL DPAULR(EF3,NFRPTS,
cwc     *               XX(LPROVEC),XX(LFOCKMA),XX(LSMAT),XX(LTMAT),
cwc     *               XX(LEPS),XX(LEPT),XX(LWRK),XX(LSIJ),XX(LTIJ),
cwc     *               XX(LFASQ),XX(LFBSQ),MXBF,MXMO,MXMO2)
cwc         CALL RETFM(NEED)
cwc  390    CONTINUE
cwc      END IF
cwcC
cwcC     ----- FRAGMENT TORQUES -----
cwcC
cwc  400 CONTINUE
cwc      DO 430 I = 1,NMTTPT+NPTTPT+NRTTPT+NTPATM+NTMO
cwc         DO 420 J = 1,3
cwc            DEF(J,I) = DEF(J,I) + EF3(J,I)
cwcC          IF(MASWRK)WRITE(6,*)'I,J,DEF(J,I)',I,J,DEF(J,I)
cwc  420    CONTINUE
cwc  430 CONTINUE
cwc      DO 440 J=1,3
cwc         EF3T(J)=ZERO
cwc  440 CONTINUE
cwc      DO 480 I=1,NMTTPT+NPTTPT+NRTTPT+NTPATM+NTMO
cwc         DO 470 J=1,3
cwc            EF3T(J)=EF3T(J)+EF3(J,I)
cwc  470    CONTINUE
cwc  480 CONTINUE
C
      IF(OUT.AND.MASWRK) THEN
cwc         WRITE(IW,9088)
cwc         CALL EGOUT(DE,NAT)
cwc         CALL EFOUT(EF3,EF3T,EF3T)
      END IF
cwc1
c         WRITE(96,9088)
c         CALL yEGOUT(DE,NAT)
cwc0
C
C     ----- DIVIDE BY NUMBER OF NODES -----
C     ALL NODES HAVE COMPUTED THE SAME QUANTITY UP TO THIS POINT,
C     SO IT MUST BE SCALED DOWN BEFORE THE EVENTUAL GLOBAL SUM.
C
cwc      IF(GOPARR) CALL XSCAL(3*NAT,ONE/NPROC,DE,1)
cwc1
      call XSCAL (3*NAT,one/nprocs,DE,1)
cwc0
      RETURN
C
 9088 FORMAT(/10X,36(1H-)/10X,'GRADIENT OF NUCLEAR REPULSION ENERGY'/
     *        10X,36(1H-))
      END
c--------------------------------------------------------------------------
C*MODULE GRD1    *DECK SDER
      SUBROUTINE xSDER(EPS,DSINT,L1,L2,NOCON)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
cwc1
#undef _PARALLEL_ACES_COM_
#include "parallel_aces.com"
cwc0
C
      LOGICAL DBG,NORM,OUT
      LOGICAL GOPARR,DSKWRK,MASWRK,NXT,SVDSKW,NOCON
C
#include "gamess.par"
C
      DIMENSION EPS(L2),DSINT(L2,3,NAT)
      DIMENSION DIJ(225), DIJINT(225), IJX(35), IJY(35),IJZ(35),
     *           XS(6,5), YS(6,5), ZS(6,5),
     *          DXS(5,5),DYS(5,5),DZS(5,5)
C
      COMMON /xDSTV  / XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,
     *                XJ,YJ,ZJ,NI,NJ,CX,CY,CZ
      COMMON /xGRAD  / DE(3,MXATM)
      COMMON /xIJPAIR/ IA(MXAO)
      COMMON /xINFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM)
      COMMON /xIOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /xOUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /xPAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /xRUNOPT/ RUNTYP,EXETYP,NEVALS
cwc      COMMON /xRELWFN/ RMETHOD,QRQMT,CLIG,CLIG2,QRTOL,IQRORD,MODQR,NESOC,
cwc     *                NRATOM,NUMU,NQMTR,NQRDAF,MORDA,NDARELB
      COMMON /xFMCOM / X(1)
C
      PARAMETER (SQRT3=1.73205080756888D+00)
      PARAMETER (SQRT5=2.23606797749979D+00, SQRT7=2.64575131106459D+00)
      PARAMETER (ONE=1.0D+00)
      PARAMETER (RLN10=2.30258D+00)
C
      DATA DEBUG/8HDEBUG   /, DBUGME/8HSDER    /, GRD1/8HGRD1    /,
     *     RESC/8HRESC    /
C
      DATA IJX / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     *           4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     *           5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     *           3, 1, 3, 2, 2/
      DATA IJY / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     *           1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     *           1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     *           1, 3, 2, 3, 2/
      DATA IJZ / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     *           1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     *           1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     *           3, 3, 2, 2, 3/
C
C     ----- CALCULATE DERIVATIVES OF THE OVERLAP MATRIX -----
C     THIS TERM RESULTS FROM DERIVATIVES OF THE DENSITY MATRIX,
C     WHICH CAN IN THE CASE OF A SYMMETRIC LAGRANGIAN MATRIX, CAN
C     BE CONVERTED INSTEAD INTO DERIVATIVES OF OVERLAP INTEGRALS.
C     RESC RUNS NEED BOTH CONTRACTION OF THE INTEGRALS WITH DENSITY
C     AND THE INTEGRALS THEMSELVES.
C     IN CASE OF INTERNALLY UNCONTRACTED RESC, ONE HAS TO TRANSFORM
C     THE OVERLAP DERIVATIVES TO CONTRACTED BASIS BEFORE
C     SUMMING (CONTRACTING) OVERLAP DERIVATIVES WITH THE LAGRANGIAN
C
C     NOCON IS TRUE IF LAGRANGIAN IS NOT TO BE SUMMED WITH.
C
c      write (96,*) 'xsder'
      DBG = (EXETYP.EQ.DEBUG  .OR.  EXETYP.EQ.DBUGME)  .AND.  MASWRK
      OUT = (EXETYP.EQ.DEBUG  .OR.  EXETYP.EQ.DBUGME  .OR.
     *       EXETYP.EQ.GRD1   .OR.  NPRINT.EQ.-3)
cwc      IF(DBG) WRITE(IW,9000)
cwc1
c      write (96,9000)
cwc0
C
      TOL = RLN10*ITOL
      NORM = NORMF .NE. 1 .OR. NORMP .NE. 1
cwc      IF(RMETHOD.EQ.RESC) CALL VCLR(DSINT,1,L2*3*NAT)
C
C     INITIALIZE PARALLEL
C
      NXT = IBTYP.EQ.1
      IPCOUNT = ME - 1
      NEXT = -1
      MINE = -1
cwc1
      nloop=-1
cwc0
C
      IF(DBG) THEN
cwc         WRITE(IW,9010)
cwc         CALL PRTRI(EPS,L1)
      END IF
cwc1
c         WRITE(96,9010)
c         CALL yPRTRI(EPS,L1)
cwc0
C
C     ----- I SHELL
C
      DO 780 II = 1,NSHELL
C
C           GO PARALLEL!
C
      IF(NXT .AND. GOPARR) THEN
         MINE = MINE + 1
cwc         IF(MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
cwc         IF(NEXT.NE.MINE) GO TO 780
      END IF
C
      IAT = KATOM(II)
      XI = C(1,IAT)
      YI = C(2,IAT)
      ZI = C(3,IAT)
      I1 = KSTART(II)
      I2 = I1+KNG(II)-1
      LIT = KTYPE(II)
      MINI = KMIN(II)
      MAXI = KMAX(II)
      LOCI = KLOC(II)-MINI
      LITDER = LIT+1
C
C     ----- J SHELL
C
      DO 760 JJ = 1,II
C
C           GO PARALLEL!
C
cwc1
cwc      IF((.NOT.NXT) .AND. GOPARR) THEN
cwc         IPCOUNT = IPCOUNT + 1
cwc         IF(MOD(IPCOUNT,NPROC).NE.0) GO TO 760
cwc      END IF
              nloop=nloop+1
              nloop=mod(nloop,nprocs)
              if (nloop.ne.irank) goto 760
cwc0
C
      JAT = KATOM(JJ)
      XJ = C(1,JAT)
      YJ = C(2,JAT)
      ZJ = C(3,JAT)
      J1 = KSTART(JJ)
      J2 = J1+KNG(JJ)-1
      LJT = KTYPE(JJ)
      MINJ = KMIN(JJ)
      MAXJ = KMAX(JJ)
      LOCJ = KLOC(JJ)-MINJ
      RR = (XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
      IF(II.EQ.JJ) GO TO 760
C
C     ----- I PRIMITIVE
C
      DO 640 IG = I1,I2
        AI = EX(IG)
        ARRI = AI*RR
        AXI = AI*XI
        AYI = AI*YI
        AZI = AI*ZI
        CSI=CS(IG)
        CPI=CP(IG)
        CDI=CD(IG)
        CFI=CF(IG)
        CGI=CG(IG)
C
C     ----- J PRIMITIVE
C
        DO 620 JG = J1,J2
          AJ = EX(JG)
          AA = AI+AJ
          AA1 = ONE/AA
          DUM = AJ*ARRI*AA1
          IF(DUM .GT. TOL) GO TO 620
          FAC = EXP(-DUM)
          CSJ = CS(JG)
          CPJ = CP(JG)
          CDJ = CD(JG)
          CFJ = CF(JG)
          CGJ = CG(JG)
          AX = (AXI+AJ*XJ)*AA1
          AY = (AYI+AJ*YJ)*AA1
          AZ = (AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR
C
          IJ = 0
          DO 580 I = MINI,MAXI
            IF(I.EQ.1) DUM1=CSI*FAC
            IF(I.EQ.2) DUM1=CPI*FAC
            IF(I.EQ.5) DUM1=CDI*FAC
            IF(I.EQ.8.AND.NORM) DUM1=DUM1*SQRT3
            IF(I.EQ.11) DUM1=CFI*FAC
            IF(I.EQ.14.AND.NORM) DUM1=DUM1*SQRT5
            IF(I.EQ.20.AND.NORM) DUM1=DUM1*SQRT3
            IF(I.EQ.21) DUM1=CGI*FAC
            IF(I.EQ.24.AND.NORM) DUM1=DUM1*SQRT7
            IF(I.EQ.30.AND.NORM) DUM1=DUM1*SQRT5/SQRT3
            IF(I.EQ.33.AND.NORM) DUM1=DUM1*SQRT3
            DO 570 J = MINJ,MAXJ
              IF(J.EQ.1) DUM2=DUM1*CSJ
              IF(J.EQ.2) DUM2=DUM1*CPJ
              IF(J.EQ.5) DUM2=DUM1*CDJ
              IF(J.EQ.8.AND.NORM) DUM2=DUM2*SQRT3
              IF(J.EQ.11) DUM2=DUM1*CFJ
              IF(J.EQ.14.AND.NORM) DUM2=DUM2*SQRT5
              IF(J.EQ.20.AND.NORM) DUM2=DUM2*SQRT3
              IF(J.EQ.21) DUM2=DUM1*CGJ
              IF(J.EQ.24.AND.NORM) DUM2=DUM2*SQRT7
              IF(J.EQ.30.AND.NORM) DUM2=DUM2*SQRT5/SQRT3
              IF(J.EQ.33.AND.NORM) DUM2=DUM2*SQRT3
              IJ=IJ+1
              NN=IA(LOCI+I)+(LOCJ+J)
cwc              IF(RMETHOD.EQ.RESC) THEN
cwc                DIJINT(IJ)=DUM2
cwc              ENDIF
              DEN = EPS(NN)
              DEN=DEN+DEN
              DIJ(IJ) = DUM2*DEN
  570       CONTINUE
  580     CONTINUE
C
C     ----- OVERLAP
C
      T = SQRT(AA1)
      X0 = AX
      Y0 = AY
      Z0 = AZ
      DO 590 J = 1,LJT
        NJ = J
        DO 590 I = 1,LITDER
          NI = I
          CALL xVINT
          XS(I,J)=XINT*T
          YS(I,J)=YINT*T
          ZS(I,J)=ZINT*T
  590 CONTINUE
C
      CALL xDERI(DXS,DYS,DZS,XS,YS,ZS,LIT,LJT,AI)
C
      IJ=0
      DO 600 I=MINI,MAXI
        IX=IJX(I)
        IY=IJY(I)
        IZ=IJZ(I)
        DO 600 J=MINJ,MAXJ
          JX=IJX(J)
          JY=IJY(J)
          JZ=IJZ(J)
          DUMX=DXS(IX,JX)* YS(IY,JY)* ZS(IZ,JZ)
          DUMY= XS(IX,JX)*DYS(IY,JY)* ZS(IZ,JZ)
          DUMZ= XS(IX,JX)* YS(IY,JY)*DZS(IZ,JZ)
          IJ=IJ+1
cwc          IF(RMETHOD.EQ.RESC) THEN
cwc            NN=IA(LOCI+I)+(LOCJ+J)
cwc            DUMINT=DIJINT(IJ)
cwc            DSINT(NN,1,IAT)=DSINT(NN,1,IAT)+DUMX*DUMINT
cwc            DSINT(NN,2,IAT)=DSINT(NN,2,IAT)+DUMY*DUMINT
cwc            DSINT(NN,3,IAT)=DSINT(NN,3,IAT)+DUMZ*DUMINT
cwc            DSINT(NN,1,JAT)=DSINT(NN,1,JAT)-DUMX*DUMINT
cwc            DSINT(NN,2,JAT)=DSINT(NN,2,JAT)-DUMY*DUMINT
cwc            DSINT(NN,3,JAT)=DSINT(NN,3,JAT)-DUMZ*DUMINT
cwc          ENDIF
          IF(.NOT.NOCON) THEN
            DE(1,IAT)=DE(1,IAT)+(DUMX*DIJ(IJ))
            DE(2,IAT)=DE(2,IAT)+(DUMY*DIJ(IJ))
            DE(3,IAT)=DE(3,IAT)+(DUMZ*DIJ(IJ))
            DE(1,JAT)=DE(1,JAT)-(DUMX*DIJ(IJ))
            DE(2,JAT)=DE(2,JAT)-(DUMY*DIJ(IJ))
            DE(3,JAT)=DE(3,JAT)-(DUMZ*DIJ(IJ))
          ENDIF
 600  CONTINUE
C
 620  CONTINUE
 640  CONTINUE
C
C     ----- END OF PRIMITIVE LOOPS -----
C
      IF(DBG) THEN
cwc         WRITE(IW,9100) II,JJ
cwc         CALL EGOUT(DE,NAT)
      END IF
cwc1
c         WRITE(96,9100) II,JJ
c         CALL yEGOUT(DE,NAT)
cwc0
  760 CONTINUE
  780 CONTINUE
C
C     ----- END OF SHELL LOOPS -----
C
      IF(OUT) THEN
cwc         IF(GOPARR) CALL DDI_GSUMF(1501,DE,3*NAT)
cwc         IF(MASWRK) THEN
cwc            WRITE(IW,9000)
cwc            CALL EGOUT(DE,NAT)
cwc         END IF
cwc         IF(GOPARR) CALL XSCAL(3*NAT,ONE/NPROC,DE,1)
cwc      END IF
cwc      IF(GOPARR .AND. NXT) CALL DDI_DLBRESET
C
cwc      IF(RMETHOD.EQ.RESC) THEN
cwc        IF(GOPARR) CALL DDI_GSUMF(1501,DSINT,L2*3*NAT)
cwc        SVDSKW=DSKWRK
cwc        DSKWRK=.TRUE.
cwc        LREC=0
cwcC       STARTING RECORD FOR THE S DERIVATIVE INTEGRALS
cwc        DO I=1,NAT
cwc          DO J=1,3
cwc            LREC=LREC+1
cwc            CALL RAWRIT(NQRDAF,X(MORDA),DSINT(1,J,I),L2,LREC,0)
cwc          ENDDO
cwc        ENDDO
cwc        DSKWRK=SVDSKW
      ENDIF
C
      RETURN
C
 9000 FORMAT(/10X,33(1H-)/10X,'GRADIENT INCLUDING DENSITY FORCES'/
     *        10X,33(1H-))
 9010 FORMAT(1X,'THE LAGRANGIAN IN THE AO BASIS IS')
 9100 FORMAT(1X,'SDER: SHELLS II,JJ=',2I5)
      END
c--------------------------------------------------------------------------
C*MODULE MTHLIB  *DECK VADD
      SUBROUTINE xVADD(A,I,B,J,C,K,N)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION A(*),B(*),C(*)
C
C     ----- COMPUTE VECTOR C = A+B -----
C     I,J,K ARE INCREMENTS THROUGH A,B,C RESPECTIVELY
C
      IF (I .NE. 1) GO TO 200
      IF (J .NE. 1) GO TO 200
      IF (K .NE. 1) GO TO 200
C
C        I = J = K = 1
C
      DO 110 L=1,N
         C(L) = A(L) + B(L)
  110 CONTINUE
      RETURN
C
  200 CONTINUE
      LA=1-I
      LB=1-J
      LC=1-K
      DO 210 L=1,N
         LA=LA+I
         LB=LB+J
         LC=LC+K
         C(LC) = A(LA) + B(LB)
  210 CONTINUE
      RETURN
      END
c--------------------------------------------------------------------------
