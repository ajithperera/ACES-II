c one-electron Gamess gradient integrals
c before using this routine, GAMDER_ini must be called
c
c on entry:
c
c nL2 = (n*n+n)/2 where n is number of orbitals
c Sder,Hder: arrays, declared at least nL2*3*nat (nat=no.of atoms)
c DE: array, declared at least 3*nat
c contr: if .true., calculate both integrals and 1-el gradient
c        if .false., calculate only integrals and nuclear repul grad
c dens: total density matrix, full triangle (only needed if contr=T)
c edens: total Lagrangian matrix, full triangle (only needed if contr=T)
c
c on exit:
c
c Sder: S derivatives; Hder: H (core Hamiltonian) derivatives
c both arrays are (nL2,3,nat);  3 is for x,y,z
c DE(3,nat): total 1-el (contr=T) or nuc-nuc repul (contr=F) gradient
c dens,edens: unchanged
c 
c all triangles follow the (1,1),(2,1),(2,2),(3,1),... order
c
      subroutine GAMDE1 (Sder,Hder,DE,nL2, dens,edens, contr)
      implicit double precision (a-h,o-z)
#include "gamess.par"
      dimension Sder(nL2,3,*),Hder(nL2,3,*)
      dimension dens(*),edens(*)
      dimension DE(3,*)
      logical contr
      dimension tab(500,500)

      common /xFMCOM / x(1)
      common /xINFOA / nat,ich,mul,num,nqmt,ne,na,nb,
     >                zan(mxatm),c(3,mxatm)
      logical goparr

      l1 = num
      l2 = (num*num+num)/2
      l3 = num*num

      call xVCLR (DE,1,3*nat)
      call xVCLR (Sder, 1, 3*l2*nat)
      call xVCLR (Hder, 1, 3*l2*nat)

      call yVALFM(loadfm)
      lderiv =loadfm+1
      lderiv1=lderiv +l2*3*nat
      last   =lderiv1 +l2*3*nat
      need = last - loadfm - 1
      call yGETFM(need)

      call GAMDER_vnn (tab,nat,DE)

      call ySDER (x(lderiv),l2,edens,DE,contr)
      call XCOPY (3*l2*nat, x(lderiv),1, Sder,1)

      call yHELFEY (x(lderiv),l2,dens,DE,contr)
      call XCOPY (3*l2*nat, x(lderiv),1, Hder,1)

      call yTVDER (x(lderiv),x(lderiv1),l2,dens,DE,contr)
      call XAXPY (3*l2*nat, 1.d0,x(lderiv),1, Hder,1)
      call XAXPY (3*l2*nat, 1.d0,x(lderiv1),1, Hder,1)

      call yRETFM(need)

cwc      if(goparr) call ddi_gsumf(1500,de,3*nat)
cwc      call ydawrit(idaf,ioda,de,3*nat,3,0)
      end
c--------------------------------------------------------------------------
      SUBROUTINE yTVDER(DTINT,DVINT,L2,DENAB,DE,contr)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical contr
      dimension DE(3,*)
      LOGICAL OUT,DBG,NORM
      LOGICAL GOPARR,DSKWRK,MASWRK,NXT,SVDSKW
#include "gamess.par"
      PARAMETER (MXCHRM=1)
      DIMENSION DENAB(L2),DTINT(L2,3,NAT),DVINT(L2,3,NAT)
      dimension dijint(225)
      DIMENSION DIJ(225), IJX(35), IJY(35), IJZ(35),
     *          XS(6,7), YS(6,7), ZS(6,7), XT(6,5), YT(6,5), ZT(6,5),
     *          DXS(5,5), DYS(5,5), DZS(5,5),
     *          DXT(5,5), DYT(5,5), DZT(5,5),
     *          XV(6,5,5), YV(6,5,5), ZV(6,5,5),
     *          DXV(5,5,5), DYV(5,5,5), DZV(5,5,5)
      COMMON /xDSTV  / XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,
     *                XJ,YJ,ZJ,NI,NJ,CX,CY,CZ
      COMMON /xIJPAIR/ IA(MXAO)
      COMMON /xINFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM)
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /xOUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /xPAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /xROOT  / XX,U(13),W(13),NROOTS
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00,TWO=2.0D+00)
      PARAMETER (RLN10=2.30258D+00)
      PARAMETER (PI212=1.1283791670955D+00)
      PARAMETER (SQRT3=1.73205080756888D+00)
      PARAMETER (SQRT5=2.23606797749979D+00,SQRT7=2.64575131106459D+00)
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

      IAZ=0
      TOL = RLN10*ITOL
      NORM = NORMF .NE. 1 .OR. NORMP .NE. 1
      CALL xVCLR(DTINT,1,L2*3*NAT)
      CALL xVCLR(DVINT,1,L2*3*NAT)
C
C INITIALIZE PARALLEL
C
      NXT = IBTYP.EQ.1
      IPCOUNT = ME - 1
      NEXT = -1
      MINE = -1

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
          IF((.NOT.NXT) .AND. GOPARR) THEN
             IPCOUNT = IPCOUNT + 1
cwc             IF(MOD(IPCOUNT,NPROC).NE.0) GO TO 1300
          END IF
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
                DEN=ONE
                IF(LOCI+I.EQ.LOCJ+J) DEN=TWO
                dijint(ij)=dum2*den
                NN = IA(MAX0(LOCI+I,LOCJ+J))+MIN0(LOCI+I,LOCJ+J)
                DEN = DENAB(NN)
                DEN = DEN + DEN
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
             NN=IA(MAX0(LOCI+I,LOCJ+J))+MIN0(LOCI+I,LOCJ+J)
             DUMINT=DIJint(IJ)
             DTINT(NN,1,IAT)=DTINT(NN,1,IAT)+DUMX*DUMINT
             DTINT(NN,2,IAT)=DTINT(NN,2,IAT)+DUMY*DUMINT
             DTINT(NN,3,IAT)=DTINT(NN,3,IAT)+DUMZ*DUMINT
             if (contr) then
                DE(1,IAT)=DE(1,IAT)+ DUMX*DIJ(IJ)
                DE(2,IAT)=DE(2,IAT)+ DUMY*DIJ(IJ)
                DE(3,IAT)=DE(3,IAT)+ DUMZ*DIJ(IJ)
             end if
  670     CONTINUE
  680   CONTINUE
C
C     ..... NUCLEAR ATTRACTION
C
  700 CONTINUE
      AAX = AA*AX
      AAY = AA*AY
      AAZ = AA*AZ
      IF(ICHARM.EQ.0) THEN
         MAXIC=NAT
      ELSE
         MAXIC=1
      END IF
C
      DO 790 IC = 1,MAXIC
               ZNUC = -ZAN(IC)
               CX = C(1,IC)
               CY = C(2,IC)
               CZ = C(3,IC)
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
             dum=dij(ij)*aa1*pi212
             DUMINT=DIJint(IJ)*AA1*PI212
                NN=IA(MAX0(LOCI+I,LOCJ+J))+MIN0(LOCI+I,LOCJ+J)
                DVINT(NN,1,IAT)=DVINT(NN,1,IAT)+DUMX*DUMINT
                DVINT(NN,2,IAT)=DVINT(NN,2,IAT)+DUMY*DUMINT
                DVINT(NN,3,IAT)=DVINT(NN,3,IAT)+DUMZ*DUMINT
                if (contr) then
                   DE(1,IAT)=DE(1,IAT)+DUMX*dum
                   DE(2,IAT)=DE(2,IAT)+DUMY*dum
                   DE(3,IAT)=DE(3,IAT)+DUMZ*dum
                end if
  780    CONTINUE
  790 CONTINUE
C
 1000 CONTINUE
 1100 CONTINUE
C
C     ----- END OF PRIMITIVE LOOPS -----
C
 1300 CONTINUE
 1400 CONTINUE
C
C     ----- END OF SHELL LOOPS -----
C
 2000 CONTINUE

cwc        IF(GOPARR) CALL DDI_GSUMF(1503,DTINT,L2*3*NAT)
cwc        IF(GOPARR) CALL DDI_GSUMF(1504,DVINT,L2*3*NAT)
      RETURN
      END
c--------------------------------------------------------------------------
      SUBROUTINE yHELFEY(HFINT,L2,DENAB,DE,contr)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical contr
      dimension DE(3,*)
      LOGICAL IANDJ,OUT,DBG,NORM,DOUBLE,GOPARR,DSKWRK,MASWRK,NXT,SVDSKW
#include "gamess.par"
      DIMENSION DENAB(L2),HFINT(L2,3,NAT)
      dimension dijint(225)
      DIMENSION IJX(35),IJY(35),IJZ(35),DIJ(225),
     *          XIN(5,5,5,2),YIN(5,5,5,2),ZIN(5,5,5,2)
      COMMON /xDSTV  / XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,
     *                XJ,YJ,ZJ,NI,NJ,CX,CY,CZ
      COMMON /xIJPAIR/ IA(MXAO)
      COMMON /xINFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM)
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /xOUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /xPAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /xROOT  / XX,U(13),W(13),NROOTS
C
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00)
      PARAMETER (RLN10=2.30258D+00)
      PARAMETER (PI212=1.1283791670955D+00)
      PARAMETER (SQRT3=1.73205080756888D+00)
      PARAMETER (SQRT5=2.23606797749979D+00)
      PARAMETER (SQRT7=2.64575131106459D+00)
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
      TOL = RLN10*ITOL
      NORM = NORMF.NE.1 .OR. NORMP.NE.1
      CALL xVCLR(HFINT,1,L2*3*NAT)
C
C     INITIALIZATION FOR PARALLEL
C
      NXT = IBTYP.EQ.1
      IPCOUNT = ME - 1
      NEXT = -1
      MINE = -1
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
        IF((.NOT.NXT) .AND. GOPARR) THEN
           IPCOUNT = IPCOUNT + 1
cwc           IF(MOD(IPCOUNT,NPROC).NE.0) GO TO 580
        END IF
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
                dijint(NN)=DUM2*PI212*AA1

                NDUM = IA(LOCI+I)+(LOCJ+J)
                DEN = DENAB(NDUM)
                IF(.NOT.IANDJ.OR.I.NE.J) DEN=DEN+DEN
                DIJ(NN)=DUM2*DEN*PI212*AA1
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

             dumi=dijint(ij)
             DUM=DIJ(IJ)
               NN=IA(LOCI+I)+(LOCJ+J)
               HFINT(NN,1,IC)=HFINT(NN,1,IC)+DUMX*DUMi
               HFINT(NN,2,IC)=HFINT(NN,2,IC)+DUMY*DUMi
               HFINT(NN,3,IC)=HFINT(NN,3,IC)+DUMZ*DUMi

               if (contr) then
               DE(1,IC) = DE(1,IC)+DUM*DUMX
               DE(2,IC) = DE(2,IC)+DUM*DUMY
               DE(3,IC) = DE(3,IC)+DUM*DUMZ
               end if

  460       CONTINUE
  480     CONTINUE
  500 CONTINUE
  520 CONTINUE
C
C     ----- END OF *PRIMITIVE* LOOPS -----
C
  580 CONTINUE
  600 CONTINUE
C
C     ----- END OF *SHELL* LOOPS -----
C

cwc        IF(GOPARR) CALL DDI_GSUMF(1502,HFINT,L2*3*NAT)

      RETURN
      END
c--------------------------------------------------------------------------
      SUBROUTINE ySDER(DSINT,L2,EPS,DE,contr)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical contr
      dimension DE(3,*)
      LOGICAL DBG,NORM,OUT
      LOGICAL GOPARR,DSKWRK,MASWRK,NXT,SVDSKW,NOCON
#include "gamess.par"
      DIMENSION EPS(L2),DSINT(L2,3,NAT)
      DIMENSION DIJ(225), DIJINT(225), IJX(35), IJY(35),IJZ(35),
     *           XS(6,5), YS(6,5), ZS(6,5),
     *          DXS(5,5),DYS(5,5),DZS(5,5)
      COMMON /xDSTV  / XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,
     *                XJ,YJ,ZJ,NI,NJ,CX,CY,CZ
      COMMON /xIJPAIR/ IA(MXAO)
      COMMON /xINFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM)
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /xOUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /xPAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      PARAMETER (SQRT3=1.73205080756888D+00)
      PARAMETER (SQRT5=2.23606797749979D+00, SQRT7=2.64575131106459D+00)
      PARAMETER (ONE=1.0D+00)
      PARAMETER (RLN10=2.30258D+00)
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
      TOL = RLN10*ITOL
c      write (*,*) 'tol=',tol
      NORM = NORMF .NE. 1 .OR. NORMP .NE. 1
      CALL xVCLR(DSINT,1,L2*3*NAT)
C
C     INITIALIZE PARALLEL
C
      NXT = IBTYP.EQ.1
      IPCOUNT = ME - 1
      NEXT = -1
      MINE = -1
C
C     ----- I SHELL
C
      DO 780 II = 1,NSHELL
      IF(NXT .AND. GOPARR) THEN
         MINE = MINE + 1
cwc         IF(MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
cwc         IF(NEXT.NE.MINE) GO TO 780
      END IF

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
      IF((.NOT.NXT) .AND. GOPARR) THEN
         IPCOUNT = IPCOUNT + 1
cwc         IF(MOD(IPCOUNT,NPROC).NE.0) GO TO 760
      END IF
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
              DIJINT(IJ)=DUM2
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
            NN=IA(LOCI+I)+(LOCJ+J)
            DUMINT=DIJINT(IJ)
            DSINT(NN,1,IAT)=DSINT(NN,1,IAT)+DUMX*DUMINT
            DSINT(NN,2,IAT)=DSINT(NN,2,IAT)+DUMY*DUMINT
            DSINT(NN,3,IAT)=DSINT(NN,3,IAT)+DUMZ*DUMINT
            DSINT(NN,1,JAT)=DSINT(NN,1,JAT)-DUMX*DUMINT
            DSINT(NN,2,JAT)=DSINT(NN,2,JAT)-DUMY*DUMINT
            DSINT(NN,3,JAT)=DSINT(NN,3,JAT)-DUMZ*DUMINT

            if (contr) then
            DE(1,IAT)=DE(1,IAT)+(DUMX*DIJ(IJ))
            DE(2,IAT)=DE(2,IAT)+(DUMY*DIJ(IJ))
            DE(3,IAT)=DE(3,IAT)+(DUMZ*DIJ(IJ))
            DE(1,JAT)=DE(1,JAT)-(DUMX*DIJ(IJ))
            DE(2,JAT)=DE(2,JAT)-(DUMY*DIJ(IJ))
            DE(3,JAT)=DE(3,JAT)-(DUMZ*DIJ(IJ))
            end if
 600  CONTINUE
C
 620  CONTINUE
 640  CONTINUE
C     ----- END OF PRIMITIVE LOOPS -----
  760 CONTINUE
  780 CONTINUE
C     ----- END OF SHELL LOOPS -----

      END
c--------------------------------------------------------------------------
      subroutine GAMDER_vnn (tab,natm,de)
      implicit double precision (a-h,o-z)
      dimension tab(natm,natm), de(3,*)
#include "gamess.par"
      COMMON /xINFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     >                ZAN(MXATM),C(3,MXATM)
      data zero,one /0.d0,1.d0/

C
C     ----- FORM DISTANCE MATRIX -----
C
      tab(1,1) = ZERO
      DO 130 K = 2,NAT
         tab(K,K) = ZERO
         K1 = K-1
         DO 120 L = 1,K1
            RKL = ZERO
            DO 100 I = 1,3
               RKL = RKL+(C(I,K)-C(I,L))**2
  100       CONTINUE
            tab(K,L) = -ONE/RKL
            tab(L,K) = SQRT(RKL)
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
               PKL = (C(KK,K)-C(KK,L))/tab(L,K)
               DE(KK,K) = DE(KK,K)+PKL*tab(K,L)*ZAK*ZAL
  230       CONTINUE
  240    CONTINUE
C
         NAT1 = NAT-1
         DO 280 K = 1,NAT1
            ZAK = ZAN(K)
            KP1 = K+1
            DO 270 L = KP1,NAT
               ZAL = ZAN(L)
               PKL = (C(KK,K)-C(KK,L))/tab(K,L)
               DE(KK,K) = DE(KK,K)+PKL*tab(L,K)*ZAK*ZAL
  270      CONTINUE
  280    CONTINUE
  290 CONTINUE

cwc      IF(GOPARR) CALL XSCAL(3*NAT,ONE/NPROC,DE,1)
      end
