c----------------------------------------------------------------------
c this routine must be called (once) before any calls to GAMDER
      subroutine GAMDER_ini
      implicit double precision (a-h,o-z)
      logical ldum,init
      dimension duma(1),iduma(1)
      data idum,ldum /0,.false./
#include "machsp.com"
      
      open (98,file='MOL')
      open (96,file='GAMDER.LOG')
      call GETINF (98,.true.,96)
      iuhf=0
      na=1
      nb=1
      memory=10 000 000
      call GAMDER_USE_SYMM(.true.)
      call GRDINI (98,.true.,96, iuhf,na,nb,nat,memory,iintfp)
      close (98)
      call ySETFM (memory)
      init=.true.
      call GAMDER (init, idum,idum,idum,idum, ldum,
     > ldum,duma,
     > duma,duma,duma, duma,duma,duma, duma,duma,duma, duma,duma,duma,
     > idum,idum,idum,idum, iduma,iduma,iduma,iduma,
     > ldum,ldum,ldum, ldum,ldum,ldum,ldum)
      write (96,*) 'GRADIENT GAMESS initialization complete'
      end
c----------------------------------------------------------------------
c interface to GAMESS first derivative Cartesian integrals
c ON ENTRY:
c init(logical) must be .false.
c ii,jj,kk,ll: Gamess shell numbers
c all_centers(logical): if .true., calculate contributions to all 4 atoms
c                if .false., let Gamess decide which atoms can be skipped
c ON EXIT:
c done(logical): if .false., the whole quad was skipped (integrals below
c               threshold or all centers equal and all_centers=F);
c               if .false., the integral buffers may contain garbage! 
c integral buffers below must be declared at least 15**4 (for g functions)
c buf(*): regular non-differentiated inegrals as a by-product
c DIx(*): x-derivatives wrt atom 1 (associated with shell ii) etc.
c nofi: No. of orbitals in shell ii etc.
c nbai(*): Gamess orbital numbers for shell ii etc.
c          (must be declared at least 15)
c IIeqJJ, KKeqLL, IJeqKL(logical): symmetry conditions for current quad
c is1(logical): atom 1 was not skipped, integrals available - etc.
c
c to retrieve the integrals from the buffers, use: 
c      ind=0
c      jlast=nofj
c      do 10 i=1,nofi
c      if (IIeqJJ) jlast=i
c      do 10 j=1,jlast
c      klast=nofk
c      if (IJeqKL) klast=i
c      do 10 k=1,klast
c      llast=nofl
c      if (KKeqLL) llast=k
c      if (IJeqKL.and.k.eq.i) llast=j
c      do 10 l=1,llast
c         ind=ind+1
c         ibas=nbai(i)
c         jbas=nbaj(j)
c         kbas=nbak(k)
c         lbas=nbal(l)
c      ..... buf(ind) (etc) contains integral (ibas,jbas|kbas,lbas)
c 10   continue

      subroutine GAMDER (init, ii,jj,kk,ll, all_centers, 
     >   done, buf, DIx,DIy,DIz, DJx,DJy,DJz, DKx,DKy,DKz, DLx,DLy,DLz,
     >   nofi,nofj,nofk,nofl, nbai,nbaj,nbak,nbal, 
     >   IIeqJJ,KKeqLL,IJeqKL, is1,is2,is3,is4)
      implicit double precision (a-h,o-z)
      save
      dimension buf(*)
      dimension DIx(*),DIy(*),DIz(*)
      dimension DJx(*),DJy(*),DJz(*)
      dimension DKx(*),DKy(*),DKz(*)
      dimension DLx(*),DLy(*),DLz(*)
      dimension nbai(*),nbaj(*),nbak(*),nbal(*)
      logical init,all_centers,done,  IIeqJJ,KKeqLL,IJeqKL
      logical is1,is2,is3,is4
      common /gms_deriv_ctrl/ want_integrals, all_cent
      logical want_integrals, all_cent

c--------- declarations from JKDER
cwc1
      COMMON/xSHLEQU/ieqj,keql,ij_kl,IJGTKL,IJLTKL, use_symm
      logical ieqj,keql,ij_kl,IJGTKL,IJLTKL, use_symm
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
      COMMON /xDSHLNO/ LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     *                NIJ,IJD,KLD,IJ,KL
      COMMON /xDSHLT / RTOL,DTOL,VTOL1,VTOL2,VTOLS,OUT,DBG
      COMMON /xFMCOM / X(1)
      COMMON /xIJPAIR/ IA(MXAO)
      COMMON /xINFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM)
      COMMON /xINTFIL/ NINTMX,NHEX,NTUPL,PACK2E,IPOPLE
      COMMON /xINTOPT/ ISCHWZ,IECP,NECP,IEFLD
      COMMON /xIOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /xNSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /xOUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /xPAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /xSHLBAS/ MAXTYP,MAXNUM
      COMMON /xWFNOPT/ SCFTYP,CITYP,DFTYPE,CCTYP,MPLEVL,MPCTYP
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
c------------------------------------------

      if (.not.init) goto 999
c code below is the modified first part (before shell loops) of the 
c routine 'JKDER' from Gamess.  It's used here to initialize variables
c and common blocks

cwc1
      dbg=.false.
      out=.false.
      want_integrals=.true.
cwc0
      SOME = MASWRK  .AND.  NPRINT.NE.-5
C
C     ----- SET STARTING PARAMETERS -----
C
      HONDO = .TRUE.
      POPLE = IPOPLE.EQ.1
cwc1
      pople=.false.
cwc0
C
C      CUTOFF IS THE SCHWARZ SCREENING CUT OFF
C      DABCUT IS THE TWO PARTICLE DENSITY CUT OFF
C
      CUTOFF=0.0d0   !  TENM9
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
      DABCUT=0.0d0   !  TENM11
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
cwc      CALL yDAREAD(IDAF,IODA,DE,3*NAT,3,0)
cwc      IF (GOPARR) CALL XSCAL(3*NAT,ONE/NPROC,DE,1)
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
      NCI = 0
C
C     FIGURE OUT THE MEMORY WE NEED FOR STORING DENSITY MATRIX
C     AND OTHER WAVEFUNCTION INFORMATION. -JKDMEM- ALLOCATES
C     MEMORY FOR DERIVATIVE COMPUTATION AND 2ND ORDER DENSITY
C     AFTER -LENGTH- WORDS.
C
cwc1
      uhftyp=.false.
      hfscf=.true.
cwc0

      LENRHF = L2
      LENUHF = L2+L2
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
           IDFT1=LOADFM+1
           IDFT2=IDFT1+L6
           IDFT3=IDFT2+L6
           IDFT4=IDFT3+L6
           IDFT5=IDFT4+L6
           IDFT6=IDFT5+L6
           IDFT7=IDFT6+L4
         NEED=IADDR-LOADFM
         CALL yGETFM(NEED)
C
C     ----- READ WAVEFUNCTION INFORMATION -----
C
cwc         CALL xDDEBUT(IA,X(LDEN),X(LVEC),L1,L2,L3,NOCORE,NOOPEN,
cwc     *               HFSCF,UHFTYP,ROGVB)
C

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

      return

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 999  continue
c calculation of one shell quadruplet; modified second part of 'JKDER'
c from Gamess
      all_cent=all_centers
      done=.false.

          IF (HONDO) THEN
            IIJJ=IA(MAX0(II,JJ))+MIN0(II,JJ)
            CALL xOEDRD(X(INIJG),NIJ,NIJ0,IIJJ)
            IF(NIJ.EQ.0) return
          END IF
          Q4 = 1.d0
              IF(IPOPLE.EQ.0) POPLE=.FALSE.
              IF(KTYPE(II).GT.2) POPLE=.FALSE.
              IF(KTYPE(JJ).GT.2) POPLE=.FALSE.
              IF(KTYPE(KK).GT.2) POPLE=.FALSE.
              IF(KTYPE(LL).GT.2) POPLE=.FALSE.
          POPLE = .false.
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
CCCCCC              IF (GMAX.LT.CUTOFF) return
C
           ISH=II
           JSH=JJ
           KSH=KK
           LSH=LL
C
C     ----- GET -KL- CHARGE DISTRIBUTION -----
C       ACTUALLY JUST THE POINTERS
C
              KKLL=IA(MAX0(KK,LL))+MIN0(KK,LL)
              CALL xOEDRD(X(INIJG),NKL,NKL0,KKLL)
              IF(NKL.EQ.0) return
C
C     ----- SELECT CENTERS FOR DERIVATIVES -----
C
              CALL xJKDATM(ISH,JSH,KSH,LSH)
              IF(SKIPI.AND.SKIPJ.AND.SKIPK.AND.SKIPL) return
C
C     ----- SET INDICES FOR SHELL BLOCK -----
C
              CALL xJKDSHL(ISH,JSH,KSH,LSH)
              CALL xJKDNDX(X(IIJKLG))
              INEW = ISH
              JNEW = JSH
              KNEW = KSH
              LNEW = LSH
C
C     ----- OBTAIN 2 BODY DENSITY FOR THIS SHELL BLOCK -----
C
cwc         CALL xDABCLU
cwc     *                        (INEW,JNEW,KNEW,LNEW,UHFTYP,X(LDEN),
cwc     *                         X(LVEC),X(IDAB),DABMAX,Q4,POPLE)
cwc1
              dabmax=1.d0
cwc0

C
C
C     ----- EVALUATE DERIVATIVE INTEGRAL, AND ADD TO THE GRADIENT -----
C
         IF(POPLE) THEN
            CALL xJKDG80(DABMAX,INEW,JNEW,KNEW,LNEW,
     *                  JTYPE,IAT,JAT,KAT,LAT)
         ELSE
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

      call extract_integrals (x(igint),x(ifint),x(iijklg),
     >   buf, DIx,DIy,DIz, DJx,DJy,DJz, DKx,DKy,DKz, DLx,DLy,DLz)   

      nofi=MAXI-MINI+1
      nofj=MAXJ-MINJ+1
      nofk=MAXK-MINK+1
      nofl=MAXL-MINL+1
      
      n=0
      do 21 i=MINI,MAXI
         n=n+1
         nbai(n)=LOCI+i
 21   continue
      n=0
      do 22 j=MINJ,MAXJ
         n=n+1
         nbaj(n)=LOCJ+j
 22   continue
      n=0
      do 23 k=MINK,MAXK
         n=n+1
         nbak(n)=LOCK+k
 23   continue
      n=0
      do 24 l=MINL,MAXL
         n=n+1
         nbal(n)=LOCL+l
 24   continue
      
      if (use_symm) then
         IIeqJJ=ieqj
         KKeqLL=keql
      else
         IIeqJJ=.FALSE. ! VFL ieqj
         KKeqLL=.FALSE. ! VFL keql
      endif
      IJeqKL=ij_kl

      is1=.not.skipi
      is2=.not.skipj
      is3=.not.skipk
      is4=.not.skipl

      done=.true.

      RETURN
      END
c-----------------------------------------------------------------------
      subroutine extract_integrals(GIJKL,FIJKL,IJKLG,
     >   buf, DIx,DIy,DIz, DJx,DJy,DJz, DKx,DKy,DKz, DLx,DLy,DLz)

c routine based on Gamess' xJKDOUT, collects the derivative integrals;
      implicit double precision (a-h,o-z)
      dimension buf(*)
      dimension DIx(*),DIy(*),DIz(*)
      dimension DJx(*),DJy(*),DJz(*)
      dimension DKx(*),DKy(*),DKz(*)
      dimension DLx(*),DLy(*),DLz(*)

      DIMENSION GIJKL(*),FIJKL(12,*),IJKLG(4,*)

      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL, use_symm

      COMMON /xDSHLNO/ LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     *                NIJ,IJD,KLD,IJ,KL
      COMMON/xSHLEQU/IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL, use_symm
      COMMON/xSHLNRM/PNRM(35)

      IJKLN=0
      DO 400 I=MINI,MAXI
        PI=PNRM(I)
        JMAX=MAXJ
        IF(use_symm .and.IIEQJJ) JMAX=I
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

              IJKLN=IJKLN+1
              INUM=LOCI+I
              JNUM=LOCJ+J
              KNUM=LOCK+K
              LNUM=LOCL+L

              if (.not. use_symm) then
                 nloci = loci 
                 nlocj = locj 
                 nlock = lock 
                 nlocl = locl 
              endif 

              ind=ijklg(1,ijkln)
              VAL=GIJKL(ind)*PIJKL
              buf(ijkln)=val
              
c              DO 10 N=1,12
c   10         F(N)=FIJKL( N,IJKLG(1,IJKLN))*PIJKL

              DIx(ijkln)=FIJKL(1,ind)*PIJKL
              DIy(ijkln)=FIJKL(2,ind)*PIJKL
              DIz(ijkln)=FIJKL(3,ind)*PIJKL
              DJx(ijkln)=FIJKL(4,ind)*PIJKL
              DJy(ijkln)=FIJKL(5,ind)*PIJKL
              DJz(ijkln)=FIJKL(6,ind)*PIJKL
              DKx(ijkln)=FIJKL(7,ind)*PIJKL
              DKy(ijkln)=FIJKL(8,ind)*PIJKL
              DKz(ijkln)=FIJKL(9,ind)*PIJKL
              DLx(ijkln)=FIJKL(10,ind)*PIJKL
              DLy(ijkln)=FIJKL(11,ind)*PIJKL
              DLz(ijkln)=FIJKL(12,ind)*PIJKL
  100       CONTINUE
  200     CONTINUE
  300   CONTINUE
  400 CONTINUE
      END
c-------------------------------------------------------------------------
      subroutine GRDINI (io,isopen,iout, iuhf,nax,nbx,natx, memory, 
     *                   iintfp)
c initialize commons not taken care of by GETINF/GAMINI
      implicit double precision (a-h,o-z)
c      integer        iflags(100)
c      common /flags/ iflags
#undef _PARALLEL_ACES_COM_
#include "parallel_aces.com"
#include "gamess.par"
      PARAMETER (mxpriL=200,mxorbL=20)
      COMMON /xHERMIT/ H11,H21,H22,H31,H32,H33,H41,H42,H43,H44,
     *                H51,H52,H53,H54,H55,H61,H62,H63,H64,H65,H66,
     *                H71,H72,H73,H74,H75,H76,H77
      COMMON /xWERMIT/ W11,W21,W22,W31,W32,W33,W41,W42,W43,W44,
     *                W51,W52,W53,W54,W55,W61,W62,W63,W64,W65,W66,
     *                W71,W72,W73,W74,W75,W76,W77
      COMMON /xIJPAIR/ IA(MXAO)
      COMMON /xINFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM)
      COMMON /xINTFIL/ NINTMX,NHEX,NTUPL,PACK2E,IPOPLE
      COMMON /xINTOPT/ ISCHWZ,IECP,NECP,IEFLD
      COMMON /xMACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
      COMMON /xPSILVL/ IPSI
      COMMON /xWFNOPT/ SCFTYP,CITYP,DFTYPE,CCTYP,MPLEVL,MPCTYP
      dimension nblks(mxpriL),xdummy(mxorbL)
      character*80 str
      logical isopen,pack2e
      DATA RHF,UHF,ROHF/8HRHF     ,8HUHF     ,8HROHF    /
      DATA RNONE/8HNONE    /
      data zero /0.d0/

      if (.not.isopen) open (io,file='MOL')
      rewind (io)

c here we use MOL only to get the atomic charges; all the other reading
c is just to skip the appropriate number of lines...
      do 1 i=1,5
 1    read (io,*)
      do 10 i=1,nat
      read (io,*) znuc, idummy, nofLs, (nblks(k),k=1,nofLs)
      read (io,*)
      zan(i)=znuc
      do 10 l=1,nofLs
      do 10 nb=1,nblks(l)
         read (io,*) nprims,ncnts
            if (nprims.gt.mxpriL) then
               write (*,'(a,i5)')
     >               '@GRDINI: MAX ALLOWED NO.OF PRIMITIVES PER BLOCK:',
     >               mxpriL
               call aces_exit(1)
            end if
            if (ncnts.gt.mxorbL) then
               write (*,'(a,i5)')
     >               '@GRDINI: MAX ALLOWED NO.OF SHELLS PER BLOCK:',
     >               mxorbL
               call aces_exit(1)
            end if
         do 10 nprim=1,nprims
         read (io,*) dummy, (xdummy(nc),nc=1,ncnts)
 10   continue

      if (.not.isopen) close (io)
 
c default values
      ipople=1
      ischwz=0
      memory=10000000

      open (97,file='ZMAT')
      str=' '
      i=1
      do while (str.ne.'*GAMESS'.and.i.ne.0)
         read (97,'(a)',end=5) str
         goto 6
 5       i=0
 6       continue
      end do
      if (i.ne.0) read (97,'(a)') str
      do while (str.ne.' '.and.str(1:1).ne.'*'.and.i.ne.0)
         i=1
         do while (str(i:i).eq.' '.and.i.lt.(len(str)-8))
            i=i+1
         end do
         if (str(i:i+6).eq.'SCHWRZ=') then
c            if (str(i+7:i+7).eq.'T') ischwz=1
         else if (str(i:i+5).eq.'FDIFF=') then
c            ! compatibility w/ vscf
         else if (str(i:i+5).eq.'POPLE=') then
            if (str(i+6:i+6).eq.'F') ipople=0
         else if (str(i:i+6).eq.'MEMORY=') then
            read(str(i+7:),*) memory
         else if (str.ne.' '.and.str(1:1).ne.'*') then
            if (irank.eq.0) then
               write (*,*) 'error in *GAMESS namelist'
               write (*,*) 'Each line should be of the form:'
               write (*,*) '(SCHWRZ|FDIFF|POPLE)=(T|F)'
            end if
            call aces_exit(1)
         endif
         read (97,'(a)',end=7) str
         goto 8
 7       i=0
 8       continue
      end do
      close (97)

      write (96,*) 'ipople=',ipople
      write (96,*) 'ischwz=',ischwz
      write (96,*) 'memory=',memory*8/1024,' kB'

      NA=NAx
      NB=NBx
      NQMT=NUM
      natx=NAT

      if (iuhf.eq.0) then
         SCFTYP=RHF
      else
         SCFTYP=UHF
      end if
      CITYP=RNONE
      MPLEVL=0
      IPSI=0

      DO 100 I=1,NUM
         IA(I)=(I*I-I)/2
 100  CONTINUE

C     ----- SET -NWDVAR- IN COMMON /MACHIN/ -----
C     THIS IS THE NUMBER OF INTEGERS IN A WORKING PRECISION
C     FLOATING POINT NUMBER, WHICH IS MOST COMMONLY 64 BITS.
C     THUS SYSTEMS WITH 64 BIT INTEGERS USE A VALUE OF 1, WHILE
C          SYSTEMS WITH 32 BIT INTEGERS USE A VALUE OF 2.
C
cMP      NWDVAR = 2
      NWDVAR = iintfp
*CRY  NWDVAR = 1
*CX1  NWDVAR = 1
*DEC  NWDVAR = 1
*ITA  NWDVAR = 1
*NEC  NWDVAR = 1
*T3E  NWDVAR = 1
      MAXFM=0
      MAXSM=0

C     ----- INITIALIZE COMMON /HERMIT/ -----
C     THE VALUES USED BELOW WERE OBTAINED BY RUNNING THE
C     NUMERICAL RECIPE'S ROUTINE IN QUADRUPLE PRECISION.
C     LET THE COMPILER TRUNCATE THEM TO THE MACHINE'S DOUBLE PRECISION.
C
      H11=  ZERO
C
      H21= -0.70710 67811 86547 52440D+00
      H22=  0.70710 67811 86547 52440D+00
C
      H31= -1.22474 48713 91589 04910D+00
      H32=  ZERO
      H33=  1.22474 48713 91589 04910D+00
C
      H41= -1.65068 01238 85784 55588D+00
      H42= -0.52464 76232 75290 31788D+00
      H43=  0.52464 76232 75290 31788D+00
      H44=  1.65068 01238 85784 55588D+00
C
      H51= -2.02018 28704 56085 63293D+00
      H52= -0.95857 24646 13818 50711D+00
      H53=  ZERO
      H54=  0.95857 24646 13818 50711D+00
      H55=  2.02018 28704 56085 63293D+00
C
      H61= -2.35060 49736 74492 22283D+00
      H62= -1.33584 90740 13696 94971D+00
      H63= -0.43607 74119 27616 50868D+00
      H64=  0.43607 74119 27616 50868D+00
      H65=  1.33584 90740 13696 94971D+00
      H66=  2.35060 49736 74492 22283D+00
C
      H71= -2.65196 13568 35233 49245D+00
      H72= -1.67355 16287 67471 44503D+00
      H73= -0.81628 78828 58964 66304D+00
      H74= ZERO
      H75=  0.81628 78828 58964 66304D+00
      H76=  1.67355 16287 67471 44503D+00
      H77=  2.65196 13568 35233 49245D+00
C
C     ----- INITIALIZE COMMON /WERMIT/ -----
C
      W11= 1.77245 38509 05516 02730D+00  ! SQRT(PI)
C
      W21= 8.86226 92545 27580 13649D-01
      W22= 8.86226 92545 27580 13649D-01
C
      W31= 2.95408 97515 09193 37883D-01
      W32= 1.18163 59006 03677 35153D+00
      W33= 2.95408 97515 09193 37883D-01
C
      W41= 8.13128 35447 24517 71430D-02
      W42= 8.04914 09000 55128 36506D-01
      W43= 8.04914 09000 55128 36506D-01
      W44= 8.13128 35447 24517 71430D-02
C
      W51= 1.99532 42059 04591 32077D-02
      W52= 3.93619 32315 22411 59828D-01
      W53= 9.45308 72048 29418 81226D-01
      W54= 3.93619 32315 22411 59828D-01
      W55= 1.99532 42059 04591 32077D-02
C
      W61= 4.53000 99055 08845 64086D-03
      W62= 1.57067 32032 28566 43916D-01
      W63= 7.24629 59522 43925 24092D-01
      W64= 7.24629 59522 43925 24092D-01
      W65= 1.57067 32032 28566 43916D-01
      W66= 4.53000 99055 08845 64086D-03
C
      W71= 9.71781 24509 95191 54149D-04
      W72= 5.45155 82819 12703 05922D-02
      W73= 4.25607 25261 01278 00520D-01
      W74= 8.10264 61755 68073 26765D-01
      W75= 4.25607 25261 01278 00520D-01
      W76= 5.45155 82819 12703 05922D-02
      W77= 9.71781 24509 95191 54149D-04
C
      end
cwc--------------------------------------------------------------------------
cwc routines existing in Gamess, but differing heavily from the originals
cwc--------------------------------------------------------------------------
      subroutine yDAREAD (idaf,ioda,tab,len,nrec,iptyp)
      implicit double precision (a-h,o-z)
#undef _PARALLEL_ACES_COM_
#include "parallel_aces.com"
      character*3 str
      common // icore(1)
      common /imems/ iEa,iEb, ica,icb, iDa,iDb, iL
      common /nbas/ nao,nmo,nat
      dimension grad1 (3000)
      real*8 tab(*)
      integer ioda(*)
c      write (96,*) 'yDAREAD:len,nrec,iptyp',len,nrec,iptyp

      if (nrec.eq.3) then
         call nr2asc (irank,str)
         open (98,file='1ELGRAD.'//str)
         do 10 i=1,3*nat
 10      read (98,*) grad1(i)
         call xcopy (3*nat, grad1,1, tab,1)
         close (98)
      else if (nrec.eq.15) then
         call xcopy (nao*nmo, icore(ica),1, tab,1)
         do i = 1, (nao-nmo)*nao
            tab(nao*nmo+i) = 0.d0
         end do
      else if (nrec.eq.16) then
         l2=(nao*(nao+1))/2
         call xcopy (l2, icore(iDa),1, tab,1)
      else if (nrec.eq.17) then
         call xcopy (nmo, icore(iEa),1, tab,1)
         do i = 1, nao-nmo
            tab(nmo+i) = 0.d0
         end do

      else if (nrec.eq.19) then
         call xcopy (nao*nmo, icore(icb),1, tab,1)
         do i = 1, (nao-nmo)*nao
            tab(nao*nmo+i) = 0.d0
         end do
      else if (nrec.eq.20) then
         l2=(nao*(nao+1))/2
         call xcopy (l2, icore(iDb),1, tab,1)
      else if (nrec.eq.21) then
         call xcopy (nmo, icore(iEb),1, tab,1)
         do i = 1, nao-nmo
            tab(nmo+i) = 0.d0
         end do

      else if (nrec.eq.36) then
         l2=(nao*(nao+1))/2
         call xcopy (l2, icore(iL),1, tab,1)

      else
         write (*,*) 'error in yDAREAD, nrec=',nrec
         call aces_exit(1)
      end if

c      write (96,*) 'DAREAD successful'
      end
c--------------------------------------------------------------------------
      subroutine yDAWRIT (idaf,ioda,tab,len,nrec,iptyp)
      implicit double precision (a-h,o-z)
#undef _PARALLEL_ACES_COM_
#include "parallel_aces.com"
      character*3 str
      real*8 tab(*)
      integer ioda(*)
      common /nbas/ nao,nmo,nat
      dimension grad1 (3000)
      if (nrec.eq.3) then
         call xcopy (3*nat, tab,1, grad1,1)
         call nr2asc (irank,str)
         open (98,file='1ELGRAD.'//str)
         do 10 i=1,3*nat
 10      write (98,'(g23.16)') grad1(i)
         close (98)
      end if
c      call DAWRIT (idaf,ioda,tab,len,nrec,iptyp)
      end
c--------------------------------------------------------------------------
      subroutine ySETFM (mwords)
      COMMON /xFMCOM / icore(1)
      COMMON /MACHSP/ IINTLN, IFLTLN, IINTFP, IALONE, IBITWD
      common /GAMmem/ istart,icurr,iend

      memalo=mwords*IINTFP
      write (96,'(a,i9,a)') 
     > 'allocating remaining    ',mwords,' words for GAMESS routines'
      call ACES_MALLOC (memalo+1,icore,istart)
      if (icore(1).eq.0) then
         write (*,'(a)') '@YSETFM: unable to allocate...'
         call aces_exit(1)
      end if
c convert the starting address from integer to words (Gamess uses words)
      if (mod(istart-1,IINTFP).eq.0) then
         istart=(istart-1)/IINTFP+1
      else
         istart=(istart-1)/IINTFP+2
      end if
      icurr=istart
      iend=istart+mwords-1
      end
c-------------------------------------------------------------
      subroutine yVALFM(loadfm)
      common /GAMmem/ istart,icurr,iend
      loadfm=icurr
      end
c-------------------------------------------------------------
      subroutine yGETFM(need)
      common /GAMmem/ istart,icurr,iend
      icurr=icurr+need
      write (96,'(a,i9,a)') 'GAMESS fills   ',need,' words'
      if (icurr.gt.iend+1) then
         write (*,'(a)') '@YGETFM: allocated memory exceeded...'
         call aces_exit(1)
      end if
      end
c-------------------------------------------------------------
      subroutine yRETFM(need)
      common /GAMmem/ istart,icurr,iend
      icurr=icurr-need
      write (96,'(a,i9,a)') 'GAMESS releases',need,' words'
      if (icurr.lt.istart) icurr=istart
      end
c-------------------------------------------------------------
      subroutine yGOTFM(ngot)
      common /GAMmem/ istart,icurr,iend
      ngot=iend-icurr+1
      end
c-------------------------------------------------------------------------
C*MODULE MTHLIB  *DECK PRTRI
      SUBROUTINE yPRTRI(D,N)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      DIMENSION D(*)
C
      COMMON /xIOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /xOUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /xPAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     ----- PRINT SYMMETRIC MATRIX -D- OF DIMENSION -N- -----
C
      IF (MASWRK) THEN
      MAX = 5
      IF (NPRINT .EQ. 6) MAX = 10
      MM1 = MAX-1
      DO 120 I0=1,N,MAX
         IL = MIN(N,I0+MM1)
cwc         WRITE(IW,9008)
cwc         WRITE(IW,9028) (I,I=I0,IL)
cwc         WRITE(IW,9008)
cwc1
         WRITE(96,9008)
         WRITE(96,9028) (I,I=I0,IL)
         WRITE(96,9008)
cwc0
         IL = -1
         DO 100 I=I0,N
            IL=IL+1
            J0=I0+(I*I-I)/2
            JL=J0+MIN(IL,MM1)
cwc            WRITE(IW,9048) I,(D(J),J=J0,JL)
cwc1
            WRITE(96,9048) I,(D(J),J=J0,JL)
cwc0
  100    CONTINUE
  120 CONTINUE
      END IF
      RETURN
 9008 FORMAT(1X)
 9028 FORMAT(6X,10(4X,I4,4X))
 9048 FORMAT(I5,1X,10F12.7)
      END
c-------------------------------------------------------------------------
C*MODULE HESS    *DECK EGOUT
      SUBROUTINE yEGOUT(EG,NAT)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      DIMENSION EG(3,NAT)
      DIMENSION DNAM(3),EGT(3),ETT(3)
C
#include "gamess.par"
      PARAMETER (MXPT=100, MXFRG=50, MXFGPT=MXPT*MXFRG)
C
      COMMON /xPAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /RUNLAB/ TITLE(10),ANAM(MXATM),BNAM(MXATM),BFLAB(MXAO)
C
      DATA DNAM /4HE'X ,4HE'Y ,4HE'Z /
C
C     ----- PRINT ENERGY GRADIENT VECTOR -----
C
      IF (MASWRK) THEN
         WRITE(96,9000) (DNAM(J),J=1,3)
         DO 220 I=1,NAT
cwc            WRITE(*,9010) I,ANAM(I),BNAM(I),(EG(J,I),J=1,3)
cwc1
            WRITE(96,9010) I,'        ','  ',(EG(J,I),J=1,3)
cwc0
  220    CONTINUE
      END IF
      RETURN
C
 9000 FORMAT(/6X,'ATOM',17X,A4,14X,A4,14X,A4)
 9010 FORMAT(1X,I4,1X,A8,A2,3F18.9)
 9020 FORMAT(6X,'TOTALS',4X,3F18.9)
 9023 FORMAT(6X,'TORQUE',4X,3F18.9)
      END
c-------------------------------------------------------------------------
      subroutine gamder_use_symm(flag)
      logical flag
cwc1
      COMMON/xSHLEQU/ieqj,keql,ij_kl,IJGTKL,IJLTKL, use_symm
      logical ieqj,keql,ij_kl,IJGTKL,IJLTKL, use_symm
cwc0
      use_symm = flag
      return
      end
