
#include "flags.h"

       PROGRAM MOLECU
C
C GAUSSIAN INTEGRAL PROGRAM, ORIGINALLY WRITTEN BY J. ALMLOF IN 1971,
C WITH EXTENSIVE MODIFICATIONS BY PETER TAYLOR.
C
C THIS VERSION DERIVES FROM ONE WHICH P. TAYLOR PROVIDED TO QTP IN THE
C SUMMER OF 1990, AND HAS BEEN MODIFIED TO ELIMINATE MACHINE DEPENDENT
C CODE AND ALSO TO ENABLE DYNAMIC ALLOCATION OF CORE.  ALSO, MUCH OF THE
C HARD-CODING IS IN THE PROCESS OF BEING REMOVED.
C
C ECP-INTEGRALS HAVE BEEN ADDED (C.HUBER, UNIVERSITY KARLSRUHE, 1992)
C
C
C     MAXPRIM  - Maximum number of primitives in a shell (30).
C     MAXFNC  - as above (apparently). Number of contracted functions in
C               a shell. Sometimes MAXPRIM is used instead.
C               [MAXPRIM and MAXFNC must keep same value the way things
C                are coded. It is desirable in any case to have them the
C                same so we can do large uncontracted calculations.
C                However, the allocation for READIN and dimensioning in
C                READIN are a bit sloppy, and should be improved.]
C     MXTNPR  - Maximum total number of primitives for all symmetry
C               inequivalent centers (400).
C     MXTNCC  - Maximum total number of contraction coefficients for
C               all symmetry inequivalent centers (800).
C     MXTNSH  - Maximum total number of shells for all symmetry
C               inequivalent centers (100).
C     MXCBF   - Maximum number of Cartesian basis functions for the
C               whole system (NOT the number of contracted functions) (500).
C
C     /INDX/  :
C
C     KMAX          - Number of shells.
C     KHKT(7)       - KHKT(I) = I*(I+1)/2
C     NHKT(MXTNSH)  - L+1 value for each shell (1 for s, 2 for p, etc).
C     NUCO(MXTNSH)  - Number of primitives in each shell.
C     NRCO(MXTNSH)  - Number of contracted functions in each shell.
C
C     /DAT/   :
C
C     ALPHA(MXTNPR)   - Exponents of symmetry inequivalent centers.
C     CONT(MXTNCC)    - Coefficients of symmetry inequivalent centers.
C     CENT(3,MXTNSH ) - Coordinates of symmetry inequivalent shells.
C     CORD(100,3)     - Coordinates of symmetry inequivalent centers.
C     CHARGE(100)     - Nuclear charges of symmetry inequivalent centers.
C     FMULT(8)        - 
C     TLA             - A cutoff 10**-INTGRL_TOL.
C     TLC             - 0.1*TLA ?
C
C     /VMTASK/ :
C     
C     ITASK           -  0  NORMAL CALCULATION (DEFAULT)
C                     -  1  READIN ONLY
C                     -  2  READIN+ONEL+ONELH ONLY
C                     - >3  NORMAL CALCULATION.
C
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C     Common block declarations (or at least some).
C-----------------------------------------------------------------------
C     /REP/
      INTEGER NEWIND,MSTOLD
C     /HARMON/
      LOGICAL HARMT
C     /ONE/
      DOUBLE PRECISION ALP,BET,S1,CSS1,PV,CAB,GAM,DEL,S2,CSS2,QV,CCD
C     /INDX/
      DOUBLE PRECISION PC
      INTEGER DSTRT,NTAP,LU2,NRSS,NUCZ,ITAG,MAXLOP,MAXLOT,KMAX,NMAX,
     &        KHKT,MULT,ISYTYP,ITYPE,AND,OR,EOR,NPARSU,NPAR,MULNUC,
     &        NHKT,MUL,NUCO,NRCO,JSTRT,NSTRT,MST,JRS
C     /DAT/
      DOUBLE PRECISION ALPHA,CONT,CENT,CORD,CHARGE,FMULT,TLA,TLC
C-----------------------------------------------------------------------
C     Parameters
C-----------------------------------------------------------------------
#include "baslims.par"
#include "mxatms.par"
C
      parameter (mxp2=maxprim*maxprim)
      parameter (khm=(nht*(nht+1))/2,nhl=mxp2*khm*khm)
      parameter (kwd=(nht+1)*(nht+2)*(nht+3)/6)
C
C Removed nh4=4*nht-3 from the parameter statement since it is alreay
C set in the ecplims.par common block, 10/07, Ajith Perera
C
      PARAMETER (MXMEMRDIN=4000000)
cch---------------------------------------------------------------------
c     iecp=1 if ecps are used, else iecp=0
cch---------------------------------------------------------------------
      parameter (maxblk = 6)
      parameter (ndilmx = 7)
      parameter (ndient = ndilmx*(ndilmx+1)/2)
      parameter (ndgout = ndient*ndient*maxblk*maxblk)
      common /ecpsts/ iecp
      common /palcali/ nshdim
      common /palcal/ toll
c      dimension gout(ndgout)
cch---------------------------------------------------------------------
      COMMON /TST/ TIM(40),IFREQ(40)
      COMMON /MMMM/ M2(12)
      COMMON /REP/ NEWIND(MXCBF) ,MSTOLD(8)
      COMMON /READINBF/ X(MXMEMRDIN)
      COMMON // ICORE(1)
      COMMON /ISTART/ I0,ICRSIZ
      COMMON /HARMON/HARMT
      COMMON /FLAGS/ IFLAGS(100)
#include "flags2.com"
C
      COMMON /VMTASK/ ITASK
C
      COMMON /ONE/ 
     &   ALP(mxp2),BET(mxp2),S1(mxp2),CSS1(mxp2),PV(mxp2,8,3),CAB(mxp2),
     &   GAM(mxp2),DEL(mxp2),S2(mxp2),CSS2(mxp2),QV(mxp2,8,3),CCD(mxp2)
C
      COMMON /INDX/ PC(512),DSTRT(8,MXCBF),NTAP,LU2,NRSS,NUCZ,ITAG,
     & MAXLOP,MAXLOT,KMAX,NMAX,KHKT(7),MULT(8),ISYTYP(3),ITYPE(7,28),
     & AND(8,8),OR(8,8),EOR(8,8),NPARSU(8),NPAR(8),MULNUC(Mxatms),
     & NHKT(MXTNSH),MUL(MXTNSH),NUCO(MXTNSH),NRCO(MXTNSH),JSTRT(MXTNSH),
     & NSTRT(MXTNSH),MST(MXTNSH),JRS(MXTNSH)
C
      COMMON /DAT/  ALPHA(MXTNPR),CONT(MXTNCC),CENT(3,MXTNSH),
     &              CORD(Mxatms,3),CHARGE(Mxatms),FMULT(8),TLA, TLC
C
      COMMON /SYMIND/ IBFS(MXCBF)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /MEMINF/ IJUNK(2),ITOTMEM,IHWM1,IHWM2
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
C-----------------------------------------------------------------------
C
      IPRINT = IFLAGS(1)
      CALL CRAPSI(X,IUHF,-1)
      INTES=0
      INTE53=0
      INTE32=0
      INTEF=4
C
C     BY DEFAULT, DO THE NORMAL CALCULATION
C
      ITask = 0
C
      CALL ZERO(TIM,40)
      CALL IZERO(IFREQ,40)
      CALL INITP
C
C ALLOCATE CORE FOR READIN ROUTINE
C
      I000=1
      I010=I000+MAXATM
      I020=I010+MXCBF*224
      I030=I020+MXTNCC
      I040=I030+MXTNCC
      I050=I040+MAXFNC*MAXFNC
      I060=I050+MAXFNC
      I070=I060+MAXFNC
      I080=I070+MAXFNC
      I090=I080+MXSHEL*28*30*8
      I100=I090+50
      I110=I100+KWD*3
      I120=I110+NHT*MAXFNC
      I130=I120+NHT*MAXFNC
      I140=I130+8
      I150=I140+100*3
      I160=I150+  8*3
      I170=I160+MXTNCC
      I180=I170+NHT
      I190=I180+8
      I200=I190+18*8
      I210=I200+MXCBF
      I220=I210+MAXFNC*2
      I230=I220+MXCBF*224
      I240=I230+MXTNPR
      I250=I240+100
      I260=I250+100
      I270=I260+3
CJDW Strange allocation for I270 --- IAD array. Only need NHT (up to 7).
      I280=I270+MXCBF
      I290=I280+MXCBF*224
      I300=I290+MXCBF*224
      ITOP=I300+MXCBF*224
cYAU      write(6,*) '  MOLECU-I, ITOP ',itop
      IF(ITOP .GT. MXMEMRDIN)THEN
       WRITE(6,*) '  @MOLECU-F, Insufficient memory for READIN. ',ITOP
       CALL ERREX
      ENDIF
      CALL READIN(X(I000),X(I010),X(I020),X(I030),
     &            X(I040),X(I050),X(I060),X(I070),
     &            X(I080), 
     &            X(I090),X(I100),X(I110),X(I120),
     &            X(I130),X(I140),X(I150),X(I160),
     &            X(I170),X(I180),X(I190),X(I200),
     &            X(I210),X(I220),X(I230),X(I240),
     &            X(I250),X(I260),X(I270),MAXFNC,
CSSS     &            KWD,NH4,X(I280),X(I290),X(I300), 
     &            KWD,X(I280),X(I290),X(I300), 
     &            NSYMOP)
      CALL TIMER(1)
C
C DIVIDE UP CORE 
C
      CALL ACES_MALLOC(ITOTMEM,ICORE,I0)
      IF(I0.EQ.-1)THEN
       WRITE(6,1111)
1111   FORMAT(T3,'@MOLECU-F, Unable to allocate memory.')
       CALL ERREX
      ENDIF

      I000=I0
      J000=I000+IHWM2+MOD(IHWM2,2)
      CALL DCOPY(50+KWD*3,X(I090),1,ICORE(J000),1)
      CALL IZERO(ICORE(I000),J000-I000)
C
C CALCULATE ECPINTS IF IECP = 1
C
c
c allocate memory for ecpvect and gout in paloop
c
#ifdef _OLDECP
      IGOUT=1
      IECPVECT=IGOUT+100000

      IF (IECP.EQ.1) THEN
        IECPSO = IECPVECT + MXCBF*MXCBF                                          
        IF (MXCBF*MXCBF*3+IECPSO>MXMEMRDIN) THEN                                    
CSSS           WRITE(6,*) 'NOT ENOUGH MEMORY FOR SO-ECP CALCULATION!'
CSSS           WRITE(6,*) 'SET IXSIZE TO :', MXCBF*MXCBF*4+IECPVECT                  
        ENDIF                                                                    
        call prepareecp
        Write(6, "(a)") "Entering new ECP module"
        call ecp_int(X(IECPVECT), MXTNPR, 1) 
     
CSSS        Write(6,"(6(1x,F10.7))") (X(iecpvect+j), j=1,50*50)
CSSS        call cvw_ecp_so(X(IECPSO), MXCBF*MXCBF)
CSSS        CALL PALOOP(NSHDIM,TOLL,X(IGOUT),NDGOUT,X(IECPVECT))
        CALL TIMER(1)
103     FORMAT(T3,'@MOLECU-I, ECP-integrals required ',F15.5,
     &            ' seconds.')
      ENDIF
#endif

C
      IG=40000
      J010=J000+50*IINTFP
      J020=J010+KWD*3*IINTFP
c      J020=J020+MOD(J020+1,2)
      J030=J020+IG*IINTFP
      J040=J030+45*IINTFP
      J050=J040+2*IINTFP
      J060=J050+600*IINTFP
      J070=J060+600*IINTFP
      J080=J070+IG*IINTFP
      J090=J080+3*IINTFP
      J100=J090+3*IINTFP
      J110=J100+3*IINTFP
      J120=J110+NHL*IINTFP
      J130=J120+IG*IINTFP
      J140=J130+IG*IINTFP
      J150=J140+2*IINTFP
      J160=J150+2*IINTFP
      J170=J160+3*IINTFP
      J180=J170+14*IINTFP
      J190=J180+3*IINTFP
C
C INTEGRAL PRINT
C
      IF(IPRINT.GE.500)THEN
       WRITE(6,500)
500    FORMAT(T3,'@MOLECU-I, Overlap integrals follow:',/,
     &        T15,'Offset',T37,'Value')
      ENDIF
      IF (ITASK .EQ. 0 .OR. ITASK .GE. 2) THEN
         CALL ONEL(ICORE(J000),ICORE(J010),ICORE(J020),ICORE(J030),
     &             ICORE(J040),ICORE(J050),ICORE(J060),ICORE(J070),
     &             ICORE(J080),ICORE(J090),ICORE(J100),ICORE(J110),
     &             ICORE(J120),ICORE(J130),ICORE(J140),ICORE(J150),
     &             ICORE(J160),ICORE(J170),ICORE(J180),IG,NHL)
      ENDIF
      J030=J020+3*IINTFP
      J040=J030+3*IINTFP
      J050=J040+3*IINTFP
c      J050=J050+MOD(J050+1,2)
      J060=J050+MXP2*8*IINTFP
      J070=J060+MXP2*IINTFP
      J080=J070+NHL*IINTFP
      J090=J080+2*IINTFP
      J100=J090+2*IINTFP
      J110=J100+2*IINTFP
      J120=J110+IG*IINTFP
      J130=J120+100*IINTFP
      J140=J130+IINTFP*IG
      J150=J140+600*IINTFP
      J160=J150+600*IINTFP
      J170=J160+IG*IINTFP
      J180=J170+IG*IINTFP
c     J190=J180+4095*8*3*IINTFP
      CALL CLCNRXLFT(NROWXLFT)
      J190=J180+NROWXLFT*8*3*IINTFP
      J200=J190+3*IINTFP 
      J210=J200+14*IINTFP
      J220=J210+3*IINTFP
      J230=J220+11*IINTFP 
      J240=J230+7*IINTFP
      J250=J240+3*IINTFP
      J260=J250+51*IINTFP
      J270=J260+3*IINTFP
      J280=J270+39*IINTFP
      J290=J280+39*IINTFP
      J300=J290+IG*IINTFP
      J310=J300+IG*IINTFP
      J320=J310+600*IINTFP
C
C INTEGRAL PRINT
C
      IF(IPRINT.GE.500)THEN
       WRITE(6,501)
501    FORMAT(T3,'@MOLECU-I, One electron integrals follow:',/,
     &        T15,'Offset',T24,'Kinetic Energy',T44,'Potential ',
     &        'Energy')
      ENDIF
      If (ITASK .EQ. 0 .OR. ITASK .GE. 2) THEN
         CALL ONELH(ICORE(J000),ICORE(J010),ICORE(J020),ICORE(J030),
     &              ICORE(J040),ICORE(J050),ICORE(J060),ICORE(J070),
     &              ICORE(J080),ICORE(J090),ICORE(J100),ICORE(J110),
     &              ICORE(J120),ICORE(J130),ICORE(J140),ICORE(J150),
     &              ICORE(J160),ICORE(J170),ICORE(J180),ICORE(J190),
     &              ICORE(J200),ICORE(J210),ICORE(J220),ICORE(J230),
     &              ICORE(J240),ICORE(J250),ICORE(J260),ICORE(J270),
     &              ICORE(J280),ICORE(J290),ICORE(J300),ICORE(J310),
     &              X(IECPVECT),IG,NHL,NROWXLFT)
      ENDIF
      CALL TIMER(1)
cYAU      WRITE(6,101)TIMENEW
 101  FORMAT(T3,'@MOLECU-I, One electron integrals required ',F15.5,
     &          ' seconds.')
C
cYAU : old
cc      J000=J000+MOD(J000+1,2)
c      ICORE(I000-1+20)=J000
cC     Next line is bug detected 8/93 while trying to port to KSR.
cC     ICORE(J000-1+50)=ITOTMEM-J000+1
cC     Next line is attempted fix, made 8/19/93.
c      ICORE(J000-1+50)=ITOTMEM-J000+I000-50+1
cYAU : new
      ITMP = J000 - I000
      ICORE(I000-1+20)=ITMP-20

C     WRITE(6,*) 
C
C The above  statement is needed to get vmol to work with  gfortran on 
C Marcel's machine. I do not know whether this is universal. most likely
C an optimization issue.
C
      ICORE(J000-1+50)=ITOTMEM-ITMP-50
cYAU : end

      if (iflags2(h_IFLAGS2_direct).eq.0) then
C
C INTEGRAL PRINT
C
      IF(IPRINT.GE.500)THEN
       WRITE(6,502)
502    FORMAT(T3,'@MOLECU-I, Two electron integrals follow:',/,
     &        T9,'I',T14,'J',T19,'K',T24,'L',T37,'Value')
      ENDIF
      IF (ITASK .EQ. 0 .OR. ITASK .GE. 3) THEN
         CALL TWOEL(ICORE(I000),ICORE(J000), NSYMOP)
      ENDIF
      CALL TIMER(1)
cYAU      WRITE(6,102)TIMENEW
 102  FORMAT(T3,'@MOLECU-I, Two electron integrals required ',F15.5,
     &          ' seconds.')

c     end if (!direct)
      end if

      CLOSE(UNIT=13,STATUS='DELETE')
      call aces_fin
      STOP
      END
