      PROGRAM TDHF
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 LAB(52)
      CHARACTER*80 ZSTR
      INTEGER ICORE(1)
      COMMON / / ICORE
      DIMENSION ITITLE(18),LCOMP(21,4)       
      DIMENSION KA(8,3)
      DIMENSION TITLE(24)
      COMMON/INFT/ITITLE,NCOMP,LCOMP,MSIZVO
      COMMON/INFOA/NBASIS,NUMSCF,NX,NMO2,NOC,NVT,NVO     
      COMMON/SWPPP/INP,IAMO,IFAMO,IINDO,IORTH                
C     COMMON/TDHFIN/FREQ,NTDHF,IFIL1,IPROP 
      COMMON/THRES/ TOLPER,DEGEN,EPSI
      COMMON/THRE1/ NITER,MAXIT
      COMMON/TDAMP/  DAMP,IDST,IDEN   
      COMMON/ITOPR/ IDCSHG,IOKE,IDCOR,IIDRI,ITHG
      COMMON/ILINEA/ISALPH,IDALPH
      COMMON/INBETA/ISBETA,ISHG,IEOPE,IOR
      COMMON/IPRNT/ IOPEV,IOPDA,IOPU,IOPFE,IOPPR,IWRPA
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /LRDUCE/ CONVI
      COMMON /LRDUC1/ KMAX
      COMMON/CONST/ EBETA,EGAMMA
      COMMON/INFSYM/NSYMHF,NSO(8),NOCS(8),NVTS(8),IDPR(8,8),NVOS(8)
     X,NIJS(8),NIJS2(8),NRDS(8),NIJSR(8)
      COMMON/MFREQ/NFREQ
      COMMON/DFREQ/VFREQ(100)
      NAMELIST /INPUTP/ IOPEV,IOPDA,IOPU,IOPFE,IOPPR,IWRPA,INP,
     X  FREQ,NITER,EPSI,CONVI,MAXIT,TOLPER,KMAX,LSIZVO,NCOMP,
     X  IDCSHG,IOKE,IDCOR,IIDRI,ITHG,NRDS1,NRDS2,NRDS3,NRDS4,
     X  NRDS5,NRDS6,NRDS7,NRDS8,NFREQ,IAMO,IFAMO,IINDO,IORTH,
     X  DEGEN,DAMP,IDST,IDEN,ISALPH,IDALPH,ISBETA,ISHG,IEOPE,
     X  IOR
      LSIZVO=0
      EBETA=0.8639D0/2.D0
      EGAMMA=0.5037D0/6.D0
C   ************************ Symmetry Table *******************
      DO 40 I=1,8
      DO 40 K=1,3
      KA(I,K)=(I-1)/(2**(K-1))-2*((I-1)/(2**K))
   40 CONTINUE
      DO 42 I=1,8
      DO 42 J=1,8
      IIAN=1
      IIOR=1
      DO 43 K=1,3
      IIAN=IIAN+MIN0(KA(I,K),KA(J,K))*2**(K-1)
   43 IIOR=IIOR+MAX0(KA(I,K),KA(J,K))*2**(K-1)
      IDPR(I,J)=IIOR-IIAN+1
   42 CONTINUE
C   ************************************************************
      NRDS1=-1
      NRDS2=-1
      NRDS3=-1
      NRDS4=-1
      NRDS5=-1
      NRDS6=-1
      NRDS7=-1
      NRDS8=-1
C  *************************************************************
      WRITE(6,*) ' WELCOME TO NON-LINEAR WORLD! ' 
C
C  Read the input stuff from the ZMAT file.  It should be appended to
C  the end of the ZMAT file, after the special  basis set section.
C
      OPEN(UNIT=7,FILE='ZMAT',STATUS='OLD',FORM='FORMATTED',
     &     ACCESS='SEQUENTIAL')
      REWIND(7)
      ZSTR=' '
      DO WHILE (INDEX(ZSTR,'INPUTP').EQ.0)
         READ(7,FMT='(A)',END=7) ZSTR
      END DO
      BACKSPACE(7)
      READ(7,INPUTP,ERR=6)
    5 GOTO 7
    6 WRITE(6,*) '@TDHF: ERROR READING INPUTP NAMELIST.'
      WRITE(6,*) '       Please use the following format in ZMAT:'
      WRITE(6,INPUTP)
      CALL ERREX
    7 CONTINUE
C
C  Now continue on with the stuff.
C
      WRITE(6,*) ' IOPEV,IOPDA,IOPU,IOPFE,IOPPR,IWRPA '
      WRITE(6,*)   IOPEV,IOPDA,IOPU,IOPFE,IOPPR,IWRPA,
     X  IAMO,IFAMO
C   ........ INP  = 0 for CRAPS
      WRITE(6,*) ' INP   = ',INP
      WRITE(6,*) ' IAMO  = ',IAMO
      WRITE(6,*) ' IFAMO = ',IFAMO
      WRITE(6,*) ' IINDO = ',IINDO
      WRITE(6,*) ' IORTH = ',IORTH                
C   ........      = 1 for EVE
      WRITE(6,*) ' NITER = ',NITER
      WRITE(6,*) ' Here NITER = 0,1,>1 stand for 
     X     Iterative, Non-iterative and Reduced Linear Equation '
      WRITE(6,*) 'ISALPH,IDALPH ',ISALPH,IDALPH    
      WRITE(6,*) 'ISBETA,ISHG,IEOPE,IOR ',ISBETA,ISHG,IEOPE,IOR
      WRITE(6,*) 'IDCSHG,IOKE,IDCOR,IIDRI,ITHG'  
      WRITE(6,*)  IDCSHG,IOKE,IDCOR,IIDRI,ITHG   
      IF(NITER.GT.1) THEN
      WRITE(6,*) ' KMAX = ',KMAX,
     X '  is maximun size of Reduced Linear Equation problem '
      WRITE(6,*) ' CONVI = ',CONVI,' if =0. conv=1.D-5'
      END IF
      IF(NITER.EQ.0) THEN
      WRITE(6,*) ' EPSI = ',EPSI,' if =0. eps=1.D-6'     
      WRITE(6,*) ' MAXIT = ',MAXIT
      WRITE(6,*) ' TOLPER = ',TOLPER
      END IF
C 
C  ****** IWRPA = 0 for No output 
C  ******       = 1 for RPA eigenvalues only
C  ******       > 1 for RPA eigenvalues and solutions
C
C  ****** NITER = 0 for Iterative Solution
C  ******       = 1 for Non-iterative by one shot diagonalization
C  ******       > 1 for Reduced Linear Equation Solution
C
C  ****** INP = 0 for CRAPS
C  ******     = 1 for EVE
C
C  ****** IAMO = 0 when one-electron integrals are on AO basis 
C  ******      = 1 when one-electron integrals are on MO basis 
C
C  ****** IORTH = 0 for non orthogonal basis
C  ******       = 1 for     orthogonal basis
C
C  ****** IFAMO = 0 for MO based fock-matrix formation
C  ******       = 1 for AO based fock-matrix formation
C 
C  ****** IINDO = 0 for Non INDO type calculation
C  ******       = 1 for     INDO type calculation
C
      IF(IINDO.NE.0.AND.IORTH.EQ.0) WRITE(6,*) ' !!!!! Come on
     X One electron part should be transformed onto OAO !!!!!'
      IF(NFREQ.NE.1.AND.NITER.NE.1) THEN
      WRITE(6,*) ' NITER.NE.1 I cannot handle multiple frequencies
     X  !!!!!! NITER = ',NITER,' and NFREQ = ',NFREQ
      STOP
      ELSE
      WRITE(6,*) ' Number of calculations with different 
     X frequencies '
      WRITE(6,*) ' NFREQ = ',NFREQ
      WRITE(6,*) ' Applied Frequencies '
      VFREQ(1)=FREQ
      DO 44 I=1,NFREQ
      READ(7,*) VFREQ(I)
      EVFREQ=VFREQ(I)*27.2113957D0
      CMIF=EVFREQ*8065.54093D0
      IF(CMIF.NE.0.) THEN
      CMFR=1.D0/CMIF
      WRITE(6,*) '(',I,') ',VFREQ(I),' a.u.',EVFREQ,' eV',CMIF,' /cm'
     X ,CMFR,' cm'
      ELSE
      WRITE(6,*) ' Applied frequency is zero '
      END IF
   44 CONTINUE
      END IF
C  >>>>>>>>>>>>>>>>>>>>>>>  Reduced Space Informations   >>>>>>>>>>>> 
      NRDS(1)=NRDS1
      NRDS(2)=NRDS2
      NRDS(3)=NRDS3
      NRDS(4)=NRDS4
      NRDS(5)=NRDS5
      NRDS(6)=NRDS6
      NRDS(7)=NRDS7
      NRDS(8)=NRDS8
      WRITE(6,*) ' Size of Reduced Space '
      WRITE(6,*) ' When size is negative, full space will be taken '
      DO 20 I=1,8
   20 WRITE(6,*) 'For ',I,' Symmetry ',NRDS(I)
 220  FORMAT(18A4)                                              
      WRITE(6,*) ' component to be calculated '
      WRITE(6,*) ' NCOMP = ',NCOMP
      IF(NCOMP.LE.21) THEN
      IF(NCOMP.EQ.21) GO TO 131
      DO 130 I=1,NCOMP
      READ(7,*) (LCOMP(I,J),J=1,4)
  130 WRITE(6,*) (LCOMP(I,J),J=1,4)
  131 CONTINUE
      ELSE
      WRITE(6,*) ' Give me a break. I do not want calculate more than
     X 21 components '
      STOP
      END IF
  500 FORMAT(1H1//1X,18a) 
C
      CLOSE(7,STATUS='KEEP')
C
C500  FORMAT(1H1//1H ,'** STEP FOCK **'//1X,18A4//             
C    %        1X,' ========== INPUT DATA =========='/         
C    1       ' NUMBER OF BASIS FUNCTION      = ',I4/         
C    2       ' NUMBER OF SCF DIMENSION       = ',I4/        
C    %       ' NUMBER OF ELECTORONS          = ',I4/       
C    4       ' FILE # OF PROPERTIES          = ',I4/      
C    4       ' FILE # OF INTEGRALS           = ',I4/     
C    *       ' SWITCH FOR P.P.P.             = ',I4///) 
C
      CALL ACES_INIT_RTE

      CALL ACES_JA_INIT
      IF(INP.EQ.0) THEN
C ....from JODA to get informations ......
      CALL GETREC(1,'JOBARC','COMPNIRR',1,NSYMHF)
      CALL GETREC(1,'JOBARC','NOCCORB ',1,NOC)
      CALL GETREC(1,'JOBARC','NVRTORB ',1,NVT)
C     CALL GETREC(1,'JOBARC','OCCUPYA ',8,NOCS)
      CALL GETREC(1,'JOBARC','SYMPOPOA',NSYMHF,NOCS)
      CALL GETREC(1,'JOBARC','SYMPOPVA',NSYMHF,NVTS)
      CALL GETREC(1,'JOBARC','NBASTOT ',1,NBASIS)
      CALL GETREC(1,'JOBARC','NBASCOMP',1,NCOMPA)
      CALL GETREC(1,'JOBARC','NUMDROPA',1,NDROP)
      write(6,*) ' NBASIS,NCOMPA,NDROP ',NBASIS,NCOMPA,NDROP
      NUMSCF=NOC+NVT
      ELSE
C... getting  NUMSCF from EVE
C  ..... For the moment NUMSCF is set to NBASIS
      IFIL1=4
      OPEN(UNIT=IFIL1,FILE='EVE',STATUS='UNKNOWN',FORM='UNFORMATTED')
      READ(IFIL1)LAB,NOC,NMO2,NBASIS
      NUMSCF=NBASIS
      NVT=NUMSCF-NOC
C .... getting symmetry informations and occupation from IIII
      NTAP=8
      OPEN (UNIT=NTAP,FILE='IIII',
     1       STATUS='UNKNOWN',FORM='UNFORMATTED')
      READ(NTAP) TITLE,NSYMHF,(NOCS(I),NVTS(I),I=1,NSYMHF),POTNUC
      CLOSE(NTAP)
      WRITE(6,*) TITLE(1)
      END IF
      WRITE(6,*) ' NSYMHF = ',NSYMHF
      WRITE(6,*) ' occupied   virtual '
      NSUM=0
      DO 50 ISSM=1,NSYMHF
      NSO(ISSM)= NOCS(ISSM)+NVTS(ISSM)
      NSUM = NSUM + NSO(ISSM)
      WRITE(6,*) '(',ISSM,')',NOCS(ISSM),NVTS(ISSM),NSO(ISSM)
   50 CONTINUE
      IF(NSUM.NE.NUMSCF) WRITE(6,*) '!!!! Sum*NSO .NE. NUMSCF 
     X !!!!!!!'
      IF(NSYMHF.GT.1.AND.IFAMO.NE.0) THEN
      WRITE(6,*) ' For the moment, AO algorithm is implemented
     X only for Non-symmetry cases ; so far'
      STOP
      END IF
      WRITE(6,*) 'NOC,NVT,NBASIS,NUMSCF ',NOC,NVT,NBASIS,NUMSCF
C .................................................
      NX=(NUMSCF*(NUMSCF+1))/2
      NMO2=NUMSCF*NUMSCF        
      NVO=NVT*NOC
C  ...... size of core .....
      NSIZ1=NUMSCF
      NSIZO=NOC
      IF(LSIZVO.NE.0) THEN
      WRITE(6,*) 'Size of Matrix to diagonalized when NITER=1'
      WRITE(6,*) 'LSIZVO = ',LSIZVO
      END IF
C   

      WRITE(6,*) ' IINTFP ',IINTFP
C .......  MCOMP = 9 ..........
      MCOMP=9
      NSIZ2=NSIZ1*NSIZ1
      NSIZ3=(NBASIS+1)*NBASIS/2
C  ... EVEC
      IO1=NBASIS*NBASIS*IINTFP+1
C     IO1=NBASIS*NUMSCF*IINTFP+1
C  ... EVAL
      IO2=IO1+NSIZ1*IINTFP
C  ... IA
      IO3=IO2+NSIZ1 + mod(nsiz1,2)
C  ... ISYMO
      IO4=IO3+NSIZ1 + mod(nsiz1,2)
C  ... DAP
      IO5=IO4+NBASIS*NBASIS*IINTFP
C  ... DAM
      IO6=IO5+NBASIS*NBASIS*IINTFP
C  ... FP
      IO7=IO6+NBASIS*NBASIS*IINTFP
C  ... FM
      IO8=IO7+NBASIS*NBASIS*IINTFP
C  ... F1
      IO9=IO8+NBASIS*NBASIS*3*IINTFP
C  ... F2
      IO10=IO9+NSIZ2*6*IINTFP
C  ... FW
      IO11=IO10+NSIZ2*3*IINTFP
C  ... FSW
      IO12=IO11+NSIZ2*9*IINTFP
C  ... FWW
      IO13=IO12+NSIZ2*9*IINTFP
C  ... FS
      IO14=IO13+NSIZ2*3*IINTFP
C  ... FS2
      IO15=IO14+NSIZ2*9*IINTFP
      NSIZO2=NSIZO*NSIZO
C  ... E1
      IO16=IO15+NSIZO2*3*IINTFP
C  ... E2
      IO17=IO16+NSIZO2*6*IINTFP
C  ... EW
      IO18=IO17+NSIZO2*3*IINTFP
C  ... ESW
      IO19=IO18+NSIZO2*9*IINTFP
C  ... EWW
      IO20=IO19+NSIZO2*9*IINTFP
C  ... ES
      IO21=IO20+NSIZO2*3*IINTFP
C  ... ES2 
      IO22=IO21+NSIZO2*6*IINTFP
C  ... HMO 
      IO23=IO22+NBASIS*NBASIS*3*IINTFP
      NSIZ1O=NSIZ1*NSIZO
C  ... UP1
      IO24=IO23+NSIZ1O*3*IINTFP
C  ... UM1
      IO25=IO24+NSIZ1O*3*IINTFP
C  ... UP2
      IO26=IO25+NSIZ1O*6*IINTFP
C  ... UM2
      IO27=IO26+NSIZ1O*6*IINTFP
C  ... UWP
      IO28=IO27+NSIZ1O*3*IINTFP
C  ... UWM
      IO29=IO28+NSIZ1O*3*IINTFP
C  ... US
      IO30=IO29+NSIZ1O*3*IINTFP
C  ... US2
      IO31=IO30+NSIZ1O*9*IINTFP
C  ... USWP
      IO32=IO31+NSIZ1O*9*IINTFP
C  ... USWM
      IO33=IO32+NSIZ1O*9*IINTFP
C  ...  UWWP
      IO34=IO33+NSIZ1O*9*IINTFP
C  ... UWWM 
      IO35=IO34+NSIZ1O*9*IINTFP
      NSIZV=NVT
      NSIZVO=NVO
      DO 9 II=1,NSYMHF
    9 NVOS(II)=0
      DO 10 II=1,NSYMHF
      DO 10 JJ=1,NSYMHF
      IJSYM= IDPR(II,JJ)
      NVOS(IJSYM)=NVOS(IJSYM) + NOCS(II)*NVTS(JJ)
   10 CONTINUE
      WRITE(6,*) '  Size of Full and Reduced Space '
      IOVO2=0
      IOVOR=0
      IORED=0
      DO 11 II=1,NSYMHF
      IOVO2=IOVO2+ NVOS(II)*NVOS(II)
      IF(NRDS(II).LT.0) NRDS(II)=NVOS(II)
      IOVOR=IOVOR+ NVOS(II)*NRDS(II)
      IORED=IORED +NRDS(II)
   11 WRITE(6,*) 'For  Symmetry ',II,'  NVOS,NRDS = ',NVOS(II),
     X NRDS(II) 
C  ...  ABP
      IABP=0
      IF(IFAMO.EQ.0) IABP=IOVO2*IINTFP
      IO36=IO35+IABP + mod(iabp,2)
C  ...  ABM
      IABM=0
      IF(IFAMO.EQ.0) IABM=IOVO2*IINTFP
      IO37=IO36+IABM + mod(iabm,2)
C  ...  IVO
C     IO38=IO37+NSIZ2
      IO38=IO37+NSIZ1*NSIZO + mod(nsiz1*nsizo,2)
C  ... IVRT
      IO39=IO38+NSIZVO + mod(nsizvo,2)
C  ...  IOCC
      IO40=IO39+NSIZVO + mod(nsizvo,2)
C  ... UPDATE
      IO41=IO40+NSIZVO*MCOMP*KMAX*2*IINTFP
C  ... ASMALL
      IO42=IO41+KMAX*KMAX*MCOMP*IINTFP
C  ... ASQUARE
      IO43=IO42+(KMAX+1)*MCOMP*IINTFP
C  ... ASCALE
      IO44=IO43+MCOMP*IINTFP
C  ... ICONV
      IO45=IO44+MCOMP + mod(mcomp,2)
C  ..........
      MSIZVO=NSIZVO
      IF(LSIZVO.NE.0) MSIZVO=LSIZVO
      IF(NITER.EQ.1) WRITE(6,*) 'Size of the Space for Diagonalization
     X = ',MSIZVO
C  ... ABPM
      IABPM=0
      IF(NITER.EQ.1) IABPM=IOVO2*IINTFP
      IO46=IO45+IABPM + mod(iabpm,2)
C  ... UVAL
      IUVAL=0
      IF(NITER.EQ.1) IUVAL=MSIZVO*IINTFP
      IO47=IO46+IUVAL + mod(iuval,2)
C  ... UVEC
      IUVEC=0
      IF(NITER.EQ.1) IUVEC=IOVO2*IINTFP
      IO48=IO47+IUVEC + mod(iuvec,2)
C  ... WVAL
      IWVAL=0
      IF(NITER.EQ.1) IWVAL=MSIZVO*IINTFP
      IO49=IO48+IWVAL + mod(iwval,2)
C  ... WVEC
      IWVEC=0
      IF(NITER.EQ.1) IWVEC=IOVO2*IINTFP
      IO50=IO49+IWVEC + mod(iwvec,2)
C  ... PORT
      IPORT=0
      IF(NITER.EQ.1) IPORT=IOVO2*IINTFP
      IO51=IO50+IPORT + mod(iport,2)
C  ... SCR
      ISCR=0
      NSCR=NSIZ1
      IF(NITER.EQ.1) NSCR=MAX(NSIZ1,MSIZVO)
      NSCR=NSCR+1
      ISCR=NSCR*5*IINTFP
      IO52=IO51+ISCR + mod(iscr,2)
C  ... ABVEC
      IO53=IO52+MCOMP*NSIZVO*IINTFP
C  ... ABVEC1
      IO54=IO53+MCOMP*NSIZVO*IINTFP
C  ... AMVEC
      IO55=IO54+MCOMP*NSIZVO*IINTFP
C  ... AMVEC1
      IO56=IO55+MCOMP*NSIZVO*IINTFP
C  ... IAS
      IO57=IO56+NSYMHF*NSIZ1 + mod(nsymhf*nsiz1,2)
C  ... IVRTS
      IO58=IO57+NSYMHF*NSIZVO + mod(nsymhf*nsizvo,2)
C  ... IOCCS
      IO59=IO58+NSYMHF*NSIZVO + mod(nsymhf*nsizvo,2)
C  ... IVOS
      IO60=IO59+NSYMHF*NSIZ1O + mod(nsymhf*nsiz1o,2)
C  ... ABPI
      IABPI=0
      IF(NITER.EQ.1) IABPI=IOVO2*IINTFP
      IO61=IO60+ IABPI + mod(iabpi,2)
C  ... SVEC 
      ISVEC=0
      IF(NITER.EQ.1) ISVEC=IOVOR*IINTFP
      IO62=IO61+ISVEC + mod(isvec,2)
C  ... TVEC 
      ITVEC=0
      IF(NITER.EQ.1) ITVEC=IOVOR*IINTFP
      IO63=IO62+ITVEC + mod(itvec,2)
C  ... ZVEC 
      IZVEC=0
      IF(NITER.EQ.1) IZVEC=IOVOR*IINTFP
      IO64=IO63+IZVEC + mod(izvec,2)
C  ... YVEC 
      IYVEC=0
      IF(NITER.EQ.1) IYVEC=IOVOR*IINTFP
      IO65=IO64+IYVEC + mod(iyvec,2)
C  ... REDVEC
      IREDV=0
      IF(NITER.EQ.1) IREDV=2*IORED*IINTFP
      IO66=IO65+IREDV + mod(iredv,2)

c Nevin
c      io66=6500000

      WRITE(6,*) ' REQUIRED CORE MEMORY ',IO66
      CALL ACES_MALLOC(IO66,ICORE,I0)
      IF (I0.EQ.-1) THEN
         WRITE(*,*) '@TDHF: Request for ',IO66,' words of memory',
     &              ' failed.'
      CALL ERREX
      END IF
      IO1=NBASIS*NBASIS*IINTFP+I0
C     IO1=NBASIS*NUMSCF*IINTFP+I0
      IO2=IO1+NSIZ1*IINTFP
      IO3=IO2+NSIZ1 + mod(nsiz1,2)
      IO4=IO3+NSIZ1 + mod(nsiz1,2)
      IO5=IO4+NBASIS*NBASIS*IINTFP
      IO6=IO5+NBASIS*NBASIS*IINTFP
      IO7=IO6+NBASIS*NBASIS*IINTFP
      IO8=IO7+NBASIS*NBASIS*IINTFP
      IO9=IO8+NBASIS*NBASIS*3*IINTFP
      IO10=IO9+NSIZ2*6*IINTFP
      IO11=IO10+NSIZ2*3*IINTFP
      IO12=IO11+NSIZ2*9*IINTFP
      IO13=IO12+NSIZ2*9*IINTFP
      IO14=IO13+NSIZ2*3*IINTFP
      IO15=IO14+NSIZ2*9*IINTFP
      IO16=IO15+NSIZO2*3*IINTFP
      IO17=IO16+NSIZO2*6*IINTFP
      IO18=IO17+NSIZO2*3*IINTFP
      IO19=IO18+NSIZO2*9*IINTFP
      IO20=IO19+NSIZO2*9*IINTFP
      IO21=IO20+NSIZO2*3*IINTFP
      IO22=IO21+NSIZO2*6*IINTFP
      IO23=IO22+NBASIS*NBASIS*3*IINTFP
      IO24=IO23+NSIZ1O*3*IINTFP
      IO25=IO24+NSIZ1O*3*IINTFP
      IO26=IO25+NSIZ1O*6*IINTFP
      IO27=IO26+NSIZ1O*6*IINTFP
      IO28=IO27+NSIZ1O*3*IINTFP
      IO29=IO28+NSIZ1O*3*IINTFP
      IO30=IO29+NSIZ1O*3*IINTFP
      IO31=IO30+NSIZ1O*9*IINTFP
      IO32=IO31+NSIZ1O*9*IINTFP
      IO33=IO32+NSIZ1O*9*IINTFP
      IO34=IO33+NSIZ1O*9*IINTFP
      IO35=IO34+NSIZ1O*9*IINTFP
      IO36=IO35+IABP + mod(iabp,2)
      IO37=IO36+IABM + mod(iabm,2)
      IO38=IO37+NSIZ1*NSIZO + mod(nsiz1*nsizo,2)
      IO39=IO38+NSIZVO + mod(nsizvo,2)
      IO40=IO39+NSIZVO + mod(nsizvo,2)
      IO41=IO40+NSIZVO*MCOMP*KMAX*2*IINTFP
      IO42=IO41+KMAX*KMAX*MCOMP*IINTFP
      IO43=IO42+(KMAX+1)*MCOMP*IINTFP
      IO44=IO43+MCOMP*IINTFP
      IO45=IO44+MCOMP + mod(mcomp,2)
      IO46=IO45+IABPM + mod(iabpm,2)
      IO47=IO46+IUVAL + mod(iuval,2)
      IO48=IO47+IUVEC + mod(iuvec,2)
      IO49=IO48+IWVAL + mod(iwval,2)
      IO50=IO49+IWVEC + mod(iwvec,2)
      IO51=IO50+IPORT + mod(iport,2)
      IO52=IO51+ISCR + mod(iscr,2)
      IO53=IO52+MCOMP*NSIZVO*IINTFP
      IO54=IO53+MCOMP*NSIZVO*IINTFP
      IO55=IO54+MCOMP*NSIZVO*IINTFP
      IO56=IO55+MCOMP*NSIZVO*IINTFP
      IO57=IO56+NSYMHF*NSIZ1 + mod(nsymhf*nsiz1,2)
      IO58=IO57+NSYMHF*NSIZVO + mod(nsymhf*nsizvo,2)
      IO59=IO58+NSYMHF*NSIZVO + mod(nsymhf*nsizvo,2)
      IO60=IO59+NSYMHF*NSIZ1O + mod(nsymhf*nsiz1o,2)
      IO61=IO60+ IABPI + mod(iabpi,2)
      IO62=IO61+ISVEC + mod(isvec,2)
      IO63=IO62+ITVEC + mod(itvec,2)
      IO64=IO63+IZVEC + mod(izvec,2)
      IO65=IO64+IYVEC + mod(iyvec,2)
      IO66=IO65+IREDV + mod(iredv,2)
      CALL NONLIN(ICORE(I0),ICORE(IO1),ICORE(IO2),ICORE(IO3),ICORE(IO4)
     X ,ICORE(IO5),ICORE(IO6),ICORE(IO7),ICORE(IO8),ICORE(IO9)
     X ,ICORE(IO10),ICORE(IO11),ICORE(IO12),ICORE(IO13),ICORE(IO14)
     X ,ICORE(IO15),ICORE(IO16),ICORE(IO17),ICORE(IO18),ICORE(IO19)
     X ,ICORE(IO20),ICORE(IO21),ICORE(IO22),ICORE(IO23),ICORE(IO24)
     X ,ICORE(IO25),ICORE(IO26),ICORE(IO27),ICORE(IO28),ICORE(IO29)
     X ,ICORE(IO30),ICORE(IO31),ICORE(IO32),ICORE(IO33),ICORE(IO34)
     X ,ICORE(IO35),ICORE(IO36),ICORE(IO37),ICORE(IO38),ICORE(IO39)
     X ,ICORE(IO40),ICORE(IO41),ICORE(IO42),ICORE(IO43),ICORE(IO44)
     X ,ICORE(IO45),ICORE(IO46),ICORE(IO47),ICORE(IO48),ICORE(IO49)
     X ,ICORE(IO50),ICORE(IO51),ICORE(IO52),ICORE(IO53),ICORE(IO54)
     X ,ICORE(IO55),ICORE(IO56),ICORE(IO57),ICORE(IO58),ICORE(IO59)
     X ,ICORE(IO60),ICORE(IO61),ICORE(IO62),ICORE(IO63),ICORE(IO64)
     X ,ICORE(IO65))
      CALL ACES_JA_FIN
      CALL ACES_EXIT(0)
      END
