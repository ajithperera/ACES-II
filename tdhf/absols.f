C                                                                       A  04710
      SUBROUTINE ABSOLS(U,AB,ABVEC,ABVEC1,UVEC,SCR,IVRT,IOCC,
     X IVRTS,IOCCS,DEN,FP,F,E,IA,HMO,EVAL,EVEC,ENRG,UPDATE,
     X ASMALL,ASQUARE,
     X ASCALE,ICONV,MCOMP,NFIRST,XX,IX,NIREP,NSIZ1,NSIZ3,NSIZO,
     X NSIZVO,NINTMX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*1 AXYZ 
      DIMENSION IVRT(NSIZVO),IOCC(NSIZVO)
     X ,IVRTS(NSIZVO,NIREP),IOCCS(NSIZVO,NIREP)
      DIMENSION U(NSIZ1,NSIZO,3),AB(1),ABVEC(NSIZVO,MCOMP),
     X ABVEC1(NSIZVO,MCOMP)
     X ,UVEC(1),SCR(1),EVAL(NSIZ1),EVEC(NSIZ1,NSIZ1),ENRG(3,6)
      DIMENSION F(NSIZ1,NSIZ1,3),E(NSIZO,NSIZO,3),HMO(NSIZ3,3)
      DIMENSION XX(NINTMX),IX(NINTMX),FP(1),IA(1),DEN(1)
C ....... for linear reduced equation solver
      DIMENSION UPDATE(1),ASMALL(1),ASQUARE(1),ASCALE(1),ICONV(1)
      DIMENSION AXYZ(3)
      COMMON /LRDUCE/ CONVI
      COMMON /LRDUC1/ KMAX
      COMMON/SWPPP/INP,IAMO,IFAMO,IINDO,IORTH                
C .........................................
      COMMON/INFOA/NBASIS,NUMSCF,NX,NMO2,NOC,NVT,NVO     
      COMMON/INFSYM/NSYMHF,NSO(8),NOCS(8),NVTS(8),IDPR(8,8),NVOS(8)
     X ,NIJS(8),NIJS2(8),NRDS(8),NIJSR(8)
      COMMON/INFPRS/IPRSYM(12),NPRSYM(8),JPRSYM(8,12)
      COMMON/THRES/ TOLPER,DEGEN,EPSI
      COMMON/THRE1/ NITER,MAXIT
      COMMON/IPRNT/ IOPEV,IOPDA,IOPU,IOPFE,IOPPR,IWRPA
      COMMON/ILINEA/ISALPH,IDALPH
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW 
      DATA ZERO,ONE,TWO,THREE,FOUR/0.D0,1.D0,2.D0,3.D0,4.D0/ 
      DATA AXYZ/'x','y','z'/
      SQRT2=SQRT(TWO)
C ... one root required 
C ..... scratch should have length of 2*NVO at least
      WRITE(6,*) ' now in ABSOLS '
      DO 30 MCOM=1,MCOMP
C     write(6,*) 'MCOM = ',MCOM
      DO 21 I=1,NUMSCF              
      DO 21 J=1,NOC
      IJ=(I-1)*I/2+J 
      IF(I.GT.NOC) THEN 
      U(I,J,MCOM)=HMO(IJ,MCOM)
      ELSE                     
      U(I,J,MCOM)=ZERO
      END IF     
   21 CONTINUE        
C     WRITE(6,*) ' U with MCOM = ',MCOM
C     CALL OUTMXD(U(1,1,MCOM),NSIZ1,NBASIS,NOC)       
C     write(6,*) 'NSIZVO ',NSIZVO
C     WRITE(6,*) ' A+B matrix in ABSOLS with IPRSYM = ',IPRSYM(MCOM)
C     CALL OUTMXD(AB(NIJS2(IPRSYM(MCOM))+1),NVOS(IPRSYM(MCOM)),
C    X NVOS(IPRSYM(MCOM)),NVOS(IPRSYM(MCOM)))               
C     write(6,*) 'i,ivrt,iocc,abvec'
      IPOSYM=IPRSYM(MCOM)
      NVOSYM=NVOS(IPOSYM)
C     write(6,*) ' IPOSYM,NVOSYM ',IPOSYM,NVOSYM
      DO 11 I=1,NVOSYM
C     ABVEC(I,MCOM)=0.D0
C     ABVEC1(I,MCOM)=0.D0
      IF(NITER.NE.0) THEN  
C     IISYM=IDPR(ISYMO(IVRT(I)),ISYMO(IOCC(I)))      
C     IF(IPRSYM(MCOM).EQ.IISYM) THEN
C     II = IVOS(IVRT(I),IOCC(I),IISYM)
      ABVEC(I,MCOM)= -U(IVRTS(I,IPOSYM),IOCCS(I,IPOSYM),MCOM) 
      ABVEC1(I,MCOM)= -U(IVRTS(I,IPOSYM),IOCCS(I,IPOSYM),MCOM) 
C     write(6,*) I,IVRTS(I,IPOSYM),IOCCS(I,IPOSYM),ABVEC(I,MCOM)
C     END IF
      END IF
   11 CONTINUE
   30 CONTINUE
      IF(NITER.NE.0) THEN     
C ... reduced linear equation solver typically set up as NITER=KMAX=20
      IF(NITER.GE.2) THEN
      CONV=1.D-5
      IF(CONVI.NE.0.) CONV=CONVI
C     write(6,*) (EVAL(I),I=1,NUMSCF)
      DO 91 MCOM=1,MCOMP
C     write(6,*) ' MCOM = ',MCOM
C     write(6,*) 'i,ivrt,iocc,abvec'
      IPOSYM=IPRSYM(MCOM)
      NVOSYM=NVOS(IPOSYM)
C     write(6,*) ' IPOSYM,NVOSYM ',IPOSYM,NVOSYM
      DO 91 I=1,NVOSYM
C     IISYM=IDPR(ISYMO(IVRT(I)),ISYMO(IOCC(I)))      
C     IF(IPRSYM(MCOM).EQ.IISYM) THEN
C     II = IVOS(IVRT(I),IOCC(I),IISYM)
      ABVEC(I,MCOM)= ABVEC(I,MCOM)/(EVAL(IVRTS(I,IPOSYM))
     X  -EVAL(IOCCS(I,IPOSYM)))
C     write(6,*) I,IVRTS(I,IPOSYM),IOCCS(I,IPOSYM),ABVEC(I,MCOM)
C     END IF
   91 CONTINUE
      CALL TIMER(0)
      CALL LINES1(IVRT,IOCC,IVRTS,IOCCS,MCOMP,0,AB,ABVEC,
     X ABVEC1,UPDATE,ASMALL,ASQUARE,ASCALE,ICONV,EVAL,EVEC,
C    X NUMSCF,NOC,NIREP,FP,DEN,CONV,NVO,KMAX)
     X NUMSCF,NOC,NIREP,FP,DEN,CONV,NVO,KMAX,XX,IX,IA,NINTMX)
      CALL TIMER(1)
      write(6,*) ' Timing for LINEQ1 for MCOMP = ',MCOMP,TIMENEW
C ....... reduced linear eq. no OWARI ...............
C here is non-iterative part temtatively commented out
      ELSE
C .................................. 
      IF(NFIRST.EQ.0) THEN 
      DO 92 MCOMS=1,NSYMHF
      NVOSYM=NVOS(MCOMS)
C .... EPS,; inverse criteria ...............
C  ....... Create AB matrix ..............
      EPS=1.D-6
      IF(EPSI.NE.0) EPS=EPSI
C ...... NSYMHF matrices will be inversed ........
      DO 10 I=1,NVOSYM
      II2=NIJS2(MCOMS)+(I-1)*NVOSYM+I
   10 AB(II2) = AB(II2)+EVAL(IVRTS(I,MCOMS))-EVAL(IOCCS(I,MCOMS))
C     WRITE(6,*) ' AB matrix to be inversed with MCOMS= ',MCOMS
      CALL MINV(AB(NIJS2(MCOMS)+1),NVOSYM,NVOSYM,SCR
     X ,DET,EPS,0,1)
C     CALL OUTMXD(AB(NIJS2(MCOMS)+1),NVOSYM,NVOSYM,NVOSYM)         
   92 CONTINUE
      END IF
C  ...... Create AB(inverse matrix)= UVEC*(UVAL)**(-1)*UVEC(t)
C  ......  from UVEC=UVEC*(UVAL)**(-1/2) 
C  ...... obtained  by diagonalization
      IF(NFIRST.EQ.2)  THEN
      DO 94 MCOMS=1,NSYMHF
      NVOSYM=NVOS(MCOMS)
      NREDUC=NRDS(MCOMS)
      CALL XGEMM('N','T',NVOSYM,NVOSYM,NREDUC,ONE
     X ,UVEC(NIJS2(MCOMS)+1),NVOSYM,UVEC(NIJS2(MCOMS)+1),NVOSYM
     X ,ZERO,AB(NIJS2(MCOMS)+1),NVOSYM)
   94 CONTINUE
      END IF
C  .........................................................
      DO 93 MCOM=1,MCOMP
      MCOMS=IPRSYM(MCOM)
      NVOSYM=NVOS(MCOMS)
C ........... To get solution of (A+B)X = b    by
C  ........   X = (A+B)**(-1)*b
      CALL MXM(AB(NIJS2(MCOMS)+1),NVOSYM,ABVEC1(1,MCOM),
     X NVOSYM,ABVEC(1,MCOM),1)
C  ...............................................
   93 CONTINUE
      END IF
      ELSE
C ..... here is ordinary iterative scheme .........
      DO 35 MCOM=1,MCOMP
      MCOMS=IPRSYM(MCOM)
C     WRITE(6,*) ' Component and its symmetry ',MCOM,MCOMS
      ITER =0
   40 IF(ITER.GE.MAXIT) GO TO 41
      ITER=ITER+1
      IDIF=0
      CALL MXM(AB(NIJS2(MCOMS)+1),NVOS(MCOMS),ABVEC(1,MCOM),NVOS(MCOMS)
     X ,ABVEC1(1,MCOM),1) 
      NVOSYM=NVOS(MCOMS)
C     write(6,*) ' MCOMS,NVOSYM ',MCOMS,NVOSYM
      DO 36 I=1,NVOSYM
      IVT= IVRTS(I,MCOMS)
      IOC= IOCCS(I,MCOMS)
C     IISYM=IDPR(ISYMO(IVT),ISYMO(IOC))      
C     IF(MCOMS.EQ.IISYM) THEN
C     II = IVOS(IVT,IOC,IISYM)
      VOLD= ABVEC(I,MCOM)
      ABVEC(I,MCOM) =(ABVEC1(I,MCOM)+U(IVT,IOC,MCOM))
     X /(EVAL(IOC)-EVAL(IVT))
      ABVEC1(I,MCOM)= VOLD
      DIF = ABS(ABVEC(I,MCOM)-ABVEC1(I,MCOM))
      IF(DIF.GT.TOLPER) IDIF=1
C     END IF
   36 CONTINUE
      IF(IDIF.EQ.0) GO TO 41
      GO TO 40
   41 CONTINUE
      WRITE(6,*) ' Converged iteration = ',ITER,' within ',TOLPER
   35 CONTINUE
      END IF
      DO 31 MCOM=1,MCOMP    
      MCOMS=IPRSYM(MCOM) 
      DO 15 I=1,NX
      FP(I)=ZERO
   15 DEN(I)=ZERO
      DO 16 I=1,NVO
      IVT=IVRT(I)
      IOC=IOCC(I)
   16 U(IVT,IOC,MCOM)=ZERO
      WRITE(6,*) ' Solution for the linear equation MCOM = ',MCOM
      NVOSYM=NVOS(MCOMS)
C     write(6,*) ' MCOMS,NVOSYM ',MCOMS,NVOSYM
      DO 20 I=1,NVOSYM
      IVT= IVRTS(I,MCOMS)
      IOC= IOCCS(I,MCOMS)
C     U(IVT,IOC,MCOM)= 0.D0
C     IISYM=IDPR(ISYMO(IVT),ISYMO(IOC))      
C     IF(IISYM.EQ.IPRSYM(MCOM)) THEN
C     II = IVOS(IVT,IOC,IISYM)
      U(IVT,IOC,MCOM)=ABVEC(I,MCOM)
C     END IF
      IJD=(IVT-1)*IVT/2+IOC
      DEN(IJD)=TWO*U(IVT,IOC,MCOM)
C     write(6,*) IVT,IOC,U(IVT,IOC,MCOM)
   20 CONTINUE
      IF(IOPU.NE.0) THEN
      WRITE(6,*) ' U1 '
      CALL OUTMXD(U(1,1,MCOM),NSIZ1,NUMSCF,NOC)
      END IF
      DO 18 M=1,3
      ENRG(M,MCOM)= ZERO
      DO 19 I=1,NVO 
      IJD=(IVRT(I)-1)*IVRT(I)/2+IOCC(I)
      ENRG(M,MCOM) = ENRG(M,MCOM) + U(IVRT(I),IOCC(I),MCOM)*HMO(IJD,M)
   19 CONTINUE
      ENRG(M,MCOM) = -ENRG(M,MCOM)*FOUR
   18 CONTINUE
      WRITE(6,*) ' Alpha '
C     WRITE(6,*) ' component = ',MCOM
C     WRITE(6,*) ' Alpha = ',(ENRG(M,MCOM),M=1,3)
      WRITE(6,200) AXYZ(MCOM),AXYZ(1),ENRG(1,MCOM),AXYZ(MCOM)
     X ,AXYZ(2),ENRG(2,MCOM),AXYZ(MCOM),AXYZ(3),ENRG(3,MCOM)
  200 FORMAT(1H ,2A1,' = ',F20.10,'  ',2A1,' = ',F20.10,'  ',
     X 2A1,' = ',F20.10)
      IF(ISALPH.EQ.0) THEN
      IF(IFAMO.EQ.0) THEN
      CALL FMOK(DEN,FP,XX,IX,NINTMX,IA)
      ELSE
      CALL DENSP(NBASIS,NVO,IVRT,IOCC,EVEC,ABVEC(1,MCOM),DEN)
      IF(IINDO.EQ.0) THEN
      CALL FNOK(DEN,EVEC,FP,NBASIS,NSIZ1,XX,IX,NINTMX,IA)
      ELSE
      CALL FINDO(NBASIS,NVO,IVRT,IOCC,EVEC,DEN,FP)
      END IF
      END IF
      IJ=0
      DO 100 I=1,NUMSCF
      DO 100 J=1,I
      F(I,J,MCOM)=ZERO
      IJ=IJ+1
      IJL=(J-1)*NUMSCF+I
      IF(IFAMO.EQ.0) IJL=IJ
      F(I,J,MCOM)= HMO(IJ,MCOM)+FP(IJL)
      F(J,I,MCOM)=F(I,J,MCOM)
  100  CONTINUE
      IF(IOPFE.NE.0) THEN
      WRITE(6,*) '  f-matrix ,component = ',MCOM
      CALL OUTMXD(F(1,1,MCOM),NSIZ1,NUMSCF,NUMSCF)
      END IF
      DO 26 I=1,NOC
      DO 26 J=1,NOC
      E(I,J,MCOM)=F(I,J,MCOM)
   26 CONTINUE
      END IF
   31 CONTINUE
      ABAR=(ENRG(1,1)+ENRG(2,2)+ENRG(3,3))/THREE
      DELA= SQRT((ENRG(1,1)-ENRG(2,2))**2
     X          +(ENRG(2,2)-ENRG(3,3))**2
     X          +(ENRG(3,3)-ENRG(1,1))**2)/SQRT2
      WRITE(6,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      WRITE(6,*) ' Average Alpha = ',ABAR
      WRITE(6,*) ' Delta Alpha = ',DELA
      WRITE(6,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      IF(ISALPH.NE.0) STOP
      RETURN
      END
