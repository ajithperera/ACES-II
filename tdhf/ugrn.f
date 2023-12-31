      SUBROUTINE UGRN(U,ABP,ABM,ABPM,ABVEC,AMVEC,ABVEC1,IVRT,IOCC,
     X IVRTS,IOCCS,PORT,UVAL,UVEC,WVAL,WVEC,SVEC,REDVEC
     X ,DEN,FP,F,EVAL,EVEC,ENRG,UPDATE,ASMALL,ASQUARE,ASCALE,ICONV,IA
     X ,HMO,E,MCOMP,XX,IX,NIREP,NSIZ1,NSIZ3,NSIZO,NSIZVO,NINTMX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*1 AXYZ 
      DIMENSION IVRT(NSIZVO),IOCC(NSIZVO),IVRTS(NSIZVO,NIREP)
     X ,IOCCS(NSIZVO,NIREP)
      DIMENSION U(NSIZ1,NSIZO,MCOMP),F(NSIZ1,NSIZ1,MCOMP)
     X ,ABP(1),ABM(1),ABPM(1),ABVEC(NSIZVO,MCOMP),ABVEC1(NSIZVO,MCOMP)
     X ,AMVEC(NSIZVO,MCOMP),E(NSIZO,NSIZO,MCOMP)
     X ,ENRG(3,MCOMP),EVAL(NSIZ1),EVEC(NSIZ1,NSIZ1)
     X ,PORT(1),UVAL(1),UVEC(1),WVAL(1),WVEC(1)
      DIMENSION DEN(1),FP(1),HMO(NSIZ3,3),IA(1)
      DIMENSION SVEC(1),REDVEC(1)
      DIMENSION XX(NINTMX),IX(NINTMX)
C ....... for linear reduced equation solver
      DIMENSION UPDATE(1),ASMALL(1),ASQUARE(1),ASCALE(1),ICONV(1)
      DIMENSION AXYZ(3),L1(9),L2(9)
      COMMON /LRDUCE/ CONVI
      COMMON /LRDUC1/ KMAX
      COMMON/SWPPP/INP,IAMO,IFAMO,IINDO,IORTH               
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/INFOA/NBASIS,NUMSCF,NX,NMO2,NOC,NVT,NVO     
      COMMON/INFSYM/NSYMHF,NSO(8),NOCS(8),NVTS(8),IDPR(8,8),NVOS(8)
     X ,NIJS(8),NIJS2(8),NRDS(8),NIJSR(8)
      COMMON/INFPRS/IPRSYM(12),NPRSYM(8),JPRSYM(8,12)
      COMMON/THRES/ TOLPER,DEGEN,EPSI
      COMMON/THRE1/  NITER,MAXIT
      COMMON/IPRNT/ IOPEV,IOPDA,IOPU,IOPFE,IOPPR,IWRPA
      DATA ZERO,TWO,FOUR/0.D0,2.D0,4.D0/ 
C ........................  (xx,yy,zz,xy,yz,zx,yx,zy,xz) ......
      DATA AXYZ/'x','y','z'/,L1/1,2,3,1,2,3,2,3,1/,
     X                       L2/1,2,3,2,3,1,1,2,3/
      WRITE(6,*) ' now in UGRN , MCOMP = ',MCOMP
C     WRITE(6,*) ' A+B matrix '
C     CALL OUTMXD(AB,NSIZVO,NVO,NVO)         
      DO 29  MCOM=1,MCOMP
C     write(6,*) ' U with MCOM = ',MCOM
C     CALL OUTMXD(U(1,1,MCOM),NSIZ1,NBASIS,NOC)
      MCOMS=IPRSYM(MCOM+3)
      NVOSYM=NVOS(MCOMS)
C     write(6,*) ' mcoms,nvosym ',MCOMS,NVOSYM
      DO 1000 I=1,NVOSYM
      IF(NITER.NE.0) THEN        
      ABVEC1(I,MCOM)=-U(IVRTS(I,MCOMS),IOCCS(I,MCOMS),MCOM)
C     write(6,*) I,IVRTS(I,MCOMS),IOCCS(I,MCOMS),ABVEC(I,MCOM)
      ELSE
      ABVEC(I,MCOM) = ZERO
      END IF
 1000 CONTINUE
   29 CONTINUE
      DO 92 MCOM=1,MCOMP
      MCOMS=IPRSYM(MCOM+3)
      NVOSYM=NVOS(MCOMS)
      NREDUC=NRDS(MCOMS)
C     CALL GRNSTAT(ABP(NIJS2(MCOMS)+1),ABM(NIJS2(MCOMS)+1),ABPM(NIJS2
C    X(MCOMS)+1),UVAL(NIJS(MCOMS)+1),UVEC(NIJS2(MCOMS)+1),WVAL(NIJS
C    X(MCOMS)+1),WVEC(NIJS2(MCOMS)+1),PORT(NIJS2(MCOMS)+1),
C    X ABVEC(1,MCOM),AMVEC(1,MCOM),ABVEC1(1,MCOM),NVOSYM,NREDUC)
      CALL RPA0(SVEC(NIJSR(MCOMS)+1),WVAL(NIJS(MCOMS)+1),REDVEC,
     X ABVEC(1,MCOM),ABVEC1(1,MCOM),NVOSYM,NREDUC)
   92 CONTINUE
      DO 31 MCOM=1,MCOMP
      MCOMS=IPRSYM(MCOM+3)
      DO 15 I=1,NX
      FP(I)=ZERO
   15 DEN(I)=ZERO
      DO 16 I=1,NVO
      IVT=IVRT(I)
      IOC=IOCC(I)
   16 U(IVT,IOC,MCOM)=ZERO
      NVOSYM=NVOS(MCOMS)
      WRITE(6,*) ' Solution for the linear equation '
      DO 20 I=1,NVOSYM
      IVT=IVRTS(I,MCOMS)
      IOC=IOCCS(I,MCOMS)
      U(IVT,IOC,MCOM)=ABVEC(I,MCOM)
      IJD=(IVT-1)*IVT/2+IOC
      DEN(IJD)=TWO*U(IVT,IOC,MCOM)
C     write(6,*) I,IVT,IOC,ABVEC(I,MCOM)
   20 CONTINUE
      DO 18 M=1,3
      ENRGS= ZERO
      DO 19 I=1,NVO
      IJD=(IVRT(I)-1)*IVRT(I)/2+IOCC(I)
   19 ENRGS = ENRGS + U(IVRT(I),IOCC(I),MCOM)*HMO(IJD,M)
      ENRGS = ENRGS*FOUR
C     write(6,*) ' delta in beta ',ENRGS
C  ... factor two for occupation another two for o-v part
C  ... there is no diagonl elments don't worry
      ENRG(M,MCOM)=ENRG(M,MCOM)+ENRGS
   18 CONTINUE
      WRITE(6,*) ' - Beta '
C     WRITE(6,*) ' component = ',MCOM
C     WRITE(6,*) ' - Beta = ',(ENRG(M,MCOM),M=1,3)
      WRITE(6,200) AXYZ(L1(MCOM)),AXYZ(L2(MCOM)),AXYZ(1),
     X  ENRG(1,MCOM),AXYZ(L1(MCOM)),AXYZ(L2(MCOM)),AXYZ(2)
     X ,ENRG(2,MCOM),AXYZ(L1(MCOM)),AXYZ(L2(MCOM)),AXYZ(3),
     X  ENRG(3,MCOM)
  200 FORMAT(1H ,3A1,' = ',F20.7,'  ',3A1,' = ',F20.7,'  ',
     X3A1,' = ',F20.7)
      IF(IOPU.NE.0) THEN
      WRITE(6,*) ' U '
      CALL OUTMXD(U(1,1,MCOM),NSIZ1,NUMSCF,NOC)               
      END IF
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
      DO 25 I=1,NUMSCF
      DO 25 J=1,I
      IJ=IJ+1
      IJL=(J-1)*NUMSCF+I
      IF(IFAMO.EQ.0) IJL=IJ
      F(I,J,MCOM)=F(I,J,MCOM)+FP(IJL)
      F(J,I,MCOM)=F(I,J,MCOM)
   25 CONTINUE
      IF(IOPFE.NE.0) THEN
      WRITE(6,*) ' F matrix ,component = ',MCOM
      CALL OUTMXD(F(1,1,MCOM),NSIZ1,NUMSCF,NUMSCF)      
      END IF         
      DO 30 I=1,NOC
      DO 30 J=1,NOC
      E(I,J,MCOM)= F(I,J,MCOM)+E(I,J,MCOM)
   30 CONTINUE
      IF(IOPFE.GE.2) THEN
      WRITE(6,*) ' E matrix '
      CALL OUTMXD(E(1,1,MCOM),NSIZO,NOC,NOC)            
      END IF   
   31 CONTINUE
      RETURN
      END
