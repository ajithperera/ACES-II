      SUBROUTINE CONST2(FI,EI,UI,UJ,E2,DEN,FP,U2,F2,EVAL,EVEC,EHF,
     X HMO,MCOMP,XX,IX,NINTMX,IA,IVRT,IOCC,NSIZ1,NSIZ3,NSIZO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/INFOA/NBASIS,NUMSCF,NX,NMO2,NOC,NVT,NVO
      COMMON/SWPPP/INP,IAMO,IFAMO,IINDO,IORTH
      COMMON/IPRNT/ IOPEV,IOPDA,IOPU,IOPFE,IOPPR,IWRPA
      DIMENSION UI(NSIZ1,NSIZO,3),U2(NSIZ1,NSIZO,MCOMP)
     X ,UJ(NSIZ1,NSIZO,3),E2(NSIZO,NSIZO,MCOMP)
      DIMENSION FI(NSIZ1,NSIZ1,3),EI(NSIZO,NSIZO,3)
     X ,F2(NSIZ1,NSIZ1,MCOMP),EHF(3,MCOMP),HMO(NSIZ3,3),ICOMP(9)
     X ,JCOMP(9),XX(NINTMX),IX(NINTMX),DEN(1),FP(1),EVAL(1)
     X ,EVEC(1),IA(1),IVRT(1),IOCC(1)
      DATA ICOMP/ 1,2,3,1,2,3,2,3,1/, JCOMP/ 1,2,3,2,3,1,1,2,3/
      DATA ZERO/0.D00/,TWO/2.D0/,HALF/5.D-01/,FOUR/4.D0/
      DATA ONE/1.D00/,OCTH/1.D-8/,THREE/3.D0/
      WRITE(6,*) ' Now in CONST2 ,MCOMP = ',MCOMP
C ...... F and E matrices with J index corespond to (-) (inversed )
C ...... in the case of static, it doesn't effect at all(hermitian)
C ...... This routine can handle MCOMP=3 for OR:D(+w,-w)
C ...... dor IRL>3, D(+w,-w) is not symmetrical any more
      DO 50 IRL=1,MCOMP
      IC = ICOMP(IRL)
      JC = JCOMP(IRL)
      WRITE(6,*) ' x,y,z ',IC,JC
      DO 10 I=1,NUMSCF
      DO 10 J=1,NOC
      IF(I.GT.NOC) THEN
      U2(I,J,IRL)=ZERO
      DO 11 K=1,NUMSCF
   11   U2(I,J,IRL)=U2(I,J,IRL)+FI(I,K,IC)*UJ(K,J,JC)
C    X +FJ(I,K,JC)*UI(K,J,IC)
     X +FI(K,I,JC)*UI(K,J,IC)
      DO 13 K=1,NOC
   13  U2(I,J,IRL)=U2(I,J,IRL)
C    X -UI(I,K,IC)*EJ(K,J,JC)-UJ(I,K,JC)*EI(K,J,IC)
     X -UI(I,K,IC)*EI(J,K,JC)-UJ(I,K,JC)*EI(K,J,IC)
      ELSE
      U2(I,J,IRL)=ZERO
      DO 12 K=1,NUMSCF
      U2(I,J,IRL)=U2(I,J,IRL)+UJ(K,I,IC)*UJ(K,J,JC)
     X                       +UI(K,I,JC)*UI(K,J,IC)
   12 CONTINUE
      U2(I,J,IRL)=-HALF*U2(I,J,IRL)
C     U2(J,I,IRL)=U2(I,J,IRL)
      END IF
   10 CONTINUE
      write(6,*) ' IRL = ',IRL
      IF(IOPU.NE.0) THEN
      write(6,*) ' U2 '
      CALL OUTMXD(U2(1,1,IRL),NSIZ1,NUMSCF,NOC)
      END IF
      DO 17 M=1,3
   17  EHF(M,IRL)=ZERO
      IJ=0
      DO 15 I=1,NUMSCF
      DO 15 J=1,I
      IJ=IJ+1
      DEN(IJ) = ZERO
      FP(IJ) = ZERO
      DENLOW= ZERO
      DO 16 K=1,NOC
      DENLOW=DENLOW+ UI(I,K,IC)*UI(J,K,JC)+UJ(I,K,JC)*UJ(J,K,IC)
   16 CONTINUE
      DEN(IJ)=DENLOW*TWO
      IF(I.LE.NOC) DEN(IJ) = DEN(IJ) + FOUR*U2(I,J,IRL)
C     write(6,*) ' den ',I,J,DEN(IJ)
   15 CONTINUE
      DO 18 M=1,3
   18 EHF(M,IRL)= TRACEP(DEN,HMO(1,M),NUMSCF)
      WRITE(6,100)(EHF(I,IRL),I=1,3)
  100 FORMAT(1H0,' CONSTANT ENERGY = ',6F15.7)
      IF(IFAMO.EQ.0) THEN
      CALL FMOK(DEN,FP,XX,IX,NINTMX,IA)
      ELSE
      CALL DENAO(NBASIS,EVEC,DEN,FP)
      IF(IINDO.EQ.0) THEN
      CALL FNOK(DEN,EVEC,FP,NBASIS,NUMSCF,XX,IX,NINTMX,IA)
      ELSE
      CALL FINDO(NBASIS,NVO,IVRT,IOCC,EVEC,DEN,FP)
      END IF
      END IF
      IJ=0
      DO 25 I=1,NUMSCF
      DO 25 J=1,I
      F2(I,J,IRL)=ZERO
      IJ=IJ+1
      IJL=(J-1)*NUMSCF+I
      IF(IFAMO.EQ.0) IJL=IJ
      F2(I,J,IRL)=FP(IJL)
      F2(J,I,IRL)=FP(IJL)
   25 CONTINUE
      IF(IOPFE.NE.0) THEN
      write(6,*) ' F2'
      CALL OUTMXD(F2(1,1,IRL),NSIZ1,NUMSCF,NUMSCF)
      END IF
      DO 30 I=1,NOC
      DO 30 J=1,NOC
      E2(I,J,IRL)=ZERO
      DO 21 K=1,NUMSCF
   21   E2(I,J,IRL)=E2(I,J,IRL)+FI(I,K,IC)*UJ(K,J,JC)
C    X +FJ(I,K,JC)*UI(K,J,IC)
     X +FI(K,I,JC)*UI(K,J,IC)
      DO 22 K=1,NOC
   22   E2(I,J,IRL)=E2(I,J,IRL)
C    X   -UI(I,K,IC)*EJ(K,J,JC)-UJ(I,K,JC)*EI(K,J,IC)
     X   -UI(I,K,IC)*EI(J,K,JC)-UJ(I,K,JC)*EI(K,J,IC)
   30 CONTINUE
      DO 40 I=1,NUMSCF
      DO 40 J=1,NOC
      IF(I.GT.NOC) THEN
      U2(I,J,IRL)=U2(I,J,IRL)+F2(I,J,IRL)
      ELSE
      E2(I,J,IRL)= E2(I,J,IRL) -U2(I,J,IRL)*(EVAL(J)-EVAL(I))
      END IF
   40 CONTINUE
      IF(IOPFE.GE.2) THEN
      write(6,*) ' E2'
      CALL OUTMXD(E2(1,1,IRL),NSIZO,NOC,NOC)
      write(6,*) ' Constant in U'
      CALL OUTMXD(U2(1,1,IRL),NSIZ1,NUMSCF,NOC)
      END IF
   50 CONTINUE
      RETURN
      END
