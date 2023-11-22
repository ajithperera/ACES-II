      SUBROUTINE LINES1(IVRT,IOCC,IVRTS,IOCCS,NPERT,IPN,AMAT,
     X BIA,BIANEW,UPDATE,ASMALL,ASQUARE,ASCALE,ICONV,EVAL,EVEC,NSIZ1,
     X NSIZO,NIREP,F,PRS,CONV,N,KMAX,XX,IX,IA,NINTMX)
C
C  THIS PROGRAM SOLVES THE LINEAR CPHF EQUATION
C
C  SUM A I Z(B,J) (A(BJ,AI)/(EA-EI) -  DELTA(BA)*DELTA(IJ))  = -X(A,I)
C
C  USING AN ITERATIVE EXPANSION IN A SUBSPACE SPANNED BY THE
C  VECTORS B(AI) A'**N (BJ,AI). THE METHOD HAS BEEN FIRST USED
C  BY POPLE AND COWORKERS FOR SOLVING THE CPHF-EQUATIONS AND
C 
C
C  NPERT ... NUMBER OF PERTURBATIONS TREATED
C  AMAT .... CONTAINS THE MODIFIED A-MATRIX AS PROVIDED BY MKARHF
C  BIA  .... B(AI) AS INPUT, U(AI) AS OUTPUT
C  BIANEW .. SCRATCH VECTOR OF LENGTH MAX(N,2*KMAX)
C  UPDATE .. CONTAINS ALL THE EXPANSION VECTOR (LENGTH KMAX*N*NPERT)
C  ASMALL .. THE A MATRIX WITHIN THE ITERATIVE SUBSPACE ( LENGTH 
C            KMAX**2)*NPERT
C  ASQUARE . A VECTOR OF LENGTH (KMAX+1)*NPERT  WHICH HOLDS THE NORMS OF THE 
C            EXPANSION VECTORS
C  ASCALE .. A VECTOR OF LENGTH NPERT WHICH HOLDS SOME SCALING FACTORS
C  ICONV ... A VECTOR HOLDING THE NUMBER OF ITERATIONS REQUIRED TO 
C            CONVERGE A SINGLE PERTURBATION (LENGTH NPERT)
C  EVAL .... EIGENVALUES
C  CONV .... CONVERGENCE CRITERION
C  N ....... NUMBER OF ELEMENTS IN B(AI) OR U(AI)
C  KMAX .... MAXIMUM NUMBER OF ITERATION ALLOWED IN THE ITERATIVE EXPANSION
C
CEND
C 
C CODED AUGUST/90  JG
C modified by HS March/91
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     INTEGER DIRPRD,POP,VRT
C     DIMENSION AMAT(N,N),BIA(N,NPERT),BIANEW(N,NPERT),
      DIMENSION AMAT(1),BIA(N,NPERT),BIANEW(N,NPERT),
     &          UPDATE(N,NPERT,KMAX),ASMALL(KMAX,KMAX,NPERT),
     &          ASQUARE(KMAX+1,NPERT),ASCALE(NPERT),
     &          EVAL(1),ICONV(NPERT)
      DIMENSION IVRTS(N,NIREP),IOCCS(N,NIREP),IVRT(N),IOCC(N)
C ........ AO algorithm ............................
      DIMENSION F(NSIZ1,NSIZ1),PRS(1),EVEC(1)
      DIMENSION XX(NINTMX),IX(NINTMX),IA(1)
      COMMON/INFOA/NBASIS,NUMSCF,NX,NMO2,NOC,NVT,NVO     
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SWPPP/INP,IAMO,IFAMO,IINDO,IORTH               
      COMMON/INFSYM/NSYMHF,NSO(8),NOCS(8),NVTS(8),IDPR(8,8),NVOS(8)
     X ,NIJS(8),NIJS2(8),NRDS(8),NIJSR(8)
      COMMON/INFPRS/IPRSYM(12),NPRSYM(8),JPRSYM(8,12)
C     COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
C     COMMON/SYM/POP(8,2),VRT(8,2),NJUNK(6)
C
      DATA AZERO,ONE,SMALL/0.D+0,1.0D+0,1.D-24/
C
      write(6,*) ' here we are in LINEQ1 '
      CUTOFF=CONV*CONV
      TOL=AZERO
      write(6,*) ' N,NPERT,KMAX ',N,NPERT,KMAX
C
C     CLEAR SMALL A MATRIX
C
      CALL IZERO(ICONV,NPERT)
      CALL ZERO(ASMALL,KMAX*KMAX*NPERT)
C
C     NORMALIZE INITIAL VECTOR BIA AND COPY INTIAL VECTOR TO UPDATE
C
      DO 10 IPERT=1,NPERT
      MCOMS=IPRSYM(IPERT+IPN)
      NVOSYM=NVOS(MCOMS)
CCC     ANORM=SDOT(N,BIA(1,IPERT),1,BIA(1,IPERT),1)
      ANORM=SDOT(NVOSYM,BIA(1,IPERT),1,BIA(1,IPERT),1)
C
C     DETERMINE SCALE FACTOR
C
      ASCALE(IPERT)=SQRT(ANORM)
C
C  IF INITIAL VECTOR EQUAL ZERO, WRITE WARNING
C
CCC     IF(ASCALE(IPERT)/N.LT.CUTOFF) THEN
      IF(ASCALE(IPERT)/NVOSYM.LT.CUTOFF) THEN
       WRITE(6,*) ' Warning from Lineq1 : The initial vector is zero.'
       ASCALE(IPERT)=ONE
      ENDIF
C
C     SCALE INITIAL VECTOR WITH ASCALE
C
      SCALE=ONE/ASCALE(IPERT)
      write(6,*) ' SCALE = ',SCALE,' for ',IPERT
CCC     CALL SSCAL(N,SCALE,BIA(1,IPERT),1)
      CALL SSCAL(NVOSYM,SCALE,BIA(1,IPERT),1)
C
      ASQUARE(1,IPERT)=ONE
   10 CONTINUE
C
C     write(6,*) ' IINTFP ',IINTFP
C     write(6,*) ' IOCC,IVRT,BIA(I,1),BIA(I,2),BIA(I,3) '
C     DO 93 I=1,N
C  93 WRITE(6,*) IOCC(I),IVRT(I),BIA(I,1),BIA(I,2),BIA(I,3)
c YAU : old
c     CALL ICOPY(IINTFP*N*NPERT,BIA(1,1),1,UPDATE(1,1,1),1)
c YAU : new
      CALL DCOPY(N*NPERT,BIA(1,1),1,UPDATE(1,1,1),1)
c YAU : end
C
C    LOOP OVER K UNTIL THE NEW VECTOR IS LINEAR DEPENDENT OR
C    KMAX IS REACHED
C
      DO 1000 K=1,KMAX
C
C     write(6,*) ' iteration = ',K
      IF(IFAMO.EQ.0) THEN
      DO 90 IP=1,NPERT
      MCOMS=IPRSYM(IP+IPN)
      NVOSYM=NVOS(MCOMS)
      CALL XGEMM('T','N',NVOS(MCOMS),1,NVOS(MCOMS),ONE,
     X AMAT(NIJS2(MCOMS)+1),NVOS(MCOMS),BIA(1,IP),NVOS(MCOMS),AZERO,
     X BIANEW(1,IP),NVOS(MCOMS))
      DO 90 I=1,NVOSYM
      IO=IOCCS(I,MCOMS)
      IV=IVRTS(I,MCOMS)
C     IISYM=IDPR(ISYMO(IV),ISYMO(IO))
C     IF(MCOMS.EQ.IISYM) THEN
C     II = IVOS(IV,IO,IISYM)
      BIANEW(I,IP)=BIANEW(I,IP)/(EVAL(IO)-EVAL(IV))
C     END IF
   90 CONTINUE
      ELSE
C   ....... Mada taishousei wa haitteinai .....
      DO 92 IP=1,NPERT
      CALL DENSP(NBASIS,N,IVRT,IOCC,EVEC,BIA(1,IP),PRS)
      IF(IINDO.EQ.0) THEN
      CALL FNOK(PRS,EVEC,F,NBASIS,NUMSCF,XX,IX,NINTMX,IA)
      ELSE
      CALL FINDO(NBASIS,N,IVRT,IOCC,EVEC,PRS,F)
      END IF
      DO 91 I=1,N
      IO=IOCC(I)
      IV=IVRT(I)
      BIANEW(I,IP)=F(IV,IO)/(EVAL(IO)-EVAL(IV))
   91 CONTINUE
   92 CONTINUE
      END IF
C     CALL FORMU(IRREP,NPERT,N,BIANEW,EVAL,POP(1,1),VRT(1,1),NOCC)
C     IF(K.LT.3) THEN
C     WRITE(6,*) ' BIA ',K
C     DO 94 I=1,N
C  94 WRITE(6,*) IOCC(I),IVRT(I),BIA(I,1),BIA(I,2),BIA(I,3)
C     WRITE(6,*) 'BIANEW '
C     DO 95 I=1,N
C  95 WRITE(6,*)IOCC(I),IVRT(I),BIANEW(I,1),BIANEW(I,2),BIANEW(I,3)
C     write(6,*) ' EVAL = ',EVAL(III),III=1,NUMSCF)
C     END IF
C
C    FORM ELEMENTS OF MATRIX ASMALL
C
      DO 20 IPERT=1,NPERT
      MCOMS=IPRSYM(IPERT+IPN)
      NVOSYM=NVOS(MCOMS)
      DO 100 L=1,K
CCC     ASMALL(L,K,IPERT)=SDOT(N,UPDATE(1,IPERT,L),1,BIANEW(1,IPERT),1)
      ASMALL(L,K,IPERT)
     X    =SDOT(NVOSYM,UPDATE(1,IPERT,L),1,BIANEW(1,IPERT),1)
100   CONTINUE
      IF(K.GT.1) ASMALL(K,K-1,IPERT)=ASQUARE(K,IPERT)
C
C    ORTHOGONALIZE XIANEW TO PREVIOUS VECTORS
C
      DO 200 L=1,K
      SCALE=-ASMALL(L,K,IPERT)/ASQUARE(L,IPERT)
CCC     CALL SAXPY(N,SCALE,UPDATE(1,IPERT,L),1,BIANEW(1,IPERT),1)
      CALL SAXPY(NVOSYM,SCALE,UPDATE(1,IPERT,L),1,BIANEW(1,IPERT),1)
200   CONTINUE
C
CCC     ASQUARE(K+1,IPERT)=SDOT(N,BIANEW(1,IPERT),1,BIANEW(1,IPERT),1)
      ASQUARE(K+1,IPERT)
     X    =SDOT(NVOSYM,BIANEW(1,IPERT),1,BIANEW(1,IPERT),1)
20    CONTINUE
      write(6,*) K+1,(ASQUARE(K+1,IPERT),IPERT=1,NPERT)
c YAU : old
c     CALL ICOPY(IINTFP*N*NPERT,BIANEW,1,BIA,1)
c     CALL ICOPY(IINTFP*N*NPERT,BIANEW,1,UPDATE(1,1,K+1),1)
c YAU : new
      CALL DCOPY(N*NPERT,BIANEW,1,BIA,1)
      CALL DCOPY(N*NPERT,BIANEW,1,UPDATE(1,1,K+1),1)
c YAU : end
C
C  CHECK FOR CONVERGENCE OF THE CREATED ITERATIVE SUBSPACE
C
      SUM=AZERO
      DO 30 IPERT=1,NPERT
      MCOMS=IPRSYM(IPERT+IPN)
      NVOSYM=NVOS(MCOMS)
CCC     TEST=ASQUARE(K+1,IPERT)/N
      TEST=ASQUARE(K+1,IPERT)/NVOSYM
      IF(TEST.LE.SMALL) THEN
       ASQUARE(K+1,IPERT)=ONE
       IF(ICONV(IPERT).EQ.0)ICONV(IPERT)=K
      ENDIF
      SUM=SUM+TEST
30    CONTINUE
      IF(SUM.LT.CUTOFF) THEN
C
C      ITERATIVE EXPANSION HAS CONVERGED, EXIT THE LOOP
C
       GO TO 300
      ENDIF
1000  CONTINUE
C
C  IF WE REACH THIS POINT, THE ITERATIVE EXPANSION HAS NOT CONVERGED
C
      WRITE(6,3000) KMAX
3000  FORMAT(' @-LINEQ1-I, CPHF did not converge after ',I4,
     X ' cycles,',' abort !')
      WRITE(6,3001) TEST
3001  FORMAT(' @-LINEQ1-I, Convergence is : ',E20.10)
CC     CALL ERREX
C
C   NOW INVERT THE ASMALL MATRIX
C
  300 CONTINUE
C
      WRITE(6,3002) K
3002  FORMAT(' @-LINEQ1-I, CPHF converged after ',I3,' iterations.')
C
C  SUBTRACT THE DIAGONAL PART 
C
C  NOTE THE ASMALL MATRIX CONSISTS OF SEVERAL PARTS
C
C  FIRST THE UPPER TRIANGLE WITH L LT K
C 
C    THAT'S SIMPLY    Y(L)  A Y(K)
C
C  THE DIAGONAL ELEMENTS
C
C   THAT'S   Y(K) A Y(K) - Y(K) Y(K)
C
C  THE LOWER TRIANGLE WITH L GT K+1
C
C  THESE ARE ZERO SINCE Y(L) SPANS A DIMENSION NOT COVERED BY A Y(K)
C
C  THE LOWER TRIANGLE WITH L EQ K+1
C
C  THESE ARE  Y(K+1) A(K) Y(K) =  Y(K+1) Y(K+1) = ASQUARE(K+1)
C
      CALL ZERO(BIA,N*NPERT)
C
      DO 40 IPERT=1,NPERT
      MCOMS=IPRSYM(IPERT+IPN)
      NVOSYM=NVOS(MCOMS)
C
      K1=K
      IF(ICONV(IPERT).NE.0) K1=ICONV(IPERT)
C
      DO 310 I=1,K1
       ASMALL(I,I,IPERT)=ASMALL(I,I,IPERT)-ASQUARE(I,IPERT)
310   CONTINUE
C
C  INVERT A (XIANEW IS HERE USED AS SCRATCH AND NO LONGER USED
C  FOR HOLDING A*Y
C
      CALL MINV(ASMALL(1,1,IPERT),K1,KMAX,BIANEW,DET,TOL,0,1)
C
      DO 330 I=1,K1
       ASQUARE(I,IPERT)=-ASCALE(IPERT)*ASMALL(I,1,IPERT)
330   CONTINUE
C
C     NOW FORM SOLUTION IN XIA
C
      DO 340 I=1,K1
CCC      CALL SAXPY(N,ASQUARE(I,IPERT),UPDATE(1,IPERT,I),1,
       CALL SAXPY(NVOSYM,ASQUARE(I,IPERT),UPDATE(1,IPERT,I),1,
     & BIA(1,IPERT),1)
340   CONTINUE
40    CONTINUE 
C
C   ALL DONE, RETURN
C
      RETURN
      END
