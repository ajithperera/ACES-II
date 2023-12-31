      SUBROUTINE VGAMMA2(F,FARG,FSCR,L,N,M)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NNN=7, TMAX=60., NSIZ=600,
     $ xinc=tmax/nsiz,rinc=1./xinc, nacc=2, nndim=4*nnn+nacc)
      DIMENSION F(1),FARG(M),FSCR(M),L(M)
c     DIMENSION F(2),L(2),TAB(NSIZ+1,NNDIM)
      DIMENSION TAB(NSIZ+1,NNDIM)
      SAVE TAB,PI
C     NACC=0,1 OR 2 DETERMINES THE ACCURACY.
C...  ROUTINE FOR INCOMPLETE GAMMA FUNCTION.
C     MEMORY REQUIREMENT IS M*MAX(8,N+4)
C     M IS THE NUMBER OF ARGUMENTS FOR WHICH THE
C     GAMMA FUNCTION IS TO BE EVALUATED.
C     N IS THE ORDER OF THE GAMMA FUNCTION.
C     THE ARGUMENTS ARE PASSED IN F(*,MAX(5,N+1))
C     THE RESULT IS RETURNED IN (F(*,I),  I=1,N+1)
C
      DATA FC /0.5/

#ifdef _TRACE_VMOL
c      write(*,*) '@VMOL: Entered file ',__FILE__
#endif

C     FC IS TO BE OPTIMIZED FOR BEST PERFORMANCE
C     ITS OPTIMUM VALUE IS INSTALLATION DEPENDENT.
C               0 < FC < 1
C     MORE EFFICIENT VECTOR PERFORMANCE AS COMPARED TO SCALAR
C     PERFORMANCE  <===> LARGER FC
C
C
      M3=M+M
      MT=M*(MAX0(5,N+1))
      ME=MT+M
      ML=ME+M
C
C...  COLLECT SMALL T-VALUES IN F(*,1)
C
      CALL COLLEC (K,FARG,FSCR,M,TMAX)
C
      IF (K.GT.0) THEN
C...  THERE ARE SOME T-VALUES <TMAX
C
      CALL GATHER(K,F,FARG,FSCR)
C
C...  GET INDICES FOR TABLE-LOOKUP, L(*), AND INCREMENTS FOR
C     TAYLOR EXPANSION, F(*,3)
CDIR$    IVDEP
*VOCL LOOP,NOVREC
      DO 2 I=1,K
      L(I) = INT(RINC*F(I)+0.5)
    2 F(I+M3) = F(I)-XINC*DFLOAT(L(I))
C            -1 < L(*) < 161
C
C     COLLECT TABULATED VALUES OF THE INCOMPLETE GAMMA FUNCTION
C     AND PERFORM THE TAYLOR EXPANSION
C
      in=MIN0(IABS(2-N)+1,2)*M
      CALL GATHER(K,F(M+1), TAB(2,N+4+NACC), L)
      do 54 KA=nacc+3,2,-1
      aka = 1./dfloat(ka)
      call gather(k,f,tab(2,n+ka),l)
CDIR$ IVDEP
*VOCL LOOP,NOVREC
      do 53 i=1,k
   53 f(m+i) = f(i)-f(m+i)*f(m3+i)*aka
   54 continue
      call gather (k,f,tab(2,n+1),l)
CDIR$ IVDEP
*VOCL LOOP,NOVREC
      do 52 i=1,k
   52 f(in+i) = f(i) - f(m+i)*f(m3+i)
C
C...  EXPAND TO FULL SET OF T VALUES
      NNMM = N*M
      CALL ZERO (F(NNMM+1),M)
      CALL SCATTER (K,F(NNMM+1),FSCR,F(IN+1))
C...  DO DOWNWARD RECURSION
CDIR$ IVDEP
*VOCL LOOP,NOVREC
      DO 1 I=1,M
    1 FSCR(I) = 0.5*EXP(-FARG(I))
      DO 6 J=1,N
      NI=N-J+1
      CD=1./(DFLOAT(NI)-0.5)
CDIR$    IVDEP
*VOCL LOOP,NOVREC
      DO 5 I=1,M
    5 F(I+NNMM-M) = (F(I+NNMM)*FARG(I)+FSCR(I))*CD
    6 NNMM = NNMM - M
C...  ARE WE DONE ?
      IF(K.EQ.M) RETURN
      ENDIF
C
      IF (K .eq. 0) then
CDIR$ IVDEP
*VOCL LOOP,NOVREC
      DO 11 I=1,M
   11 FSCR(I) = 0.5*EXP(-FARG(I))
      ENDIF
C
C...  ASYMPTOTIC EXPANSION OF ZEROTH ORDER GAMMA FUNCTION
C     FOR LARGE T VALUES
C
C     Deal with values of T which are > TMAX
C
      DO 91 I=1,M
      IF(FARG(I) .LT. TMAX) GO TO 91
      F(I) = (0.5*SQRT(PI*FARG(I)) - FSCR(I))/FARG(I)
      DO 99 J=1,N
   99 F(I+J*M) = (F(I+J*M-M)*(J-0.5) - FSCR(I))/FARG(I)
   91 CONTINUE
C
      RETURN
C
C
      ENTRY VPTAG
C     INITIALIZES VALUES OF THE INCOMPLETE GAMMA FUNCTION
C     TO BE USED IN THE TAYLOR EXPANSION.
C     NNN=5  GENERATES VALUES FOR BASIS FUNCTIONS UP TO G-TYPE.
C
      PI = 2.*DACOS(0.0D00)
c
      DO 260 K=1,NNDIM
  260 TAB(1,K)=1.0/DFLOAT(2*K-1)
      DO 300 LL=1,NSIZ
      TT=XINC*LL
      AM1=NNDIM-0.5
      TERM=TAB(1,NNDIM)
      SUM=TERM
  280 AM1=AM1+1.
      TERM=TERM*TT/AM1
      SUM=SUM+TERM
      IF(TERM .GE. 1.E-16*SUM) GO TO 280
      COMEX= EXP(-TT)
      TAB(LL+1,NNDIM)=COMEX*SUM
      DO 290 K=1,NNDIM-1
  290 TAB(LL+1,NNDIM-K) =
     $ (TAB(LL+1,NNDIM+1-K)*2.*TT+COMEX)/(2*(NNDIM-K)-1.)
  300 CONTINUE

#ifdef _TRACE_VMOL
c      write(*,*) '@VMOL: Leaving file ',__FILE__
#endif

      RETURN
      END
