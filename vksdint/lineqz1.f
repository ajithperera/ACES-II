      subroutine lineqz1(amat,xia,xianew,xupdate,asmall,
     &                  asquare,eval,conv,n,kmax,nocc)

      implicit double precision(a-h,o-z)
      INTEGER DIRPRD,POP,VRT
      DIMENSION DET(2)
      dimension amat(n,n),xia(n), xianew(n),xupdate(n,kmax),
     &          asmall(kmax,kmax),asquare(n),eval(1)
      common/machsp/iintln,ifltln,iintfp,ialone,ibitwd
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
c same problem as before below common block is stupid
c  sym2 is not compatible with the rest of common block 
       common /sym/ pop(8,2),vrt(8,2),nt(2),nfmi(2),nfea(2)
c      COMMON/SYM2/POP(8,2),VRT(8,2),NJUNK(6)
      
      data azero,one/0.d+0,1.0d+0/

      cutoff=conv*conv
      tol=azero
      call zero(asmall,kmax*kmax)
      anorm=sdot(n,xia,1,xia,1)

      ascale=sqrt(anorm)
      if(ascale/n.lt.cutoff) then
       write(6,*) ' wARNING FROM lINEQZ1 : the initial vector is zero.'
       ascale=one
      endif

      scale=one/ascale
      call sscal(n,scale,xia,1)
      write(*,*) 'xiascla'
c      call kkk(9,xia)    ! that seems right
      asquare(1)=one

      call dcopy(n,xia,1,xupdate(1,1),1)

      do 1000 k=1,kmax
      call xgemm('n','n',1,n,n,one,xia,1,amat,n,azero,xianew,1)
      call formz(xianew,eval,pop(1,1),vrt(1,1),nocc)
      DO 100 L=1,K
      ASMALL(L,K)=SDOT(N,XUPDATE(1,L),1,XIANEW,1)
100   CONTINUE
      IF(K.GT.1) ASMALL(K,K-1)=ASQUARE(K)
C
C    ORTHOGONALIZE XIANEW TO PREVIOUS VECTORS
C
      DO 200 L=1,K
      SCALE=-ASMALL(L,K)/ASQUARE(L)
      CALL SAXPY(N,SCALE,XUPDATE(1,L),1,XIANEW,1)
200   CONTINUE
C
      ASQUARE(K+1)=SDOT(N,XIANEW,1,XIANEW,1)
      CALL DCOPY(N,XIANEW,1,XIA,1)
      CALL DCOPY(N,XIANEW,1,XUPDATE(1,K+1),1)
C
C  CHECK FOR CONVERGENCE OF THE CREATED ITERATIVE SUBSPACE
C
      TEST=ASQUARE(K+1)/N
      IF(TEST.LE.CUTOFF) THEN
       ASQUARE(K+1)=ONE
C
C      ITERATIVE EXPANSION HAS CONVERGED, EXIT THE LOOP
C
       GO TO 300
      ENDIF
1000  CONTINUE
C
C  IF WE REACH THIS POINT, THE ITERATIVE EXPANSION HAS NOT CONVERGED
C
      WRITE(6,3000)
3000  FORMAT( '  The iterative expansion of D(ai) did not ',
     &           'converge, abort !')
      WRITE(6,3001) TEST
3001  FORMAT( '  The convergence is : ',E20.10)
      CALL ERREX
C
C   NOW INVERT THE ASMALL MATRIX
C
  300 CONTINUE
C
      WRITE(6,3002) K
3002  FORMAT('  The iterative expansion of D(ai)',
     &           ' converged after ',I3,' iterations.')
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
      DO 310 I=1,K
       ASMALL(I,I)=ASMALL(I,I)-ASQUARE(I)
310   CONTINUE
C
C  INVERT A (XIANEW IS HERE USED AS SCRATCH AND NO LONGER USED
C  FOR HOLDING A*Y
C
CSSS      CALL DGEFA(ASMALL,KMAX,K,XIANEW,I)
CSSS      IF (I.NE.0) THEN
CSSS         PRINT *, '@LINEQZ1: ASMALL cannot be inverted'
CSSS         CALL ERREX
CSSS      END IF
CSSS      CALL DGEDI(ASMALL,KMAX,K,XIANEW,0.d0,XIA,1)
C
C Eleminated above call to make the code compatible with 
C scs/mkl matheatical libraries, 08/08, Ajith Perera
C
      CALL MINV(ASMALL,K,KMAX,XIANEW,DET,AZERO,0,1)
C
C     NOW FORM SOLUTION IN XIA
C
      CALL DGEMV('N',N,K,
     &           -ASCALE,XUPDATE,N,
     &                   ASMALL, 1,
     &           0.d0,   XIA,    1)
C
C   ALL DONE, RETURN
C
      RETURN
      END
