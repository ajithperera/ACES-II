      SUBROUTINE CALCZETA(SCR,MAXCOR,IUHF,ASMALL,ASQUARE,KMAX,NSIZEC,
     &                    LISTZ1,LISTZ1OFF,LISTZ2)
C
C THIS ROUTINE SOLVES THE LINEAR EQUATION FOR THE ZETA AMPLITUDES
C NEEDED IN EOM-CC GRADIENTS:
C
C
C   <Q|exp(-T) H exp(T)|Q><Q|ZETA|0> = <Q|XI|0>
C
C USING A STANDARD REDUCED SUBSPACE ALGORITHM.  THIS IS HACKED FROM
C THE ROUTINE CALCZETA IN THE VCCEH PROGRAM.
C
C LIST 1,470 - OLD EXPANSION (E) VECTORS, ONE VECTOR PER LOGICAL RECORD
C LIST 2,470 - OLD H*E VECTORS, ONE VECTOR PER LOGICAL RECORD
C LIST 1,471 - M VECTOR
C LIST 1,88 - DIAGONAL PART OF Hd, STORED ON ONE LOGICAL RECORD
C
C HERE H IS PARTITIONED INTO [Hd + D] WHERE Hd IS A DIAGONAL
C PART OF H AND D IS THE REMAINDER OF THE MATRIX
C
C EQUATION TO BE SOLVED FOR Y IS:
C
C         H * Y = M
C
C      [Hd + D]Y = M
C
C    [1 + Hd(-1)*D] Y - Hd(-1) * M = 0
C
C THE LAST EQUATION IS THE FORM USED IN THIS ROUTINE.
C
C SUBSPACE BUILT UP BY SUCCESSIVE MULTIPLICATION OF TRIAL VECTORS BY
C Hd(-1)*D FOLLOWED BY ORTHOGONALIZATION
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL PRINT,PRINT2
C
      DIMENSION SCR(MAXCOR),ASMALL(KMAX,KMAX),ASQUARE(KMAX)
C
      integer IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/FLAGS/IFLAGS(100)
C
      DATA ONE /1.0D0/
C
      PRINT2=IFLAGS(1).GE.100
      PRINT=IFLAGS(1).GE.0
C
      CUTOFF=10.0D0**(-IFLAGS(88))
C
      I000=1
      I010=I000+NSIZEC
      I020=I010+NSIZEC
      I030=I020+NSIZEC
C
C PUT Hd ELEMENTS ON LIST 88
C
c      CALL HBARDIAG(1,SCR,MAXCOR,IUHF)
C
C PUT EVENTUAL SOLUTION VECTOR ON LIST 471
C
      CALL GETLST(SCR(I000),1,1,1,1,474)
      CALL PUTLST(SCR(I000),1,1,1,1,471)
C
C INITIATE CONSTRUCTION OF SOLUTION SPACE
C
      CALL ZERO(ASMALL,KMAX*KMAX)
C
C SCALE INPUT VECTOR WITH DIAGONAL ELEMENTS OF HBAR - Hd-1 * M
C
      CALL GETLST(SCR(I010),1,1,1,1,88)
      DO 901 I=1,NSIZEC
       SCR(I000-1+I)=SCR(I000-1+I)/(SCR(I010-1+I))
901   CONTINUE
C
C NORMALIZE VECTOR - THIS IS FIRST ELEMENT OF EXPANSION SPACE
C
      CALL SCOPY(NSIZEC,SCR(I000),1,SCR(I010),1)
      IF(IUHF.EQ.0)THEN
       CALL SPNTSING(1,SCR(I010),SCR(I020),MAXCOR-I010+1)
      ENDIF
      ZNORM=SQRT(SDOT(NSIZEC,SCR(I000),1,SCR(I010),1))
      X0=ONE/ZNORM
      CALL SSCAL(NSIZEC,X0,SCR(I000),1)
      CALL UPDATES(1,SCR(I000),444,0,490,IUHF)
      CALL PUTLST(SCR(I000),1,1,1,1,470)
C
      ASQUARE(1)=ONE
C
C BUILD REMAINING EXPANSION SPACE
C
      DO 1000 K=1,KMAX
C
C CALCULATE [Hd(-1) * D]*C(k) = Hd(-1) * [H*C(k) - Hd*C(k)]
C
       CALL HBARXC(SCR(I000),MAXCOR*IINTFP,IUHF,2,1)
       CALL FETCHVEC(1,SCR(I000),IUHF,2,490,460)
C
       CALL GETLST(SCR(I010),1,1,1,1,88)
       CALL GETLST(SCR(I020),K,1,1,1,470)
C
       DO 902 I=1,NSIZEC
        SCR(I000-1+I)=(SCR(I000-1+I)-SCR(I010-1+I)*SCR(I020-1+I))/
     &                 SCR(I010-1+I)
902    CONTINUE
C
       CALL PUTLST(SCR(I000),K,1,1,2,470)
C
       DO 100 I=1,K
        CALL GETLST(SCR(I000),I,1,1,1,470)
        DO 101 J=1,K
         CALL GETLST(SCR(I010),J,1,1,2,470)
         IF(IUHF.EQ.0)THEN
          CALL SPNTSING(1,SCR(I010),SCR(I020),MAXCOR-I010+1)
         ENDIF
         ASMALL(I,J)=SDOT(NSIZEC,SCR(I000),1,SCR(I010),1)
101     CONTINUE
100    CONTINUE
C
C ORTHOGONALIZE H'*C(k) TO EXISTING SUBSPACE
C
       CALL GETLST(SCR(I010),K,1,1,2,470)
       DO 102 I=1,K
        CALL GETLST(SCR(I000),I,1,1,1,470)
        X=-ASMALL(I,K)/ASQUARE(I)
        CALL SAXPY(NSIZEC,X,SCR(I000),1,SCR(I010),1)
102    CONTINUE
C
       CALL SCOPY(NSIZEC,SCR(I010),1,SCR(I020),1)
       IF(IUHF.EQ.0)THEN
        CALL SPNTSING(1,SCR(I020),SCR(I030),MAXCOR-I020+1)
       ENDIF
       ASQUARE(K+1)=SDOT(NSIZEC,SCR(I010),1,SCR(I020),1)
       RESID=ASQUARE(K+1)
C
       IF(PRINT)THEN
        WRITE(6,5002)K,RESID
5002    FORMAT(T3,'@CALCZETA-I, Residual after ',I5,
     &         ' iterations is ',D15.10,'.')
       ENDIF
       CALL PUTLST(SCR(I010),K+1,1,1,1,470)
       CALL UPDATES(1,SCR(I010),444,0,490,IUHF)
       IF(RESID.LT.CUTOFF)THEN
        ASQUARE(K+1)=ONE
        WRITE(6,5003)K
5003    FORMAT(T3,'@CALCZETA-I, Derivative amplitudes converged',
     &            ' after ',I5,
     &            ' iterations.')
        GOTO 300
       ENDIF
C
1000  CONTINUE
C
      GOTO 999
C
300   WRITE(6,5000)
5000  FORMAT(T3,'@CALCZETA-I, Solving linear equation in reduced ',
     &          'subspace.')
C
C FIRST WE MUST MODIFY SUBSPACE MATRIX, WHICH IS CURRENTLY A REPRESENTATION
C OF THE OPERATOR [Hd(-1) * D] TO [1 + Hd(-1) * D].  THIS IS SIMPLY
C ACCOMPLISHED BY ADDING THE NORM OF THE APPROPRIATE EXPANSION VECTORS
C TO THE DIAGONAL ELEMENTS.
C
      CALL SAXPY(K,ONE,ASQUARE,1,ASMALL,KMAX+1)
C
C INVERT SUBSPACE REPRESENTATION OF [1 + Hd(-1) * D]
C
      CALL MINV(ASMALL,K,KMAX,SCR(I000),DET,TOL,0,1)
C
C GENERATE SOLUTION VECTOR AS LINEAR COMBINATION OF BASIS VECTORS
C
      CALL ZERO(SCR(I000),NSIZEC)
      DO 340 I=1,K
       CALL GETLST(SCR(I010),I,1,1,1,470)
       CALL SAXPY (NSIZEC,ASMALL(I,1)/X0,SCR(I010),1,SCR(I000),1)
340   CONTINUE
C
C WRITE SOLUTION VECTOR TO DISK - OVERWRITES LAMBDA!
C
      CALL VMINUS(SCR(I000),NSIZEC)
      CALL PUTLST(SCR(I000),1,1,1,1,474)
      call UPDATES(1,scr(i000),LISTZ2,LISTZ1OFF,LISTZ1,IUHF)
C
      IF(PRINT2)THEN
       WRITE(6,1001)
1001   FORMAT(T3,'@CALCZETA-I, Derivative T amplitudes : ')
       CALL PRVECR(SCR(I000),NSIZEC)
      ENDIF
C
      RETURN
999   WRITE(6,5001)KMAX
      CALL ERREX
5001  FORMAT(T3,'@CALCZETA-F, Convergence failed after ',I5,
     &          ' iterations.')
      RETURN
      END
