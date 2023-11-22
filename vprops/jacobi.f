      SUBROUTINE JACOBI(A,B,C,N)
C.10/09/72
C     DIAGONALIZE SYMMETRIC MATRIX BY THE JACOBI METHOD
C     A = LOWER TRIANGLE OF GIVEN SYMMETRIC MATRIX
C     B = MATRIX OF EIGENVECTORS
C     C = ARRAY OF EIGENVALUES
C     N = ORDER OF MATRIX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1),B(1),C(1)
      DATA EPS/1.0E-8/,LOOPMX/5000/
      LOOPC=0
      NA=N
      NN=(NA*(NA+1))/2
      K=1
      DO 115 I=1,NA
      DO 115 J=1,NA
      B(K)=0.
      IF (I.EQ.J) B(K)=1.0
115   K=K+1
      SUM=0.
      NNA=NA-1
      IF (NNA) 325,310,125
125   K=1
      AMAX=0.
      DO 155 I=1,NA
      DO 150 J=1,I
      IF (I.EQ.J) GO TO 145
      IF (ABS(A(K)).GT.AMAX) AMAX=ABS(A(K))
145   TERM=A(K)*A(K)
      SUM=SUM+TERM+TERM
150   K=K+1
155   SUM=SUM-TERM
      SUM=SQRT(SUM)
      THRESH=SUM/SQRT(DFLOAT(NA))
      THRSHG=THRESH*EPS
      IF(THRSHG.GE.AMAX) GO TO 310
      THRESH=AMAX/3.
      IF(THRESH.LT.THRSHG) THRESH=THRSHG
180   K=2
      M=0
      JD=1
      KJ=0
      DO 270 J=2,NA
      KJ=KJ+NA
      JD=JD+J
      JJ=J-1
      JJJ=JJ-1
      ID=0
      KI=-NA
      DO 265 I=1,JJ
      KI=KI+NA
      ID=ID+I
      IF (ABS(A(K)).LE.THRESH) GO TO 265
      M=M+8
      ALPHA=(A(JD)-A(ID))/(A(K)+A(K))
      BETA=0.25/(1.0+ALPHA**2)
      ROOT=0.5+ABS(ALPHA)*SQRT(BETA)
      IF (ALPHA.GE.0.) GO TO 205
      CSQ=BETA/ROOT
      SSQ=ROOT
      GO TO 210
205   SSQ=BETA/ROOT
      CSQ=ROOT
210   CC=SQRT(CSQ)
      S=-SQRT(SSQ)
      TWOSC=CC*(S+S)
      TEMPA=CSQ*A(ID)+TWOSC*A(K)+SSQ*A(JD)
      A(JD)=CSQ*A(JD)-TWOSC*A(K)+SSQ*A(ID)
      A(ID)=TEMPA
      A(K)=0.
      KA=JD-J
      KB=ID-I
      KC=KI
      KD=KJ
      II=I-1
      IF (II.EQ.0) GO TO 230
      DO 220 L=1,II
      KC=KC+1
      KD=KD+1
      TEMPA=CC*B(KC)+S*B(KD)
      B(KD)=-S*B(KC)+CC*B(KD)
      B(KC)=TEMPA
      KB=KB+1
      KA=KA+1
      TEMPA=CC*A(KB)+S*A(KA)
      A(KA)=-S*A(KB)+CC*A(KA)
220   A(KB)=TEMPA
230   KC=KC+1
      KD=KD+1
      TEMPA=CC*B(KC)+S*B(KD)
      B(KD)=-S*B(KC)+CC*B(KD)
      B(KC)=TEMPA
      KB=KB+1
      KA=KA+1
      IF (I.EQ.JJ) GO TO 250
      DO 240 L=I,JJJ
      KC=KC+1
      KD=KD+1
      TEMPA=CC*B(KC)+S*B(KD)
      B(KD)=-S*B(KC)+CC*B(KD)
      B(KC)=TEMPA
      KB=KB+L
      KA=KA+1
      TEMPA=CC*A(KB)+S*A(KA)
      A(KA)=-S*A(KB)+CC*A(KA)
240   A(KB)=TEMPA
250   KC=KC+1
      KD=KD+1
      TEMPA=CC*B(KC)+S*B(KD)
      B(KD)=-S*B(KC)+CC*B(KD)
      B(KC)=TEMPA
      KB=KB+JJ
      KA=KA+1
      IF (J.EQ.NA) GO TO 265
      DO 260 L=J,NNA
      KC=KC+1
      KD=KD+1
      TEMPA=CC*B(KC)+S*B(KD)
      B(KD)=-S*B(KC)+CC*B(KD)
      B(KC)=TEMPA
      KB=KB+L
      KA=KA+L
      TEMPA=CC*A(KB)+S*A(KA)
      A(KA)=-S*A(KB)+CC*A(KA)
260   A(KB)=TEMPA
265   K=K+1
270   K=K+1
      LOOPC=LOOPC+1
      IF(LOOPC.GT.LOOPMX) GO TO 330
      IF (M.GT.NN) GO TO 180
      IF(THRESH.EQ.THRSHG) GO TO 300
      THRESH=THRESH/3.
      IF(THRESH.GE.THRSHG) GO TO 180
      THRESH=THRSHG
      GO TO 180
300   IF (M.NE.0) GO TO 180
310   LL=0
      DO 320 L=1,NA
      LL=LL+L
320   C(L)=A(LL)
325   RETURN
330   PRINT 335,LOOPMX
335   FORMAT (/18H NO CONVERGENCE IN,I4,11H ITERATIONS)
      STOP
      END
