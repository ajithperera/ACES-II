      SUBROUTINE FNOK(D,C,F,NBAS,NSIZ1,BUF,IBUF,NINTMX,IA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 LABEL
      DIMENSION LABEL(52)
      DIMENSION D(NBAS,NBAS),F(NSIZ1,NSIZ1),C(NBAS,NBAS),IA(1)
      DIMENSION BUF(NINTMX),IBUF(NINTMX)
      COMMON/SWPPP/INP,IAMO,IFAMO,IINDO,IORTH                
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/INFOA/NBASIS,NUMSCF,NX,NMO2,NOC,NVT,NVO     
      DATA ZERO,ONE,TWO/0.D0,1.D0,2.D0/ 
C     WRITE(6,*) ' Now in FNOK '
      NTAPE=10
C  ...... F takes constant part from lower order .
C  ...... At the end, non-diagonal elements will
C  ...... be divided to half .......................
C  ...... you need this manipulation .............
C  ... Now this routine calculates only two-electron part  .....
      DO 1 I=1,NBASIS
      DO 1 J=1,NBASIS
    1 F(I,J) = ZERO
C  .............................................................
C     IF(IFAMO.EQ.0) THEN                      
C     OPEN(UNIT=NTAPE,FORM='UNFORMATTED',FILE='HF2',ACCESS='SEQUENTIAL')
C     READ(NTAPE) LABEL(1)
C     READ(NTAPE) LABEL(1)
C     READ(NTAPE) LABEL(1)
C     READ(NTAPE) LABEL(1)
C     ELSE
      OPEN(UNIT=NTAPE,FORM='UNFORMATTED',FILE='IIII',
     X ACCESS='SEQUENTIAL')
      CALL SEARCH('TWOELSUP',NTAPE)
C     END IF
C     WRITE(6,*) LABEL
C     WRITE(6,*) ' # , INTEGRAL I , J , K , L '
      NTOT=0
    5 READ(NTAPE,END=50) BUF,IBUF,NUT
      IF(NUT.LE.0) GO TO 30
      DO 10 IX=1,NUT
      NTOT = NTOT + 1
      L = AND(IBUF(IX),IALONE)
      K = AND(ISHFT(IBUF(IX),-IBITWD),IALONE)
      J = AND(ISHFT(IBUF(IX),-2*IBITWD),IALONE)
      I = AND(ISHFT(IBUF(IX),-3*IBITWD),IALONE)
      IF(J.GT.I) THEN
      IJE=I
      I=J
      J=IJE
      END IF
      IF(L.GT.K) THEN
      KLE=K
      K=L
      L=KLE
      END IF
C  I>J, K>L and (I,J)>(K,L) 
      VAL=BUF(IX)
      IF(I.EQ.J) VAL= VAL/2.D0
      IF(K.EQ.L) VAL= VAL/2.D0
      NIJ=IA(I)+J
      NKL=IA(K)+L
      IF(NKL.GT.NIJ) THEN
      IKE=I
      I=K
      K=IKE
      JLE=J
      J=L
      L=JLE
      END IF
      IF(NIJ.EQ.NKL) VAL= VAL/2.D0
C     write(6,*) I,J,K,L,VAL
      VAL4 =(VAL+VAL)+(VAL+VAL)
      VALD=VAL4*D(K,L)
      F(I,J)=F(I,J)+VALD
      F(J,I)=F(J,I)+VALD
      VALD=VAL4*D(I,J)
      F(K,L)=F(K,L)+VALD
      F(L,K)=F(L,K)+VALD
      VALD=VAL*D(J,L)
      F(I,K)=F(I,K)-VALD
      F(K,I)=F(K,I)-VALD
      VALD=VAL*D(J,K)
      F(I,L)=F(I,L)-VALD
      F(L,I)=F(L,I)-VALD
      VALD=VAL*D(I,L)
      F(J,K)=F(J,K)-VALD
      F(K,J)=F(K,J)-VALD
      VALD=VAL*D(I,K)
      F(J,L)=F(J,L)-VALD
      F(L,J)=F(L,J)-VALD
   10 CONTINUE
      GO TO  5
   30 CONTINUE
      DO 40  M=1,NBASIS
      DO 40  N=1,NBASIS
      F(M,N)=F(M,N)/TWO  
   40 CONTINUE  
C     CALL MATMUL(F,C,D,NBASIS,NBASIS,NBASIS,1,0)
      CALL XGEMM('N','N',NBASIS,NSIZ1,NBASIS,ONE,F,NBASIS,C,NBASIS
     X ,ZERO,D,NBASIS)
      CALL TRANSQ(C,NBASIS)
C     CALL MATMUL(C,D,F,NBASIS,NBASIS,NBASIS,1,0)
      CALL XGEMM('N','N',NSIZ1,NSIZ1,NBASIS,ONE,C,NBASIS,D,NBASIS
     X ,ZERO,F,NBASIS)
      CALL TRANSQ(C,NBASIS)
      CLOSE(NTAPE) 
      RETURN 
   50 CONTINUE
      WRITE(6,*) ' END of the file encountered ',NTAPE
      CLOSE(NTAPE)
      RETURN
      END
