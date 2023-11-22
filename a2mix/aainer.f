      FUNCTION AAINER(R,S,T,L,M,N,D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION D(*)
      COMMON /VARS/FNU(20),DP(20),FN(15),FD(15),PI,DFTR(11)
      INTEGER R,S,T,U,V,W,UVWT,RSTT,UVWTH
      RSTT=R+S+T
      LMNT=L+M+N
      LMNRST=RSTT+LMNT
      LH=L/2
      MH=M/2
      NH=N/2
      AAINER=0.0
      LAST1 = LH+1
      LAST2 = MH + 1
      LAST3 = NH + 1
      DO 6 II=1,LAST1
      I = II-1
      IL=L-2*I+1
      ILR=R+IL
      ILRH=(ILR-1)/2
      FI=FN(ILR)*FD(IL)*FD(I+1)
      LAST4 = ILRH + 1
      DO 5 JJ=1,LAST2
      J = JJ-1
      JM=M-2*J+1
      JMS=S+JM
      JMSH=(JMS-1)/2
      FIJ=FN(JMS)*FD(JM)*FD(J+1)*FI
      LAST5 = JMSH + 1
      DO 4 KK=1,LAST3
      K = KK-1
      KN=N-2*K+1
      KNT=T+KN
      KNTH=(KNT-1)/2
      IJKT=I+J+K
      FIJK=FN(KNT)*FD(KN)*FD(K+1)*DP(IJKT+1)*FIJ
      LMRSIJ=LMNRST-2*IJKT
      LAST6 = KNTH + 1
      DO 3 IU=1,LAST4
      U = IU-1
      ILRU=ILR-2*U
      FU=FD(U+1)*FD(ILRU)
      IF( ABS(D(1))-1.E-10) 10,10,11
10    IF(ILRU-1) 12,12, 3
11    FU=FU*D(1)**(ILRU-1)
12    CONTINUE
      DO 2 IV=1,LAST5
      V = IV-1
      JMSV=JMS-2*V
      FUV=FU*FD(V+1)*FD(JMSV)
      IF( ABS(D(2))-1.E-10) 20,20,21
20    IF(JMSV-1) 22,22,2
21    FUV=FUV*D(2)**(JMSV-1)
22    CONTINUE
      DO 1 IW=1,LAST6
      W = IW-1
      KNTW=KNT-2*W
      FUVW=FUV*FD(W+1) *FD(KNTW)
      IF( ABS(D(3))-1.E-10) 30,30,31
30    IF(KNTW-1) 32,32,1
31    FUVW=FUVW*D(3)**(KNTW-1)
32    UVWT=U+V+W
      UVWTH=UVWT/2
      IF(2*UVWTH-UVWT) 33,34,33
33    FUVW=-FUVW
34    NUINDX=LMRSIJ-UVWT
      FUVW=FUVW*FNU(NUINDX  +1)*DP(UVWT+1)
      AAINER=FIJK*FUVW+AAINER
1     CONTINUE
2     CONTINUE
3     CONTINUE
4     CONTINUE
5     CONTINUE
6     CONTINUE
      RETURN
      END
