      SUBROUTINE TRANS380(A,M,N,MN,MOVE,IWRK,IOK)
      double precision A(MN)
      integer MOVE(IWRK)
      IF(M.LT.2.OR.N.LT.2) GO TO 60
      IF(MN.NE.M*N) GO TO 92
      IF(IWRK.LT.1) GO TO 93
      IF(M.EQ.N) GO TO 70
      NCOUNT=2
      M2=M-2
      DO 10 I=1,IWRK
   10   MOVE(I)=0
      IF(M2.LT.1) GO TO 12
      DO 11 IA=1,M2
        IB=IA*(N-1)/(M-1)
        IF(IA*(N-1).NE.IB*(M-1)) GO TO 11
        NCOUNT=NCOUNT+1
        I=IA*N+IB
        IF(I.GT.IWRK) GO TO 11
        MOVE(I)=1
   11   CONTINUE
   12 K=MN-1
      KMI=K-1
      MAX=MN
      I=1
      GO TO 30
   20        MAX=K-I
             I=I+1
             KMI=K-I
             IF(I.GT.MAX) GO TO 90
             IF(I.GT.IWRK) GO TO 21
             IF(MOVE(I).LT.1) GO TO 30
             GO TO 20
   21        IF(I.EQ.M*I-K*(I/N)) GO TO 20
             I1=I
   22        I2=M*I1-K*(I1/N)
             IF(I2.LE.I .OR. I2.GE.MAX)  GO TO 23
             I1=I2
             GO TO 22
   23        IF(I2.NE.I) GO TO 20
   30             I1=I
   31             B=A(I1+1)
   32             I2=M*I1-K*(I1/N)
                  IF(I1.LE.IWRK) MOVE(I1)=2
   33             NCOUNT=NCOUNT+1
                  IF(I2.EQ.I .OR. I2.GE.KMI) GO TO 35
   34             A(I1+1)=A(I2+1)
                  I1=I2
                  GO TO 32
   35             IF(MAX.EQ.KMI .OR. I2.EQ.I) GO TO 41
                  MAX=KMI
                  GO TO 34
   41        A(I1+1)=B
             IF(NCOUNT.GE.MN) GO TO 60
             IF(I2.EQ.MAX .OR. MAX.EQ.KMI) GO TO 20
             MAX=KMI
             I1=MAX
             GO TO 31
   60 IOK=0
      RETURN
   70 N1=N-1
      DO 71 I=1,N1
        J1=I+1
        DO 71 J=J1,N
          I1=I+(J-1)*N
          I2=J+(I-1)*M
          B=A(I1)
          A(I1)=A(I2)
          A(I2)=B
   71     CONTINUE
      GO TO 60
   90 IOK=I
   91 RETURN
   92 IOK=-1
      GO TO 91
   93 IOK=-2
      GO TO 91
      END
