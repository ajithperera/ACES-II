      SUBROUTINE TRANS513(A, M, N, MN, MOVE, IWRK, IOK)
      DIMENSION A(MN), MOVE(IWRK)
      IF (M.LT.2 .OR. N.LT.2) GO TO 120
      IF (MN.NE.M*N) GO TO 180
      IF (IWRK.LT.1) GO TO 190
      IF (M.EQ.N) GO TO 130
      NCOUNT = 2
      K = MN - 1
      DO 10 I=1,IWRK
        MOVE(I) = 0
   10 CONTINUE
      IF (M.LT.3 .OR. N.LT.3) GO TO 30
      IR2 = M - 1
      IR1 = N - 1
   20 IR0 = MOD(IR2,IR1)
      IR2 = IR1
      IR1 = IR0
      IF (IR0.NE.0) GO TO 20
      NCOUNT = NCOUNT + IR2 - 1
   30 I = 1
      IM = M
      GO TO 80
   40 MAX = K - I
      I = I + 1
      IF (I.GT.MAX) GO TO 160
      IM = IM + M
      IF (IM.GT.K) IM = IM - K
      I2 = IM
      IF (I.EQ.I2) GO TO 40
      IF (I.GT.IWRK) GO TO 60
      IF (MOVE(I).EQ.0) GO TO 80
      GO TO 40
   50 I2 = M*I1 - K*(I1/N)
   60 IF (I2.LE.I .OR. I2.GE.MAX) GO TO 70
      I1 = I2
      GO TO 50
   70 IF (I2.NE.I) GO TO 40
   80 I1 = I
      KMI = K - I
      B = A(I1+1)
      I1C = KMI
      C = A(I1C+1)
   90 I2 = M*I1 - K*(I1/N)
      I2C = K - I2
      IF (I1.LE.IWRK) MOVE(I1) = 2
      IF (I1C.LE.IWRK) MOVE(I1C) = 2
      NCOUNT = NCOUNT + 2
      IF (I2.EQ.I) GO TO 110
      IF (I2.EQ.KMI) GO TO 100
      A(I1+1) = A(I2+1)
      A(I1C+1) = A(I2C+1)
      I1 = I2
      I1C = I2C
      GO TO 90
  100 D = B
      B = C
      C = D
  110 A(I1+1) = B
      A(I1C+1) = C
      IF (NCOUNT.LT.MN) GO TO 40
  120 IOK = 0
      RETURN
  130 N1 = N - 1
      DO 150 I=1,N1
        J1 = I + 1
        DO 140 J=J1,N
          I1 = I + (J-1)*N
          I2 = J + (I-1)*M
          B = A(I1)
          A(I1) = A(I2)
          A(I2) = B
  140   CONTINUE
  150 CONTINUE
      GO TO 120
  160 IOK = I
  170 RETURN
  180 IOK = -1
      GO TO 170
  190 IOK = -2
      GO TO 170
      END
