      SUBROUTINE TRANS513(A, M, N, MN, MOVE, IWRK, IOK)
      double precision A(MN)
      integer MOVE(IWRK)

      IF (M.LT.2 .OR. N.LT.2) THEN
         IOK = 0
         RETURN
      END IF
      IF (MN.NE.M*N) THEN
         IOK = -1
         RETURN
      END IF
      IF (IWRK.LT.1) THEN
         IOK = -2
         RETURN
      END IF

      IF (M.EQ.N) THEN
      DO I=1,n-1
        DO J=i+1,N
          B              = A(I + (J-1)*N)
          A(I + (J-1)*N) = A(J + (I-1)*N)
          A(J + (I-1)*N) = B
        END DO
      END DO
      IOK = 0
      RETURN
      END IF

      NCOUNT = 2
      K = MN - 1
      DO I=1,IWRK
         MOVE(I) = 0
      END DO

      IF (M.GE.3 .AND. N.GE.3) THEN
      IR2 = M - 1
      IR1 = N - 1
      IR0 = 1
      DO WHILE (IR0.NE.0)
         IR0 = MOD(IR2,IR1)
         IR2 = IR1
         IR1 = IR0
      END DO
      NCOUNT = NCOUNT + IR2 - 1
      END IF

      I = 1
      IM = M
      GO TO 80
   40 CONTINUE
      MAX = K - I
      I = I + 1
      IF (I.GT.MAX) THEN
         IOK = I
         RETURN
      END IF
      IM = IM + M
      IF (IM.GT.K) IM = IM - K
      I2 = IM
      IF (I.EQ.I2) GO TO 40
      IF (I.GT.IWRK) THEN
         IF (I2.LE.I .OR. I2.GE.MAX) THEN
            IF (I2.NE.I) GO TO 40
         ELSE
            GO TO 80
         END IF
         I1 = I2
         DO
            I2 = M*I1 - K*(I1/N)
            IF (I2.LE.I .OR. I2.GE.MAX) THEN
               IF (I2.NE.I) GO TO 40
            ELSE
               GO TO 80
            END IF
            I1 = I2
         END DO
      END IF
      IF (MOVE(I).NE.0) THEN
         GO TO 40
      ELSE
         GO TO 80
      END IF
   80 CONTINUE
      I1 = I
      KMI = K - I
      B = A(I1+1)
      I1C = KMI
      C = A(I1C+1)
      DO
         I2 = M*I1 - K*(I1/N)
         I2C = K - I2
         IF (I1.LE.IWRK)  MOVE(I1)  = 2
         IF (I1C.LE.IWRK) MOVE(I1C) = 2
         NCOUNT = NCOUNT + 2
         IF (I2.EQ.I) THEN
            A(I1+1) = B
            A(I1C+1) = C
            IF (NCOUNT.LT.MN) GO TO 40
            IOK = 0
            RETURN
         END IF
         IF (I2.EQ.KMI) THEN
            D = B
            B = C
            C = D
            A(I1+1) = B
            A(I1C+1) = C
            IF (NCOUNT.LT.MN) GO TO 40
            IOK = 0
            RETURN
         END IF
         A(I1+1) = A(I2+1)
         A(I1C+1) = A(I2C+1)
         I1 = I2
         I1C = I2C
      END DO
      END

