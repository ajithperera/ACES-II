C
      SUBROUTINE OVPP(AA,X1,X2,E,I2,I1,GAMA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /VARS/ FNU(20),DP(20),FN(15),FD(15),PI,DFTR(11)
      DIMENSION S(20),S2(20)
      AX=E-X1
      GO TO (1,2,3,4,5,6,7,8,9),I1
1     CONTINUE
      S(1)=1.0
      GO TO 10
2     CONTINUE
      S(1)=AX
      S(2)=1.0
      GO TO 10
3     CONTINUE
      S(1)=AX**2
      S(2)=2.*AX
      S(3)=1.
      GO TO 10
4     CONTINUE
      S(1)=AX**3
      S(2)=3.*AX**2
      S(3)=3.*AX
      S(4)=1.
      GO TO 10
5     CONTINUE
      S(1)=AX**4
      S(2)=4.*AX**3
      S(3)=6.*AX**2
      S(4)=4.*AX
      S(5)=1.
      GO TO 10
6     CONTINUE
      S(1)=AX**5
      S(2)=5.*AX**4
      S(3)=10.*AX**3
      S(4)=10.*AX**2
      S(5)=5.*AX
      S(6)=1.
      GO TO 10
7     CONTINUE
      S(1)=AX**6
      S(2)=6.*AX**5
      S(3)=15.*AX**4
      S(4)=20.*AX**3
      S(5)=15.*AX**2
      S(6)=6.*AX
      S(7)=1.
      GOTO 10
8     CONTINUE
      S(1) = AX**7
      S(2) = 7.*AX**6
      S(3) = 21.*AX**5
      S(4) = 35.*AX**4
      S(5) = 35.*AX**3
      S(6) = 21.*AX**2
      S(7) = 7.*AX
      S(8) = 1.
      GOTO 10
9     CONTINUE
      S(1) = AX**8
      S(2) = 8.*AX**7
      S(3) = 28.*AX**6
      S(4) = 56.*AX**5
      S(5) = 70.*AX**4
      S(6) = 56.*AX**3
      S(7) = 28.*AX**2
      S(8) = 8.*AX
      S(9) = 1.
10    CONTINUE
      AX=E-X2
      GO TO (11,12,13,14,15,16,17,18,19),I2
11    CONTINUE
      S2(1)=1.0
      GO TO 20
12    CONTINUE
      S2(1)=AX
      S2(2)=1.0
      GO TO 20
13    CONTINUE
      S2(1)=AX**2
      S2(2)=2.*AX
      S2(3)=1.
      GO TO 20
14    CONTINUE
      S2(1)=AX**3
      S2(2)=3.*AX**2
      S2(3)=3.*AX
      S2(4)=1.
      GO TO 20
15    CONTINUE
      S2(1)=AX**4
      S2(2)=4.*AX**3
      S2(3)=6.*AX**2
      S2(4)=4.*AX
      S2(5)=1.
      GO TO 20
16    CONTINUE
      S2(1)=AX**5
      S2(2)=5.*AX**4
      S2(3)=10.*AX**3
      S2(4)=10.*AX**2
      S2(5)=5.*AX
      S2(6)=1.
      GO TO 20
17    CONTINUE
      S2(1)=AX**6
      S2(2)=6.*AX**5
      S2(3)=15.*AX**4
      S2(4)=20.*AX**3
      S2(5)=15.*AX**2
      S2(6)=6.*AX
      S2(7)=1.
      GOTO 20
18    CONTINUE
      S2(1) = AX**7
      S2(2) = 7.*AX**6
      S2(3) = 21.*AX**5
      S2(4) = 35.*AX**4
      S2(5) = 35.*AX**3
      S2(6) = 21.*AX**2
      S2(7) = 7.*AX
      S2(8) = 1.
      GOTO 20
19    CONTINUE
      S2(1) = AX**8
      S2(2) = 8.*AX**7
      S2(3) = 28.*AX**6
      S2(4) = 56.*AX**5
      S2(5) = 70.*AX**4
      S2(6) = 56.*AX**3
      S2(7) = 28.*AX**2
      S2(8) = 8.*AX
      S2(9) = 1.
20    CONTINUE
      AA=0.0
      DO 30 I=1,I1
      DO 31 J=1,I2
      NH=(I+J)/2
      IF(NH*2.NE.(I+J))GO TO 31
      AA=AA+(.5/GAMA)**(NH-1)*DFTR(NH)*S(I)*S2(J)
31    CONTINUE
30    CONTINUE
      RETURN
      END
