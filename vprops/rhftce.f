      SUBROUTINE RHFTCE(C,A,E,ITYP,ITM,Q,LMN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE IP
      DIMENSION LMN(27),IP(20)
      DIMENSION FAC(9,9),XH(9),YH(9),ZH(9)
      COMMON /HIGHL/ LMNVAL(3,84), ANORM(84)
      DIMENSION C(27)
      DIMENSION A(3),E(3)
      IBTAND(I,J) = IAND(I,J)
      IBTOR(I,J)  = IOR(I,J)
      IBTXOR(I,J) = IEOR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTNOT(I)   = NOT(I)
      AEX = E(1)-A(1)
      AEY = E(2)-A(2)
      AEZ = E(3)-A(3)
      IF(ITYP.GT.20)GO TO 40
C
C....    GOTO LIST IS REARRANGED TO ACCOMMODATE DICTIONARY ORDER
C
      GOTO (11, 12,13,14,
     1      15,18,19,16,20,17,
     2      21,24,25,26,30,28,22,27,29,23),ITYP
11    CONTINUE
C     ITYP=1  S
      C(1)=Q
      LMN(1)=IP(1)
      ITM=1
      GO TO 100
12    CONTINUE
C     ITYP=2 X
      C(1)=AEX*Q
      C(2)=Q
      ITM=2
      LMN(1)=IP(1)
      LMN(2)=IP(2)
      GO TO 100
13    CONTINUE
C     ITYP=3 Y
      C(1)=AEY*Q
      C(2)=Q
      ITM=2
      LMN(1)=IP(1)
      LMN(2)=IP(3)
      GO TO 100
14    CONTINUE
C     ITYP=4 Z
      C(1)=AEZ*Q
      C(2)=Q
      ITM=2
      LMN(1)=IP(1)
      LMN(2)=IP(4)
      GO TO 100
15    CONTINUE
C     ITYP=5 XX
      C(1)=AEX*AEX*Q
      C(2)=2.0*AEX*Q
      C(3)=Q
      ITM=3
      LMN(1)=IP(1)
      LMN(2)=IP(2)
      LMN(3)=IP(5)
      GO TO 100
16    CONTINUE
C     ITYP=8 YY
      C(1)=AEY*AEY*Q
      C(2)=2.0*AEY*Q
      C(3)=Q
      ITM=3
      LMN(1)=IP(1)
      LMN(2)=IP(3)
      LMN(3)=IP(6)
      GO TO 100
17    CONTINUE
C     ITYP=10 ZZ
      C(1)=AEZ*AEZ*Q
      C(2)=2.0*AEZ*Q
      C(3)=Q
      ITM=3
      LMN(1)=IP(1)
      LMN(2)=IP(4)
      LMN(3)=IP(7)
      GO TO 100
18    CONTINUE
C     ITYP=6 XY
      C(1)=AEX*AEY*Q
      C(2)=AEY*Q
      C(3)=AEX*Q
      C(4)=Q
      ITM=4
      LMN(1)=IP(1)
      LMN(2)=IP(2)
      LMN(3)=IP(3)
      LMN(4)=IP(8)
      GO TO 100
19    CONTINUE
C     ITYP=7 XZ
      C(1)=AEX*AEZ*Q
      C(2)=AEZ*Q
      C(3)=AEX*Q
      C(4)=Q
      ITM=4
      LMN(1)=IP(1)
      LMN(2)=IP(2)
      LMN(3)=IP(4)
      LMN(4)=IP(9)
      GO TO 100
20    CONTINUE
C     ITYP=9 YZ
      C(1)=AEY*AEZ*Q
      C(2)=AEZ*Q
      C(3)=AEY*Q
      C(4)=Q
      ITM=4
      LMN(1)=IP(1)
      LMN(2)=IP(3)
      LMN(3)=IP(4)
      LMN(4)=IP(10)
      GO TO 100
21    CONTINUE
C     ITYP=11 XXX
      C(1)=AEX*AEX*AEX*Q
      C(2)=3.0*AEX*AEX*Q
      C(3)=3.0*AEX*Q
      C(4)=Q
      ITM=4
      LMN(1)=IP(1)
      LMN(2)=IP(2)
      LMN(3)=IP(5)
      LMN(4)=IP(11)
      GO TO 100
22    CONTINUE
C     ITYP=17 YYY
      C(1)=AEY*AEY*AEY*Q
      C(2)=3.0*AEY*AEY*Q
      C(3)=3.0*AEY*Q
      C(4)=Q
      ITM=4
      LMN(1)=IP(1)
      LMN(2)=IP(3)
      LMN(3)=IP(6)
      LMN(4)=IP(12)
      GO TO 100
23    CONTINUE
C     ITYP=20 ZZZ
      C(1)=AEZ*AEZ*AEZ*Q
      C(2)=3.0*AEZ*AEZ*Q
      C(3)=3.0*AEZ*Q
      C(4)=Q
      ITM=4
      LMN(1)=IP(1)
      LMN(2)=IP(4)
      LMN(3)=IP(7)
      LMN(4)=IP(13)
      GO TO 100
24    CONTINUE
C     ITYP=12 XXY
      C(1)=AEX*AEX*AEY*Q
      C(2)=2.0*AEX*AEY*Q
      C(3)=AEX*AEX*Q
      C(4)=AEY*Q
      C(5)=2.0*AEX*Q
      C(6)=Q
      ITM=6
      LMN(1)=IP(1)
      LMN(2)=IP(2)
      LMN(3)=IP(3)
      LMN(4)=IP(5)
      LMN(5)=IP(8)
      LMN(6)=IP(14)
      GO TO 100
25    CONTINUE
C     ITYP=13 XXZ
      C(1)=AEX*AEX*AEZ*Q
      C(2)=2.0*AEX*AEZ*Q
      C(3)=AEX*AEX*Q
      C(4)=AEZ*Q
      C(5)=2.0*AEX*Q
      C(6)=Q
      ITM=6
      LMN(1)=IP(1)
      LMN(2)=IP(2)
      LMN(3)=IP(4)
      LMN(4)=IP(5)
      LMN(5)=IP(9)
      LMN(6)=IP(15)
      GO TO 100
26    CONTINUE
C     ITYP=14 XYY
      C(1)=AEY*AEY*AEX*Q
      C(2)=AEY*AEY*Q
      C(3)=2.0*AEY*AEX*Q
      C(4)=AEX*Q
      C(5)=2.0*AEY*Q
      C(6)=Q
      ITM=6
      LMN(1)=IP(1)
      LMN(2)=IP(2)
      LMN(3)=IP(3)
      LMN(4)=IP(6)
      LMN(5)=IP(8)
      LMN(6)=IP(16)
      GO TO 100
27    CONTINUE
C     ITYP=18 YYZ
      C(1)=AEY*AEY*AEZ*Q
      C(2)=2.0*AEY*AEZ*Q
      C(3)=AEY*AEY*Q
      C(4)= AEZ*Q
      C(5)=2.0*AEY*Q
      C(6)=Q
      ITM=6
      LMN(1)=IP(1)
      LMN(2)=IP(3)
      LMN(3)=IP(4)
      LMN(4)=IP(6)
      LMN(5)=IP(10)
      LMN(6)=IP(17)
      GO TO 100
28    CONTINUE
C     ITYP=16 XZZ
      C(1)=AEZ*AEZ*AEX*Q
      C(2)=AEZ*AEZ*Q
      C(3)=2.0*AEZ*AEX*Q
      C(4)=AEX*Q
      C(5)=2.0* AEZ*Q
      C(6)=Q
      ITM=6
      LMN(1)=IP(1)
      LMN(2)=IP(2)
      LMN(3)=IP(4)
      LMN(4)=IP(7)
      LMN(5)=IP(9)
      LMN(6)=IP(18)
      GO TO 100
29    CONTINUE
C     ITYP=19 YZZ
      C(1)=AEZ*AEZ*AEY*Q
      C(2)=AEZ*AEZ*Q
      C(3)=2.0*AEZ* AEY*Q
      C(4)= AEY*Q
      C(5)=2.0*AEZ*Q
      C(6)=Q
      ITM=6
      LMN(1)=IP(1)
      LMN(2)=IP(3)
      LMN(3)=IP(4)
      LMN(4)=IP(7)
      LMN(5)=IP(10)
      LMN(6)=IP(19)
      GO TO 100
30    CONTINUE
C     ITYP=15 XYZ
      C(1)=AEX*AEY*AEZ*Q
      C(2)=AEY*AEZ*Q
      C(3)=AEX*AEZ*Q
      C(4)=AEX*AEY*Q
      C(5)=AEZ*Q
      C(6)=AEY*Q
      C(7)=AEX*Q
      C(8)=Q
      ITM=8
      LMN(1)=IP(1)
      LMN(2)=IP(2)
      LMN(3)=IP(3)
      LMN(4)=IP(4)
      LMN(5)=IP(8)
      LMN(6)=IP(9)
      LMN(7)=IP(10)
      LMN(8)=IP(20)
100   CONTINUE
      RETURN
40    CONTINUE
C     GENERAL CODE
      L=LMNVAL(1,ITYP)+1
      M=LMNVAL(2,ITYP)+1
      N=LMNVAL(3,ITYP)+1
      X=1.0
      DO 150 L1=1,L
      XH(L1)=X
      X=X*AEX
150   CONTINUE
      Y=1.
      DO 151 M1=1,M
      YH(M1)=Y
      Y=Y*AEY
151   CONTINUE
      Z=1.
      DO 152 N1=1,N
      ZH(N1)=Z
      Z=Z*AEZ
152   CONTINUE
      ITM=0
      DO 130 L1=1,L
      X=FAC(L,L1)*XH(L1)*Q
      DO 131 M1=1,M
      Y=FAC(M,M1)*YH(M1)*X
      DO 132 N1=1,N
      Z=FAC(N,N1)*ZH(N1)*Y
      ITM=ITM+1
      LMN(ITM)=IBTOR(IBTSHL((L-L1),20),IBTOR(IBTSHL((M-M1),10),(N-N1)))
      C(ITM)=Z
132   CONTINUE
131   CONTINUE
130   CONTINUE
      RETURN
      ENTRY SETRHF
      FAC(1,1)=1.0
      FAC(2,1)=1.0
      FAC(2,2)=1.0
      DO 110 I=3,9
      FAC(I,1)=1.0
      FAC(I,I)=1.0
      JE=I-1
      DO 115 J=2,JE
      FAC(I,J)=FAC(I-1,J-1)+FAC(I-1,J)
115   CONTINUE
110   CONTINUE
      IP(1)=0
      IP(2)=IBTSHL(1,20)
      IP(3)=IBTSHL(1,10)
      IP(4)=1
      IP(5)=IBTSHL(2,20)
      IP(6)=IBTSHL(2,10)
      IP(7)=2
      IP(8)=IBTOR(IBTSHL(1,20),IBTSHL(1,10))
      IP(9)=IBTOR(IBTSHL(1,20),1)
      IP(10)=IBTOR(IBTSHL(1,10),1)
      IP(11)=IBTSHL(3,20)
      IP(12)=IBTSHL(3,10)
      IP(13)=3
      IP(14)=IBTOR(IBTSHL(2,20),IBTSHL(1,10))
      IP(15)=IBTOR(IBTSHL(2,20),1)
      IP(16)=IBTOR(IBTSHL(1,20),IBTSHL(2,10))
      IP(17)=IBTOR(IBTSHL(2,10),1)
      IP(18)=IBTOR(IBTSHL(1,20),2)
      IP(19)=IBTOR(IBTSHL(1,10),2)
      IP(20)=IBTOR(IBTSHL(1,20),IBTOR(IBTSHL(1,10),1))
      RETURN
      END
