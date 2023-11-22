      SUBROUTINE WRGEN(NTAP,IREC,N,A)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N)
      call ienter(89)
      WRITE(NTAP,REC=IREC)A
      call iexit(89)
      RETURN
      END
      SUBROUTINE WROV4(INO,NO,NU,TI,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION TI(NU,NU,NU),V(NU,NU,NU,NU)
      NNO=NO
      NNU=NU
      IF(INO.EQ.1)THEN
      NNO=NU
      NNU=NO
      ENDIF
      INO1=1-INO
      NLAST=6*NNO+3*NNU+INO1*NNO+7
      DO 1 A=1,NU
      DO 2 B=1,NU
      DO 2 C=1,NU
      DO 2 D=1,NU
      TI(B,C,D)= V(A,B,C,D)
 2    CONTINUE
      IASV=NLAST+A
      WRITE(NALL4,REC=IASV)TI
 1    CONTINUE
      RETURN
      END
      SUBROUTINE WROV4A(INO,NO,NU,TI,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      DIMENSION TI(NU,NU,NU),V(NU,NU,NU,NU)
      NNO=NO
      NNU=NU
      IF(INO.EQ.1)THEN
      NNO=NU
      NNU=NO
      ENDIF
      INO1=1-INO
      NLAST=INO1*NNO
      NALL4=45
      DO 1 A=1,NU
      DO 2 B=1,NU
      DO 2 C=1,NU
      DO 2 D=1,NU
      TI(B,C,D)= V(A,B,C,D)
 2    CONTINUE
      IASV=NLAST+A
      WRITE(NALL4,REC=IASV)TI
 1    CONTINUE
      RETURN
      END

      subroutine wro4(iasv,nu,v)
      implicit double precision(a-h,o-z)
      common/newt4/nt4,no4,lt4,nall4,ll4
      dimension v(nu,nu,nu,nu)
      write(no4,rec=iasv)v
      return
      end
      SUBROUTINE WRT4(IASV,NU,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NEWT4/NT4,NO4,LT4,NALL4,LL4
      DIMENSION V(NU,NU,NU,NU)
      call ienter(90)
      WRITE(NT4,REC=IASV)V
      call iexit(90)
      RETURN
      END
      SUBROUTINE WRVA(IAS,NU,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A
      DIMENSION V(NU,NU,NU,NU)
      call ienter(91)
      DO 10 A=1,NU
      IASV=IAS+A
      CALL WRVT3(IASV,NU,V(1,1,1,A))
   10 CONTINUE
      call iexit(91)
      RETURN
      END
      SUBROUTINE WRVEM(IS,NO,NU,TI,VE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION TI(NU,NU,NU),VE(NU,NU,NU,NO)
      call ienter(92)
      NNO=NO
      NNU=NU
      IF(IS.EQ.11.OR.IS.EQ.21)THEN
      NNO=NU
      NNU=NO
      ENDIF
      IF (IS.GE.20) THEN
      IS1=IS-20
      NLAST=9*NNO+4*NNU+7+IS1*NU
      ELSE
      IS1=IS-10
      NLAST=5*NNO+2*NNU+7+IS1*NU
      ENDIF
      DO 20 I=1,NO
      DO 21 A=1,NU
      DO 21 B=1,NU
      DO 21 C=1,NU
      TI(A,B,C)=VE(A,B,C,I)
 21   CONTINUE
      IASV=NLAST+I
      WRITE(NALL4,REC=IASV)TI
 20   CONTINUE
      call iexit(92)
      RETURN
      END
      SUBROUTINE WRVEM1Q2(IN,NO,NU,TI,VEM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C
      DIMENSION TI(NU,NU,NU),VEM(NO,NU,NU,NU)
      NALL4=45
      nu3=nu*nu*Nu
      nou3=no*nu3
      NNO=NO
      NNU=NU
      IF (IN.EQ.1.OR.IN.EQ.11)THEN
      NNO=NU
      NNU=NO
      ENDIF
      I11=0
      IF(IN.EQ.10.OR.IN.EQ.11)I11=NO+NU
      NLAST=3*NNO+NNU+I11
      INO=1
      IF(IN.EQ.1.OR.IN.EQ.11)INO=0
      DO 2 I=1,NO
      DO 1 A=1,NU
      DO 1 B=1,NU
      DO 1 C=1,NU
      TI(A,B,C)=VEM(I,C,A,B)
 1    CONTINUE
      IASV=NLAST+NU*INO+I
      WRITE(NALL4,REC=IASV)TI
 2    CONTINUE
      RETURN
      END
      SUBROUTINE WRVEMQ2(IN,NO,NU,TI,VEM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C
      DIMENSION TI(NU,NU,NU),VEM(NU,NU,NU,NO)
      NALL4=45
      nu3=nu*nu*Nu
      nou3=no*nu3
      NNO=NO
      NNU=NU
      IF (IN.EQ.1.OR.IN.EQ.11)THEN
      NNO=NU
      NNU=NO
      ENDIF
      I11=0
      IF(IN.EQ.10.OR.IN.EQ.11)I11=NO+NU
      NLAST=3*NNO+NNU+I11
      INO=1
      IF(IN.EQ.1.OR.IN.EQ.11)INO=0
      DO 2 I=1,NO
      DO 1 A=1,NU
      DO 1 B=1,NU
      DO 1 C=1,NU
      TI(A,B,C)=VEM(A,B,C,I)
 1    CONTINUE

      IASV=NLAST+NU*INO+I
      WRITE(NALL4,REC=IASV)TI
 2    CONTINUE
      RETURN
      END
      SUBROUTINE WRVO(NO,NU,TI,VOE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION TI(NU,NU,NO),VOE(NO,NU,NU,NO)
      call ienter(93)
      NLAST=8*NO+4*NU+7
      nou2=no*nu*nu
      no2u2=no*nou2
      DO 2 I=1,NO
      DO 1 A=1,NU
      DO 1 B=1,NU
      DO 1 J=1,NO
      TI(A,B,J)=VOE(I,A,B,J)
 1    CONTINUE
      IASV=NLAST+I
      WRITE(NALL4,REC=IASV)TI
 2    CONTINUE
      call iexit(93)
      RETURN
      END
      SUBROUTINE WRVOA(NO,NU,TI,VOE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      DIMENSION TI(NU,NU,NO),VOE(NO,NU,NU,NO)
      NALL4=45
      NLAST=2*NO+NU
      DO 2 I=1,NO
      DO 1 A=1,NU
      DO 1 B=1,NU
      DO 1 J=1,NO
      TI(A,B,J)=VOE(I,A,B,J)
 1    CONTINUE
      IASV=NLAST+I
      WRITE(NALL4,REC=IASV)TI
 2    CONTINUE
      RETURN
      END
      SUBROUTINE WRVOE(NO,NU,TI,VOE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION TI(NU,NU,NO),VOE(NO,NU,NU,NO)
      call ienter(94)
      NLAST=7*NO+4*NU+7
      nou2=no*nu*nu
      no2u2=no*nou2
      DO 2 I=1,NO
      DO 1 A=1,NU
      DO 1 B=1,NU
      DO 1 J=1,NO
      TI(A,B,J)=VOE(I,A,B,J)
 1    CONTINUE
      IASV=NLAST+I
      WRITE(NALL4,REC=IASV)TI
 2    CONTINUE
      call iexit(94)
      RETURN
      END
      SUBROUTINE WRVOEA(NO,NU,TI,VOE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      DIMENSION TI(NU,NU,NO),VOE(NO,NU,NU,NO)
      NALL4=45
      NLAST=NO+NU
      DO 2 I=1,NO
      DO 1 A=1,NU
      DO 1 B=1,NU
      DO 1 J=1,NO
      TI(A,B,J)=VOE(I,A,B,J)
 1    CONTINUE
      IASV=NLAST+I
      WRITE(NALL4,REC=IASV)TI
 2    CONTINUE
      RETURN
      END
      SUBROUTINE WRVOEP(IS,NO,NU,TI,VOE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION TI(NU,NU,NO),VOE(NO,NU,NU,NO)
      call ienter(95)
      NLAST=10*NO+5*NU+7
      NLAST=NLAST+(IS-1)*NO
      nou2=no*nu*nu
      no2u2=no*nou2
      DO 2 I=1,NO
      DO 1 A=1,NU
      DO 1 B=1,NU
      DO 1 J=1,NO
      TI(A,B,J)=VOE(I,A,B,J)
 1    CONTINUE
      IASV=NLAST+I
      WRITE(NALL4,REC=IASV)TI
 2    CONTINUE
      call iexit(95)
      RETURN
      END
      SUBROUTINE WRVT3(IASV,NU,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NEWTAPE/NT2T4,NT4INT
      DIMENSION V(NU,NU,NU)
      call ienter(100)
      WRITE(NT4INT,REC=IASV)V
      call iexit(100)
      RETURN
      END
      SUBROUTINE WRVT3I(IASV,NU,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NEWT3/NT3INT,LT3INT
      DIMENSION V(NU,NU,NU)
      call ienter(96)
      WRITE(NT3INT,REC=IASV)V
      call iexit(96)
      RETURN
      END
      SUBROUTINE WRVT3ID(IASV,NU,V,VOLD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NEWT3/NT3INT,LT3INT
      DIMENSION V(NU,NU,NU),VOLD(NU,NU,NU)
      NU3=NU*NU*NU
      READ(NT3INT,REC=IASV)VOLD
      CALL VECADD(V,VOLD,NU3)
      WRITE(NT3INT,REC=IASV)V
      RETURN
      END
      SUBROUTINE WRVT3IMN(NO,NU,TI,VT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A
      COMMON/NEWT3/NT3INT,LT3INT
      DIMENSION TI(NO,NO,NO),VT(NU,NO,NO,NO,NO,NO)
      call ienter(97)
      NO2=NO*NO
      NO3=NO2*NO
      DO 10 A=1,NU
      DO 10 M=1,NO
      DO 10 N=1,NO
      DO 11 I=1,NO
      DO 11 J=1,NO
      DO 11 K=1,NO
      TI(I,J,K)=VT(A,M,N,I,J,K)
 11   CONTINUE
      IAS=2*NO3+NO2*(A-1)+NO*(M-1)+N
      WRITE(NT3INT,REC=IAS)TI
 10   CONTINUE
      call iexit(97)
      RETURN
      END
      SUBROUTINE WRVT3N(IASV,NU,V)
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      COMMON/NT3T3/NO3,NT3,LT3
      DIMENSION V(NU,NU,NU)
      WRITE(NT3,REC=IASV)V
      RETURN
      END
      SUBROUTINE WRVT3O(IASV,NU,V)
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      COMMON/NT3T3/NO3,NT3,LT3
      DIMENSION V(NU,NU,NU)
      WRITE(NO3,REC=IASV)V
      RETURN
      END
