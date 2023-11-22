      subroutine desm21(A,N1,N2,N3,N4,IJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N1,N2,N3,N4)
      DATA THIRD /.33333333333333D+0/,HALF/0.5D+0/
      IF(IJ.EQ.23)GOTO 23
      IF(IJ.EQ.14)GOTO 14
      IF(IJ.EQ.12)GOTO 12
      IF(IJ.EQ.13)GOTO 13
      IF(IJ.EQ.34)GOTO 34
 23   CONTINUE
      DO 10 I=1,N1
      DO 10 J=1,N2
      DO 10 K=1,J
      DO 10 L=1,N4
      A1=A(I,J,K,L)
      A2=A(I,K,J,L)
      X=A1+A2
      Y=(A1-A2)*THIRD
      A(I,J,K,L)=(X+Y)*HALF
      A(I,K,J,L)=(X-Y)*HALF
   10 CONTINUE
      GOTO 100
 14   CONTINUE
      DO 101 I=1,N1
      DO 101 J=1,N2
      DO 101 K=1,N3
      DO 101 L=1,I
      A1=A(I,J,K,L)
      A2=A(L,J,K,I)
      X=A1+A2
      Y=(A1-A2)*THIRD
      A(I,J,K,L)=(X+Y)*HALF
      A(L,J,K,I)=(X-Y)*HALF
 101  CONTINUE
      GOTO 100
 34   CONTINUE
      DO 102 I=1,N1
      DO 102 J=1,N2
      DO 102 K=1,N3
      DO 102 L=1,K
      A1=A(I,J,K,L)
      A2=A(I,J,L,K)
      X=A1+A2
      Y=(A1-A2)*THIRD
      A(I,J,K,L)=(X+Y)*HALF
      A(I,J,L,K)=(X-Y)*HALF
 102  CONTINUE
      GOTO 100
 12   CONTINUE
      DO 103 I=1,N1
      DO 103 J=1,I
      DO 103 K=1,N3
      DO 103 L=1,N4
      A1=A(I,J,K,L)
      A2=A(J,I,K,L)
      X=A1+A2
      Y=(A1-A2)*THIRD
      A(I,J,K,L)=(X+Y)*HALF
      A(J,I,K,L)=(X-Y)*HALF
 103  CONTINUE
      GOTO 100
 13   CONTINUE
      DO 113 I=1,N1
      DO 113 J=1,N2
      DO 113 K=1,I
      DO 113 L=1,N4
      A1=A(I,J,K,L)
      A2=A(K,J,I,L)
      X=A1+A2
      Y=(A1-A2)*THIRD
      A(I,J,K,L)=(X+Y)*HALF
      A(K,J,I,L)=(X-Y)*HALF
 113  CONTINUE
      GOTO 100
 100  CONTINUE
      RETURN
      END
      SUBROUTINE ienter(n)
      IMPLICIT REAL*8 (a-h,o-z)
      common/timeinfo/timein,timenow,timetot,timenew
      COMMON /etim/time(150),rn(150),ient(150),timex(200)
      if(n.eq.1)then
      CALL IZERO(IENT,150)
      CALL ZEROMA(TIME,1,150)
      CALL ZEROMA(TIMEX,1,200)
      endif
      call timer(1)
      ient(n) = ient(n) + 1
      time(n)=timenow
      RETURN
      END
      SUBROUTINE iexit(n)
      IMPLICIT REAL*8 (a-h,o-z)
      common/timeinfo/timein,timenow,timetot,timenew
      COMMON /etim/time(150),rn(150),ient(150),timex(200)
      call timer(1)
      time(n) = timenow - time(n)
      timex(n) = timex(n) + time(n)
      RETURN
      END
      integer function ign(i,j,k,no)
      implicit double precision (a-h,o-z)
      ign=(i-1)*no*no+(j-1)*no+k
      return
      end
      SUBROUTINE IPERM(K,L,M,N,IPRM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      call ienter(101)
      IF (L.GE.M.AND.M.GE.N)THEN
      IPRM=1
      ELSE
      IF (L.GE.N.AND.N.GE.M)THEN
      IPRM=2
      ELSE
      IF (K.GE.M.AND.L.GE.N)THEN
      IPRM=3
      ELSE
      IF (K.GE.N.AND.L.GE.M)THEN
      IPRM=4
      ELSE
      IF (K.GE.M.AND.M.GE.N)THEN
      IPRM=5
      ELSE
      IF (K.GE.N.AND.N.GE.M)THEN
      IPRM=6
      ELSE
      IF (M.GE.K.AND.L.GE.N)THEN
      IPRM=7
      ELSE
      IF (N.GE.K.AND.L.GE.M)THEN
      IPRM=8
      ELSE
      IF (M.GE.K.AND.K.GE.N)THEN
      IPRM=9
      ELSE
      IF (N.GE.K.AND.K.GE.M)THEN
      IPRM=10
      ELSE
      IF (M.GE.N.AND.N.GE.K)THEN
      IPRM=11
      ELSE
      IPRM=12
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      call iexit(101)
      RETURN
      END
      SUBROUTINE IPERM24(K,L,M,N,IPRM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      call ienter(102)
      IF (K.GT.L)THEN
      IF (L.GE.M.AND.M.GE.N)THEN
      IPRM=1
      ELSE
      IF (L.GE.N.AND.N.GE.M)THEN
      IPRM=2
      ELSE
      IF (K.GE.M.AND.L.GE.N)THEN
      IPRM=3
      ELSE
      IF (K.GE.N.AND.L.GE.M)THEN
      IPRM=4
      ELSE
      IF (K.GE.M.AND.M.GE.N)THEN
      IPRM=5
      ELSE
      IF (K.GE.N.AND.N.GE.M)THEN
      IPRM=6
      ELSE
      IF (M.GE.K.AND.L.GE.N)THEN
      IPRM=7
      ELSE
      IF (N.GE.K.AND.L.GE.M)THEN
      IPRM=8
      ELSE
      IF (M.GE.K.AND.K.GE.N)THEN
      IPRM=9
      ELSE
      IF (N.GE.K.AND.K.GE.M)THEN
      IPRM=10
      ELSE
      IF (M.GE.N.AND.N.GE.K)THEN
      IPRM=11
      ELSE
      IPRM=12
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ELSE
      IF (K.GE.M.AND.M.GE.N)THEN
      IPRM=13
      ELSE
      IF (K.GE.N.AND.N.GE.M)THEN
      IPRM=14
      ELSE
      IF (L.GE.M.AND.K.GE.N)THEN
      IPRM=15
      ELSE
      IF (L.GE.N.AND.K.GE.M)THEN
      IPRM=16
      ELSE
      IF (L.GE.M.AND.M.GE.N)THEN
      IPRM=17
      ELSE
      IF (L.GE.N.AND.N.GE.M)THEN
      IPRM=18
      ELSE
      IF (M.GE.L.AND.K.GE.N)THEN
      IPRM=19
      ELSE
      IF (N.GE.L.AND.K.GE.M)THEN
      IPRM=20
      ELSE
      IF (M.GE.L.AND.L.GE.N)THEN
      IPRM=21
      ELSE
      IF (N.GE.L.AND.L.GE.M)THEN
      IPRM=22
      ELSE
      IF (M.GE.N.AND.N.GE.L)THEN
      IPRM=23
      ELSE
      IPRM=24
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      RETURN
      call iexit(102)
      END
      SUBROUTINE IPERM4(J,K,L,M,IPRM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      call ienter(103)
      IF (L.GE.M)THEN
      IPRM=1
      ELSE
      IF (K.GE.M)THEN
      IPRM=2
      ELSE
      IF (J.GE.M)THEN
      IPRM=3
      ELSE
      IPRM=4
      ENDIF
      ENDIF
      ENDIF
      call iexit(103)
      RETURN
      END
      FUNCTION IT2(I,J)
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      if(j.le.i)then
         IT2=i*(i-1)/2+j
      else
         IT2=j*(j-1)/2+i
      endif
      RETURN
      END
      FUNCTION IT3(I,J,K)
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      I2=(I-1)*(I-2)
      IT3=I2+I2*(I-3)/6+J*(J-1)/2+K
      RETURN
      END
      FUNCTION IT4(I,J,K,L)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      I1=I-1
      J1=J-1
      K1=K-1
      I2=I1*I1
      IT=(I2*I2-13*I2+6*I1*(I2+1))/24
      J2=J1*J1
      JT=(J1*J2+3*J2-4*J1)/6
      KT=K1*(K1+1)/2
      IT4=IT+JT+KT+L
      RETURN
      END
      SUBROUTINE   SET   (I,J,K,L,MU)
C
      IF (J-K) 10,20,30
   10 IF (J-L) 11,14,17
   11 IF(K-L) 99,12,13
   12 MU=8
      RETURN
   13 MU=14
      RETURN
   14 IF (I-K) 99,15,16
   15 MU=2
      RETURN
   16 MU=11
      RETURN
   17 IF (I-K) 99,18,19
   18 MU=10
      RETURN
   19 MU=13
      RETURN
   20 IF (K-L) 99,21,24
   21 IF (I-J) 99,22,23
   22 MU=1
      RETURN
   23 MU=6
      RETURN
   24 IF (I-J) 99,25,26
   25 MU=4
      RETURN
   26 MU=9
      RETURN
   30 IF (I-J) 99,31,34
   31 IF (K-L) 99,32,33
   32 MU=3
      RETURN
   33 MU=5
      RETURN
   34 IF (K-L) 99,35,36
   35 MU=7
      RETURN
   36 MU=12
      RETURN
   99 MU=9999
      RETURN
      END
      SUBROUTINE SUMT1(T1HP,T1PH,NH,NP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A
      DIMENSION T1HP(NH,NP),T1PH(NP,NH)
      DO 10 I=1,NH
      DO 10 A=1,NP
      T1HP(I,A)=T1HP(I,A)+T1PH(A,I)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE SYMETR(T2HPPH,NH,NP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      DIMENSION T2HPPH(NH,NP,NP,NH)
      DO 10 I=1,NH
      DO 10 J=1,I
      DO 10 A=1,NP
      LIMB=NP
      IF(I.EQ.J)LIMB=A
      DO 10 B=1,LIMB
      X=T2HPPH(I,A,B,J)
      Y=T2HPPH(J,B,A,I)
      T2HPPH(I,A,B,J)=X+Y
      T2HPPH(J,B,A,I)=X+Y
   10 CONTINUE
      RETURN
      END  

      
      SUBROUTINE SYMT21(A,N1,N2,N3,N4,IJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N1,N2,N3,N4)
      DATA TWO/2.0D+0/
      IF(IJ.EQ.23)GOTO 23
      IF(IJ.EQ.14)GOTO 14
      IF(IJ.EQ.12)GOTO 12
      IF(IJ.EQ.34)GOTO 34
      IF(IJ.EQ.13)GOTO 13
 23   CONTINUE
      DO 10 I=1,N1
      DO 10 J=1,N2
      DO 10 K=1,J
      DO 10 L=1,N4
      X=A(I,J,K,L)
      Y=A(I,K,J,L)
      A(I,J,K,L)=TWO*X-Y
      A(I,K,J,L)=TWO*Y-X
   10 CONTINUE
      GOTO 1000
 14   CONTINUE
      DO 11 I=1,N1
      DO 11 J=1,N2
      DO 11 K=1,N3
      DO 11 L=1,I
      X=A(I,J,K,L)
      Y=A(L,J,K,I)
      A(I,J,K,L)=TWO*X-Y
      A(L,J,K,I)=TWO*Y-X
 11   CONTINUE
      GOTO 1000
 12   CONTINUE
      DO 112 I=1,N1
      DO 112 J=1,I
      DO 112 K=1,N3
      DO 112 L=1,N4
      X=A(I,J,K,L)
      Y=A(J,I,K,L)
      A(I,J,K,L)=TWO*X-Y
      A(J,I,K,L)=TWO*Y-X
 112  CONTINUE
      GOTO 1000
 13   CONTINUE
      DO 113 I=1,N1
      DO 113 J=1,N2
      DO 113 K=1,I
      DO 113 L=1,N4
      X=A(I,J,K,L)
      Y=A(K,J,I,L)
      A(I,J,K,L)=TWO*X-Y
      A(K,J,I,L)=TWO*Y-X
 113  CONTINUE
      GOTO 1000
 34   CONTINUE
      DO 134 I=1,N1
      DO 134 J=1,N2
      DO 134 K=1,N3
      DO 134 L=1,K
      X=A(I,J,K,L)
      Y=A(I,J,L,K)
      A(I,J,K,L)=TWO*X-Y
      A(I,J,L,K)=TWO*Y-X
 134  CONTINUE
 1000 CONTINUE
      RETURN
      END
      SUBROUTINE SYMV(V,NU)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      DIMENSION V(NU,NU,NU,NU)
      data half/0.5d0/,two/2.0d+0/
      nu4=nu*Nu*nu*nu
 10   continue
      DO 200 I=1,NU
      DO 200 J=1,I
      DO 200 A=1,NU
      LIMB=NU
      IF(I.EQ.J)LIMB=A
      DO 200 B=1,LIMB
      Z=V(a,i,j,b)+v(b,j,i,a)
      V(a,i,j,b)=z
      V(b,j,i,a)=z
 200  CONTINUE
      do 20 i=1,nu
      do 20 j=1,nu
 20   continue
 123  format('v:',3f15.10)
      RETURN
      END
      subroutine symvt(no,nu,vtij,vtji)
      implicit double precision (a-h,o-z)
      COMMON/NVTAP/NTITER,NQ2,NVT
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      integer a,b
      dimension vtij(nu,nu,no,no),vtji(nu,nu,no,no)
      do 11 i=1,no
      do 11 j=1,i
      do 12 k=1,no
      kk=no2*(i-1)+no*(j-1)+k
      call rdgen(nvt,kk,nou2,vtij(1,1,1,k))
      kk=no2*(j-1)+no*(i-1)+k
      call rdgen(nvt,kk,nou2,vtji(1,1,1,k))
 12   continue
      do 10 k=1,no
      do 10 l=1,no
      do 10 a=1,nu
      limb=nu
      if(i.eq.j)limb=a
      do 10 b=1,limb
      x=vtij(b,a,l,k)
      y=vtji(a,b,l,k)
      z=x+y
 777  format('i,j,k,l,a,b,zyn:',6i3,3f15.10)
      vtij(b,a,l,k)=z
      if (i.eq.j)then
      vtij(a,b,l,k)=z
      else
      vtji(a,b,l,k)=z
      endif
 10   continue
      do 13 k=1,no
      kk=no2*(j-1)+no*(i-1)+k
      call wrgen(nvt,kk,nou2,vtij(1,1,1,k))
      kk=no2*(i-1)+no*(j-1)+k
      if(i.ne.j)call wrgen(nvt,kk,nou2,vtji(1,1,1,k))
 13   continue
 11   continue
 788  format('i,j,k:',3i3)
 789  format('vtnew:',3f15.10)
      return
      end
      DOUBLE PRECISION FUNCTION timt(elapse)
      IMPLICIT REAL*8 (a-h,o-z)
C=================================================================
C=                                                               =
C=    Time Routine for the RS/6000 based on the MCLOCK intrinsic.=
C=                                                               =
C=    Note:  The routine MCLOCK returns an I*4 which gives       =
C=           process time in hundreths of a second.              =
C=                                                               =
C=================================================================
      LOGICAL*1 first
      DATA first/.TRUE./

      IF (first) THEN
         base = dfloat(mclock())/100.0D0
         first = .FALSE.
      END IF
      timt = dfloat(mclock())/100.0D0 - base
      RETURN
      END
      SUBROUTINE TRANSQ(A,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N,N)
      DO 10 I=1,N
      DO 10 J=1,I
      X=A(I,J)
      A(I,J)=A(J,I)
      A(J,I)=X
   10 CONTINUE
      RETURN
      END
      SUBROUTINE VECADD(A,B,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*),B(*)
      call ienter(76)
      DO 10 I=1,N
      A(I)=A(I)+B(I)
   10 CONTINUE
      call iexit(76)
      RETURN
      END  
      SUBROUTINE VECCOP(N,A,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N),B(N)
      call ienter(77)
      call dcopy(n,b,1,a,1)
      call iexit(77)
      RETURN
      END
      SUBROUTINE VECMUL(A,N,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N)
      call ienter(78)
      call dscal(n,x,a,1)
      call iexit(78)
      RETURN
      END
      SUBROUTINE VECSUB(A,B,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*),B(*)
      call ienter(79)
      DO 10 I=1,N
      A(I)=A(I)-B(I)
   10 CONTINUE
      call iexit(79)
      RETURN
      END  
      SUBROUTINE VMINUS(V,LEN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(LEN)
      DATA X /-1.D0/
      call ienter(82)
      DO 10 I=1,LEN
       V(I)=-V(I)
10    CONTINUE
      call iexit(82)
      RETURN
      END
      SUBROUTINE ZEROMA(WORK,NF,NL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WORK(NL)
      DO 1 I=NF,NL
1     WORK(I)=0.0d0
      RETURN
      END
      SUBROUTINE symt3four(v,t3,nu)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION t3(nu,nu,nu),v(nu,nu,nu)
      data two/2.0d+0/
      do 1 j=1,nu
      do 1 k=1,nu
      do 1 i=1,nu
         v(i,j,k)=(t3(i,j,k)-t3(j,i,k))*two-t3(i,k,j)+t3(k,i,j)
 1    continue
      RETURN
      END
      subroutine icopymm(ia,ib,n)
      implicit double precision (a-h,o-z)
      dimension ia(n),ib(n)
      do 10 i=1,n
      ib(i)=ia(i)
 10   continue
      return
      end
