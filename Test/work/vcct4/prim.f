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
      integer function ign(i,j,k,no)
      implicit double precision (a-h,o-z)
      ign=(i-1)*no*no+(j-1)*no+k
      return
      end
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
c      IT3=I2+I2*(I-3)/6+J*(J-1)/2+K+i-1
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
c      write (6,*)'ijkl,itjtkt:',i,j,k,l,it,jt,kt
      IT4=IT+JT+KT+L
c      write(6,11)i,j,k,l,it4
c 11   format('i,j,k,l,it4:',5i3)
      RETURN
      END
      SUBROUTINE   SET   (I,J,K,L,MU)
C
      IF (J-K) 10,20,30
C..........IJIJ IJKK IJIL IJKJ IJKL IJKL
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
C..........IIII IIIL IJJJ IJJL
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
C..........IIKL IIKK IJKK IJKL
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
      DIMENSION A(N),B(N)
      DO 10 I=1,N
      A(I)=A(I)+B(I)
   10 CONTINUE
      RETURN
      END  
      SUBROUTINE VECCOP(N,A,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N),B(N)
      call dcopy(n,b,1,a,1)
c      DO 10 I=1,N
c      A(I)=B(I)
c 10   CONTINUE
      RETURN
      END
      SUBROUTINE IVECCOP(N,iA,iB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IA(N),IB(N)
      DO 10 I=1,N
      IA(I)=IB(I)
 10   CONTINUE
      RETURN
      END
      SUBROUTINE VECMUL(A,N,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N)
      DO 10 I=1,N
      A(I)=X*A(I)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE VECSUB(A,B,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N),B(N)
      DO 10 I=1,N
      A(I)=A(I)-B(I)
   10 CONTINUE
      RETURN
      END  
      SUBROUTINE VECSUBA(A,B,N,alpha)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N),B(N)
      DO 10 I=1,N
      A(I)=alpha*A(I)-B(I)
   10 CONTINUE
      RETURN
      END  
      SUBROUTINE ZEROMA(WORK,NF,NL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WORK(NL)
      DO 1 I=NF,NL
1     WORK(I)=0.0D0
      RETURN
      END
      SUBROUTINE MATMUL(A,B,C,NI,NJ,NK,IACC,IFUN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision mone
      DIMENSION A(NI,NK),B(NK,NJ),C(NI,NJ)
      DATA ONE,ZERO,MONE/1.0D+0,0.0D+0,-1.0D+0/
C
C
C     THE ARGUMENT TOL IS NOT USED. IN THE VAX VERSION IT
C     IS A TOLERANCE FOR SAXPYS.
C
C     DEPENDING ON THE VALUES OF IACC AND IFUN, FOUR SORTS
C     OF MATRIX MULTIPLICATIONS MAY BE CARRIED OUT :
C
C     IACC=  0, IFUN= 0    C = C + A*B
C     IACC=  0, IFUN><0    C = C - A*B
C     IACC=><0, IFUN= 0    C =     A*B
C     IACC=><0, IFUN><0    C =   - A*B
C     
      IF (IACC.EQ.0.AND.IFUN.EQ.0) THEN
      ALPHA=ONE
      BETA=ONE
      ELSE
      IF (IACC.EQ.0.AND.IFUN.NE.0) THEN
      ALPHA=MONE
      BETA=ONE
      ELSE
      IF (IACC.NE.0.AND.IFUN.EQ.0) THEN
      ALPHA=ONE
      BETA=ZERO
      ELSE
      ALPHA=MONE
      BETA=ZERO
      ENDIF
      ENDIF
      ENDIF
c      CALL SGEMM('N','N',NI,NJ,NK,ALPHA,A,NI,B,NK,BETA,C,NI)
      CALL XGEMM('N','N',NI,NJ,NK,ALPHA,A,NI,B,NK,BETA,C,NI)
C FPS 164
C      CALL PMMUL(A,B,C,NI,NK,NI,NI,NJ,NK,IACC,IFUN,8192,IERR)
      RETURN
      END
      SUBROUTINE MATMULSK(A,B,C,NI,NJ,NK,IACC,IFUN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision mone
      DIMENSION A(NI,NK),B(NK,NJ),C(NI,NJ)
      DATA ONE,ZERO,MONE/1.0D+0,0.0D+0,-1.0D+0/
C
C
C     THE ARGUMENT TOL IS NOT USED. IN THE VAX VERSION IT
C     IS A TOLERANCE FOR SAXPYS.
C
C     DEPENDING ON THE VALUES OF IACC AND IFUN, FOUR SORTS
C     OF MATRIX MULTIPLICATIONS MAY BE CARRIED OUT :
C
C     IACC=  0, IFUN= 0    C = C + A*B
C     IACC=  0, IFUN><0    C = C - A*B
C     IACC=><0, IFUN= 0    C =     A*B
C     IACC=><0, IFUN><0    C =   - A*B
C     
      IF (IACC.EQ.0.AND.IFUN.EQ.0) THEN
      ALPHA=ONE
      BETA=ONE
      ELSE
      IF (IACC.EQ.0.AND.IFUN.NE.0) THEN
      ALPHA=MONE
      BETA=ONE
      ELSE
      IF (IACC.NE.0.AND.IFUN.EQ.0) THEN
      ALPHA=ONE
      BETA=ZERO
      ELSE
      ALPHA=MONE
      BETA=ZERO
      ENDIF
      ENDIF
      ENDIF
c      CALL SGEMM('N','N',NI,NJ,NK,ALPHA,A,NI,B,NK,BETA,C,NI)
      CALL XGEMM('N','N',NI,NJ,NK,ALPHA,A,NI,B,NK,BETA,C,NI)
C FPS 164
C      CALL PMMUL(A,B,C,NI,NK,NI,NI,NJ,NK,IACC,IFUN,8192,IERR)
      RETURN
      END
      SUBROUTINE TRT1(NO,NU,TI,T1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A
      DIMENSION TI(NO,NU),T1(NU,NO)
      NOU=NO*NU
      CALL VECCOP(NOU,TI,T1)
      DO 10 I=1,NO
      DO 10 A=1,NU
      T1(A,I)=TI(I,A)
 10   CONTINUE
      RETURN
      END
      SUBROUTINE RINTHH(IP,N,TI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/newio/npp1,nhh,npp,ntvoe,nto3,no1,nu1
      INTEGER A,B
      DIMENSION TI(N,N)
c      NHH=58
      IF(IP.EQ.1)THEN
         n2=n*n
         IR=1
         call rdd(nhh,ir,n2,ti)
c         READ(NHH,REC=IR)TI
      ELSE
         DO 10 I=1,N
            IR=I+1
         call rdd(nhh,ir,n,ti(1,i))
 10      CONTINUE
      ENDIF
      RETURN
      END
      SUBROUTINE RINTPP(IP,N,TI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/newio/npp1,nhh,npp,ntvoe,nto3,no1,nu1
      INTEGER A,B
      DIMENSION TI(N,N)
c      NPP=59
      n2=n*n
      IF(IP.EQ.1)THEN
         IR=1
         call rdd(npp,ir,n2,ti)
      ELSE
         DO 10 I=1,N
            IR=I+1
         call rdd(npp,ir,n,ti(1,i))
 10      CONTINUE
      ENDIF
      RETURN
      END
      SUBROUTINE RINTPPI(Nu,N,TI,ia)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/newio/npp1,nhh,npp,ntvoe,nto3,no1,nu1
      INTEGER A,B
      DIMENSION TI(N,Nu)
      DO 10 I=1,Nu
         IR=(ia-1)*nu+I+1
         call rdd(npp,ir,n,ti(1,i))
 10   CONTINUE
      RETURN
      END
      SUBROUTINE RINTVO(IP,NO,NU,TI,T2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      common/newio/npp1,nhh,npp,ntvoe,nto3,no1,nu1
      DIMENSION TI(NU,NU,NO),T2(NO,NU,NU,NO)
      data zero/0.0d+0/
c      NTVOE=60
      IVO=IABS(IP)
      NLAST=NO*(IVO-1)
      DO 1 I=1,NO
      IASV=NLAST+I
      READ(NTVOE,REC=IASV)TI
      if(ip.eq.1)then 
         endif
      DO 1 J=1,NO
      DO 1 A=1,NU
      DO 1 B=1,NU
         X=ZERO
         IF (IP.LT.0)X=T2(I,A,B,J)
      T2(I,A,B,J)=X+TI(A,B,J)
 1    CONTINUE
      RETURN
      END
      SUBROUTINE RDVT3ONW(I,J,K,NU,V)
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      COMMON/newio/ni1,ni2,ni3,ni4,no3,no1,nu1
      DIMENSION V(NU,NU,NU)
c      write(6,*)'uwaga************************** rdvt3onw start'
      nu3=nu*nu*nu
 489  format('i,j,k:',i3,5x,3i6,3x,3i3)
      IF (I.GT.J)THEN
      IF(J.GT.K) THEN
         IAS=IT3(I,J,K)
         call rpakt3(ias,nu3,v)
c      READ(NO3,REC=IAS)V
      ELSE
      IF(I.GT.K)THEN
      IAS=IT3(I,K,J)
      call rpakt3(ias,nu3,v)
c      READ(NO3,REC=IAS)V
      CALL TRANMD(V,NU,NU,NU,1,23)
      ELSE
      IAS=IT3(K,I,J)
      call rpakt3(ias,nu3,v)
c      READ(NO3,REC=IAS)V
      CALL TRANMD(V,NU,NU,NU,1,312)
      ENDIF
      ENDIF
      ELSE
      IF (I.GT.K)THEN
      IAS=IT3(J,I,K)
      call rpakt3(ias,nu3,v)
c      READ(NO3,REC=IAS)V
      CALL TRANMD(V,NU,NU,NU,1,12)
      ELSE
      IF(J.GT.K) THEN
      IAS=IT3(J,K,I)
      call rpakt3(ias,nu3,v)
c      READ(NO3,REC=IAS)V
      CALL TRANMD(V,NU,NU,NU,1,231)
      ELSE
      IAS=IT3(K,J,I)
      call rpakt3(ias,nu3,v)
c      READ(NO3,REC=IAS)V
      CALL TRANMD(V,NU,NU,NU,1,13)
      ENDIF
      ENDIF
      ENDIF
c      write(6,*)'uwaga************************** rdvt3onw end'
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
      SUBROUTINE VRINTHH(IP,N,TI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/newio/npp1,nhh,npp,ntvoe,nto3,no1,nu1
      INTEGER A,B
      DIMENSION TI(N,N)
c      NPP=58
      IF(IP.EQ.1)THEN
         n2=n*n
         IR=1
            call wrt(nhh,ir,n2,ti)
      ELSE
         DO 10 I=1,N
            IR=I+1
            call wrt(nhh,ir,n,ti(1,i))
 10      CONTINUE
      ENDIF
      RETURN
      END
      SUBROUTINE VRINTPP(IP,N,TI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/newio/npp1,nhh,npp,ntvoe,nto3,no1,nu1
      INTEGER A,B
      DIMENSION TI(N,N)
c      NPP=59
      n2=n*n
      IF(IP.EQ.1)THEN
         IR=1
            call wrt(npp,ir,n2,ti)
         WRITE(NPP,REC=IR)TI
      ELSE
         DO 10 I=1,N
            IR=I+1
            call wrt(npp,ir,n,ti(1,i))
c            WRITE(NPP,REC=IR)TI(1,I)
 10      CONTINUE
      ENDIF
      RETURN
      END
      SUBROUTINE VRINTVO(IP,NO,NU,TI,T2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      common/newio/npp1,nhh,npp,ntvoe,nto3,no1,nu1
      DIMENSION TI(NU,NU,NO),T2(NO,NU,NU,NO)
c      NTVOE=60
      NLAST=NO*(IP-1)
c      write(6,*)'vrintvo,no,nu,nlast,ip',no,nu,nlast,ip
      DO 1 I=1,NO
      IASV=NLAST+I
      DO 2 J=1,NO
      DO 2 A=1,NU
      DO 2 B=1,NU
      TI(A,B,J)=T2(I,A,B,J)
 2    CONTINUE
      WRITE(NTVOE,REC=IASV)TI
 1    CONTINUE
      RETURN
      END
      subroutine rdd(npp,ir,n,t)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension t(n)
      read(npp,rec=ir)t
      return
      end
      subroutine wrt(npp,ir,n,t)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension t(n)
      write(npp,rec=ir)t
      return
      end
      
      SUBROUTINE SYMT3(V,NU,ID)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      DIMENSION V(NU,NU,NU)
      DATA TWO/2.0D+0/
      GO TO (1,2,3,4,5,6),ID
    1 CONTINUE
      DO 100 A=1,NU
      DO 100 B=1,A
      DO 100 C=1,B
      XAB=V(A,B,C)
      XAC=V(A,C,B)
      XBA=V(B,A,C)
      XBC=V(B,C,A)
      XCA=V(C,A,B)
      XCB=V(C,B,A)
      V(A,B,C)=TWO*XAB-XAC-XBA
      V(A,C,B)=TWO*XAC-XAB-XCA
      V(B,A,C)=TWO*XBA-XBC-XAB
      V(B,C,A)=TWO*XBC-XBA-XCB
      V(C,A,B)=TWO*XCA-XCB-XAC
      V(C,B,A)=TWO*XCB-XCA-XBC
  100 CONTINUE
      GO TO 1000
 3    CONTINUE
      DO 300 A=1,NU
      DO 300 B=1,A
      DO 300 C=1,B
      XAB=V(A,B,C)
      XAC=V(A,C,B)
      XBA=V(B,A,C)
      XBC=V(B,C,A)
      XCA=V(C,A,B)
      XCB=V(C,B,A)
      V(A,B,C)=TWO*XAB-XCB-XBA
      V(A,C,B)=TWO*XAC-XBC-XCA
      V(B,A,C)=TWO*XBA-XCA-XAB
      V(B,C,A)=TWO*XBC-XAC-XCB
      V(C,A,B)=TWO*XCA-XBA-XAC
      V(C,B,A)=TWO*XCB-XAB-XBC
  300 CONTINUE
      GO TO 1000
 6    CONTINUE
      DO 101 A=1,NU
      DO 101 B=1,A
      DO 101 C=1,B
      XAB=V(A,B,C)
      XAC=V(A,C,B)
      XBA=V(B,A,C)
      XBC=V(B,C,A)
      XCA=V(C,A,B)
      XCB=V(C,B,A)
      V(A,B,C)=TWO*XAB-XAC-XCB
      V(A,C,B)=TWO*XAC-XAB-XBC
      V(B,A,C)=TWO*XBA-XBC-XCA
      V(B,C,A)=TWO*XBC-XBA-XAC
      V(C,A,B)=TWO*XCA-XCB-XBA
      V(C,B,A)=TWO*XCB-XCA-XAB
 101  CONTINUE
      GOTO 1000
    2 CONTINUE
      DO 200 A=1,NU
      DO 200 B=1,A
      DO 200 C=1,NU
      X=V(A,B,C)
      Y=V(B,A,C)
      V(A,B,C)=TWO*X-Y
      V(B,A,C)=TWO*Y-X
  200 CONTINUE
      GO TO 1000
    4 CONTINUE
      DO 400 B=1,NU
      DO 400 C=1,B
      DO 400 A=1,C
      X=V(A,B,C)
      V(A,B,C)=V(B,C,A)
      V(B,C,A)=V(C,A,B)
      V(C,A,B)=X
      IF(B.EQ.C.OR.C.EQ.A)GO TO 400
      X=V(B,A,C)
      V(B,A,C)=V(A,C,B)
      V(A,C,B)=V(C,B,A)
      V(C,B,A)=X
  400 CONTINUE
      GO TO 1000
    5 CONTINUE
      DO 500 A=1,NU
      DO 500 C=1,A
      DO 500 D=1,C
      X=V(C,D,A)
      V(C,D,A)=V(A,C,D)
      V(A,C,D)=V(D,A,C)
      V(D,A,C)=X
      IF (A.EQ.C.OR.C.EQ.D)GO TO 500
      X=V(A,D,C)
      V(A,D,C)=V(C,A,D)
      V(C,A,D)=V(D,C,A)
      V(D,C,A)=X
  500 CONTINUE
      GO TO 1000
 1000 CONTINUE
      RETURN
      END
      SUBROUTINE T2DEN(T2HPPH,TI,EH,EP,NH,NP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      COMMON/ITRAT/ITER,MAXIT,ICCNV
      DIMENSION T2HPPH(NH,NP,NP,NH),EH(NH),EP(NP),TI(NP,NP,NH)
      DATA ZERO/0.0D+0/,TRES/1.0D-10/,ONE/1.0D+0/
      DO 10 I=1,NH
      DO 10 A=1,NP
      DO 10 B=1,NP
      DO 10 J=1,NH
      DEN=EH(I)+EH(J)-EP(A)-EP(B)
      XT2=T2HPPH(I,A,B,J)
      T2HPPH(I,A,B,J)=XT2/DEN
   10 CONTINUE
      RETURN
      END
