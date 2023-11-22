      SUBROUTINE CFCVEC(N,F,C,T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION F(N,N),C(N,N),T(N,N)
      CALL MATMUL(F,C,T,N,N,N,1,0)
      CALL TRANSQ(C,N)
      CALL MATMUL(C,T,F,N,N,N,1,0)
      RETURN
      END
      subroutine copt3(i,j,k,no,nu,ti,t3)
      implicit double precision (a-h,o-z)
      integer a,b,c
      dimension ti(nu,nu,nu),t3(no,no,no,nu,nu,nu)
      do 10 a=1,nu
      do 10 b=1,nu
      do 10 c=1,nu
      t3(i,j,k,a,b,c)=ti(a,b,c)
 10   continue
      return
      end
      subroutine drdequa(i,j,k,l,no,nu,o1,oeh,oep)
      implicit double precision (a-h,o-z)
      logical print
      common/flags/iflags(100)
      dimension o1(1),oeh(no),oep(nu)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/itrat/icycle,mx,icn
      common /totmem/mem
      print=iflags(1).gt.10
      i1=1          !t2new
      i2=i1+no2u2   !t4
      i3=i2+nu4     !voe
      i4=i3+no2u2   !ti
      it=i4+nou2
      if(icycle.eq.1.and.i.eq.2.and.print)write(6,98)mem,it
 98   format('Space usage in dequa available - ',i10,'   used - ',i10)
      call dequa(i,j,k,l,no,nu,o1(i1),o1(i3),o1(i2),o1(i4),oeh,oep)
      return
      end
      SUBROUTINE DEQUA(I,J,K,L,NO,NU,T2NEW,VOE,T4,TI,OEH,OEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      DIMENSION T2NEW(NO,NU,NU,NO),VOE(NO,NU,NU,NO),T4(NU,NU,NU,NU),
     *TI(NU,NU,NO),OEH(NO),OEP(NU)
      DATA ZERO/0.0D+0/,TRESH/1.0D-12/,TWO/2.0D+00/,FOUR/4.0D+00/,
     *HALF/0.5D+00/
      DABS(XXX)=ABS(XXX)
      call ienter(5)
      call ro2hpp(1,no,nu,ti,voe)
      DO 95 A=1,NU
      DO 95 B=1,NU
      XIJ=ZERO
      XIK=ZERO
      XIL=ZERO
      XJK=ZERO
      XJL=ZERO
      XKL=ZERO
      DO 92 C=1,NU
      DO 92 D=1,NU
      ABCD=T4(A,B,C,D)
      ABDC=T4(A,B,D,C)
      ACBD=T4(A,C,B,D)
      ACDB=T4(A,C,D,B)
      ADBC=T4(A,D,B,C)
      ADCB=T4(A,D,C,B)
      BACD=T4(B,A,C,D)
      BADC=T4(B,A,D,C)
      BCAD=T4(B,C,A,D)
      BCDA=T4(B,C,D,A)
      BDAC=T4(B,D,A,C)
      BDCA=T4(B,D,C,A)
      CABD=T4(C,A,B,D)
      CADB=T4(C,A,D,B)
      CBAD=T4(C,B,A,D)
      CBDA=T4(C,B,D,A)
      CDAB=T4(C,D,A,B)
      CDBA=T4(C,D,B,A)
      DABC=T4(D,A,B,C)
      DACB=T4(D,A,C,B)
      DBAC=T4(D,B,A,C)
      DBCA=T4(D,B,C,A)
      DCAB=T4(D,C,A,B)
      DCBA=T4(D,C,B,A)
      VKCDL=VOE(K,C,D,L)
      VKDCL=VOE(K,D,C,L)
      VJCDL=VOE(J,C,D,L)
      VJDCL=VOE(J,D,C,L)
      VJCDK=VOE(J,C,D,K)
      VJDCK=VOE(J,D,C,K)
      VICDL=VOE(I,C,D,L)
      VIDCL=VOE(I,D,C,L)
      VICDK=VOE(I,C,D,K)
      VIDCK=VOE(I,D,C,K)
      VICDJ=VOE(I,C,D,J)
      VIDCJ=VOE(I,D,C,J)
      X1=ZERO
      X2=ZERO
      X3=ZERO
      X4=ZERO
      X5=ZERO
      X6=ZERO
      X7=ZERO
      X8=ZERO
      X9=ZERO
      X10=ZERO
      X11=ZERO
      X12=ZERO
      IF(DABS(VKCDL).GT.TRESH) 
     *X1=(TWO*(ABCD-ACBD-CBAD)-ABDC+ACDB+CBDA+HALF*(CDAB+DCBA))*VKCDL
      IF (J.NE.K.AND.K.NE.L.AND.DABS(VKDCL).GT.TRESH)
     *X2=(TWO*(ABDC-ACDB-CBDA)-ABCD+ACBD+CBAD+HALF*(CDBA+DCAB))*VKDCL
      XIJ=XIJ+X1+X2
      IF (J.NE.K.AND.K.NE.L.AND.DABS(VJCDL).GT.TRESH)
     *X3=(TWO*(ACBD-ABCD-CABD)-ADBC+ADCB+CDBA+HALF*(CADB+DBCA))*VJCDL
      IF(DABS(VJDCL).GT.TRESH) 
     *X5=(TWO*(ADBC-ADCB-CDBA)-ACBD+ABCD+CABD+HALF*(CBDA+DACB))*VJDCL
      XIK=XIK+X3+X5
      IF(DABS(VJCDK).GT.TRESH) 
     *X4=(TWO*(ACDB-ABDC-CADB)-ADCB+ADBC+CDAB+HALF*(CABD+DBAC))*VJCDK
      IF (J.NE.K.AND.K.NE.L.AND.DABS(VJDCK).GT.TRESH)
     *X6=(TWO*(ADCB-ADBC-CDAB)-ACDB+ABDC+CADB+HALF*(CBAD+DABC))*VJDCK
      XIL=XIL+X4+X6
      IF (I.NE.J.AND.K.NE.L.AND.DABS(VICDL).GT.TRESH)
     *X7=(TWO*(CABD-BACD-ACBD)-DABC+DACB+DCBA+HALF*(ACDB+BDCA))*VICDL
      IF (I.NE.J.AND.DABS(VIDCL).GT.TRESH)
     *X8=(TWO*(DABC-DACB-DCBA)-CABD+BACD+ACBD+HALF*(BCDA+ADCB))*VIDCL
      XJK=XJK+X7+X8
      IF (I.NE.J.AND.J.NE.K.AND.DABS(VICDK).GT.TRESH)
     *X9=(TWO*(CADB-BADC-ACDB)-DACB+DABC+DCAB+HALF*(ACBD+BDAC))*VICDK
      IF (I.NE.J.AND.K.NE.L.AND.DABS(VIDCK).GT.TRESH)
     *X10=(TWO*(DACB-DABC-DCAB)-CADB+BADC+ACDB+HALF*(BCAD+ADBC))*VIDCK
      XJL=XJL+X9+X10
      IF (DABS(VICDJ).GT.TRESH)
     *X11=(TWO*(CDAB-BDAC-ADCB)-DCAB+DBAC+DACB+HALF*(ABCD+BADC))*VICDJ
      IF (I.NE.J.AND.J.NE.K.AND.DABS(VIDCJ).GT.TRESH)
     *X12=(TWO*(DCAB-DBAC-DACB)-CDAB+BDAC+ADCB+HALF*(BACD+ABDC))*VIDCJ
      XKL=XKL+X11+X12
      XSUM=X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12
 1236 FORMAT(6F12.8)
   92 CONTINUE
      T2NEW(I,A,B,J)=T2NEW(I,A,B,J)+XIJ
      T2NEW(I,A,B,K)=T2NEW(I,A,B,K)+XIK
      T2NEW(I,A,B,L)=T2NEW(I,A,B,L)+XIL
      T2NEW(J,A,B,K)=T2NEW(J,A,B,K)+XJK
      T2NEW(J,A,B,L)=T2NEW(J,A,B,L)+XJL
      T2NEW(K,A,B,L)=T2NEW(K,A,B,L)+XKL
   95 CONTINUE
      call iexit(5)
      RETURN
      END
      subroutine drenrcon(no,nu,o1,oeh,oep,esdtq,x0)
      implicit double precision (a-h,o-z)
      logical print
      common/flags/iflags(100)
      dimension o1(1),oeh(no),oep(nu)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common /totmem/mem
      common/itrat/icycle,mx,icn
      print=iflags(1).gt.10
      i1=1          !t2new
      i2=i1+no2u2   !ti
      i3=i2+nou2    !t2vo
      i4=i3+no2u2   !voe
      i5=i4+no2u2   !vo
      it=i5+no2u2
      if(icycle.eq.1.and.print)write(6,98)mem,it
 98   format('Space usage in enrcon available - ',i10,'   used - ',i10)
      call enrcon(no,nu,o1(i2),o1,o1(i3),o1(i4),o1(i5),oeh,oep,esdtq,x0)
      return
      end
      SUBROUTINE ENRCON(NO,NU,TI,T2NEW,T2VO,VOE,VO,OEH,OEP,ESDTQTOT,x0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*8 met,ccopt
      INTEGER A,B
      LOGICAL Q1
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      common/timeinfo/timein,timenow,timetot,timenew
      COMMON/NEWOPT/NOPT(6)
      COMMON/ECCCONT/CONECC(20)
      COMMON/OPTSD/NOPTSD(6)
      COMMON/ENRSDT/ECCSDT
      common/rstrt/irest,ccopt(10),t32
      common/finerg/ergfin0,ergfin1
      COMMON/SCFEN/ESCF
      COMMON/ENERGIES/ENRES(150)
      COMMON/ITRAT/IT,MAXIT,ICCNV
      COMMON/ENERGY/ECORR(500,2),IXTRLE(500)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION T2NEW(NO,NU,NU,NO),VOE(NO,NU,NU,NO),VO(NU,NU,NO,NO),
     *T2VO(NO,NU,NU,NO),OEH(NO),OEP(NU),ti(1)
      DATA ZERO/0.0D+0/,TRESH/1.0D-12/,TWO/2.0D+00/,FOUR/4.0D+00/,
     *     HALF/0.5D+00/,ONE/1.0D+0/
      CALL IENTER(6)
      met=ccopt(nopt(1)+2)
      Q1=NOPT(1).GT.2.AND.NOPT(6).NE.9.
     &          OR.NOPT(6).EQ.8.AND.NOPTSD(2).GE.8
      IF(.NOT.Q1)THEN
      XE=ZERO
      CALL ZEROMA(T2NEW,1,NO2U2)
      GOTO 111
      ENDIF
      DO 49 I=1,NO
      I1=I-1
      DO 49 A=1,NU
      DO 49 B=1,NU
      DO 49 J=1,I1
      T2NEW(J,A,B,I)=T2NEW(I,B,A,J)
 49   CONTINUE
      CALL RO2HPP(1,NO,NU,TI,VOE)
      CALL RO2HPP(0,NO,NU,TI,T2VO)
      CALL RO2PPH(1,NO,NU,TI,vo)
      XE=ZERO
      XE1=ZERO
      XE2=ZERO
      DO 150 I=1,NO
      DO 150 J=1,NO
      DO 150 A=1,NU
      DO 150 B=1,NU
      DENOM=OEH(I)+OEH(J)-OEP(A)-OEP(B)
      T2AB=T2NEW(I,A,B,J)
      T2BA=T2NEW(I,B,A,J)
      XE =XE + VOE(I,A,B,J)*(TWO*T2AB-T2BA)/DENOM
      XE1=XE1+T2VO(I,A,B,J)*(TWO*T2AB-T2BA)
      XE2=XE2+  VO(A,B,J,I)*(TWO*T2AB-T2BA)
 150  CONTINUE
      CALL T2DEN(T2NEW,TI,OEH,OEP,NO,NU)
 111  CONTINUE
      EXT4=XE
      if(nopt(6).eq.8)ext4=xe1
      ESDTQTOT=EXT4+ECCSDT
      ETOT=ESDTQTOT+ESCF
      CALL PUTREC(20,'JOBARC','TOTENERG',IINTFP,ETOT)
      ecorr(it+1,2)=ecorr(it+1,2)+xe
 152  FORMAT('It=',I3,a8,'energy=',F15.10,' Total=',F15.10,' cpu=',I5's'
     *)
      call timer(1)
      ii0=timenow-x0
      WRITE(6,152)IT,met,ESDTQTOT,(ESDTQTOT+ESCF),ii0
      ergfin0=esdtqtot
      ergfin1=esdtqtot+escf
      ENRES(IT)=ESDTQTOT
 2000 FORMAT('   **',T5,'W   :',F14.10,'               **')
 2001 FORMAT('   **',T5,'T2  :',F14.10,'               **')
 2002 FORMAT('   **',T5,'T(2):',F14.10,'               **')
      EQ=XE
      CONECC(4)=XE1
      CALL ADT2T4D(NO,NU,TI,T2NEW,VOE)
      CALL IEXIT(6)
      RETURN
      END
      SUBROUTINE MATMUL(A,B,C,NI,NJ,NK,IACC,IFUN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      call ienter(16)
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
      call iexit(16)
      RETURN
      END
      SUBROUTINE MTRANS(V,NU,ID)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      DIMENSION V(NU,NU,NU,NU)
      call ienter(19)
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,
     *22,23),ID
    1 CONTINUE
      DO 100 K=1,NU
      DO 100 L=1,K
      DO 100 J=1,NU
      DO 100 I=1,NU
      X=V(I,J,K,L)
      V(I,J,K,L)=V(I,J,L,K)
      V(I,J,L,K)=X
  100 CONTINUE
      GO TO 1000
    2 CONTINUE
      DO 200 B=1,NU
      DO 200 C=1,B
      DO 200 D=1,C
      DO 200 A=1,NU
      X=V(A,D,B,C)
      V(A,D,B,C)=V(A,B,C,D)
      V(A,B,C,D)=V(A,C,D,B)
      V(A,C,D,B)=X
      IF(B.EQ.C.OR.C.EQ.D)GO TO 200
      X=V(A,B,D,C)
      V(A,B,D,C)=V(A,D,C,B)
      V(A,D,C,B)=V(A,C,B,D)
      V(A,C,B,D)=X
  200 CONTINUE
      GO TO 1000
    3 CONTINUE
      DO 300 B=1,NU
      DO 300 C=1,B
      DO 300 D=1,C
      DO 300 A=1,NU
      X=V(A,C,D,B)
      V(A,C,D,B)=V(A,B,C,D)
      V(A,B,C,D)=V(A,D,B,C)
      V(A,D,B,C)=X
      IF (B.EQ.C.OR.C.EQ.D)GO TO 300
      X=V(A,B,D,C)
      V(A,B,D,C)=V(A,C,B,D)
      V(A,C,B,D)=V(A,D,C,B)
      V(A,D,C,B)=X
  300 CONTINUE
      GO TO 1000
    4 CONTINUE
	call mtranrec(v,nu,3)
	call mtranrec(v,nu,8)
	goto 1000
      DO 400 A=1,NU
      DO 400 B=1,A
      DO 400 C=1,A
      LC=C
      IF (A.EQ.C)LC=B
      DO 400 D=1,LC
      X=V(D,C,A,B)
      V(D,C,A,B)=V(A,B,C,D)
      V(A,B,C,D)=V(C,D,B,A)
      V(C,D,B,A)=V(B,A,D,C)
      V(B,A,D,C)=X
      IF (C.EQ.D.OR.A.EQ.C.AND.B.EQ.D.OR.A.EQ.B)GO TO 400
      X=V(C,D,A,B)
      V(C,D,A,B)=V(A,B,D,C)
      V(A,B,D,C)=V(D,C,B,A)
      V(D,C,B,A)=V(B,A,C,D)
      V(B,A,C,D)=X
  400 CONTINUE
      GO TO 1000
    5 CONTINUE
	call mtranrec(v,nu,8)
	call mtranrec(v,nu,1)
	goto 1000
      DO 500 A=1,NU
      DO 500 C=1,A
      DO 500 D=1,C
      DO 500 B=1,NU
      X=V(C,B,D,A)
      V(C,B,D,A)=V(A,B,C,D)
      V(A,B,C,D)=V(D,B,A,C)
      V(D,B,A,C)=X
      IF (A.EQ.C.OR.C.EQ.D)GO TO 500
      X=V(A,B,D,C)
      V(A,B,D,C)=V(C,B,A,D)
      V(C,B,A,D)=V(D,B,C,A)
      V(D,B,C,A)=X
  500 CONTINUE
      GO TO 1000
    6 CONTINUE
      DO 600 D=1,NU
      DO 600 C=1,NU
      DO 600 A=1,NU
      DO 600 B=1,A
      X=V(B,A,C,D)
      V(B,A,C,D)=V(A,B,C,D)
      V(A,B,C,D)=X
  600 CONTINUE
c      do 600 d=1,nu
c      do 600 c=1,nu
c      call transq(v(1,1,c,d),nu)
c 600  continue
      GO TO 1000
    7 CONTINUE
      DO 700 D=1,NU
      DO 700 B=1,NU
      DO 700 C=1,B
      DO 700 A=1,NU
      X=V(A,B,C,D)
      V(A,B,C,D)=V(A,C,B,D)
      V(A,C,B,D)=X
  700 CONTINUE
      GO TO 1000
    8 CONTINUE
      DO 800 D=1,NU
      DO 800 B=1,NU
      DO 800 A=1,NU
      DO 800 C=1,A
      X=V(A,B,C,D)
      V(A,B,C,D)=V(C,B,A,D)
      V(C,B,A,D)=X
  800 CONTINUE
      GO TO 1000
    9 CONTINUE
	call mtra1(v,nu,10)
	call mtranrec(v,nu,8)
	goto 1000
      DO 900 A=1,NU
      DO 900 B=1,NU
      DO 900 C=1,A
      LIMD=NU
      IF (A.EQ.C)LIMD=B
      DO 900 D=1,LIMD
      X=V(A,B,C,D)
      V(A,B,C,D)=V(C,D,A,B)
      V(C,D,A,B)=X
  900 CONTINUE
c      nu2=nu*nu
c      call transq1(v,nu2)
      GO TO 1000
 10   CONTINUE
	call mtra1(v,nu,10)
	goto 1000
      DO 950 A=1,NU
      DO 950 B=1,NU
      DO 950 C=1,NU
      DO 950 D=1,B
      X=V(A,B,C,D)
      V(A,B,C,D)=V(A,D,C,B)
      V(A,D,C,B)=X
 950  CONTINUE
      GOTO 1000
 11   CONTINUE
	call mtranrec(v,nu,1)
	call mtranrec(v,nu,8)
	call mtranrec(v,nu,1)
	goto 1000
      DO 960 A=1,NU
      DO 960 B=1,NU
      DO 960 C=1,NU
      DO 960 D=1,A
      X=V(A,B,C,D)
      V(A,B,C,D)=V(D,B,C,A)
      V(D,B,C,A)=X
 960  CONTINUE
      GOTO 1000
 12   CONTINUE
	call mtranrec(v,nu,7)
	call mtranrec(v,nu,6)
	goto 1000
      call tranmd(v,nu,nu,nu,nu,231)      
      GOTO 1000
 13   CONTINUE
	call mtranrec(v,nu,6)
	call mtranrec(v,nu,7)
	goto 1000
      call tranmd(v,nu,nu,nu,nu,312)
      GOTO 1000
 14   CONTINUE
	call mtranrec(v,nu,1)
	call mtranrec(v,nu,8)
      GO TO 1000
      DO 970 A=1,NU
      DO 970 B=1,NU
      DO 970 C=1,A
      DO 970 D=1,C
      X=V(D,B,A,C)
      V(D,B,A,C)=V(A,B,C,D)
      V(A,B,C,D)=V(C,B,D,A)
      V(C,B,D,A)=X
      IF (A.EQ.C.OR.C.EQ.D)GO TO 970
      X=V(D,B,C,A)
      V(D,B,C,A)=V(C,B,A,D)
      V(C,B,A,D)=V(A,B,D,C)
      V(A,B,D,C)=X
  970 CONTINUE
      GO TO 1000
 15   CONTINUE
      call mtra1(v,nu,10)
      call mtranrec(v,nu,6)
      GO TO 1000
      DO 980 A=1,NU
      DO 980 B=1,A
      DO 980 C=1,NU
      DO 980 D=1,B
      X=V(D,A,C,B)
      V(D,A,C,B)=V(A,B,C,D)
      V(A,B,C,D)=V(B,D,C,A)
      V(B,D,C,A)=X
      IF (A.EQ.B.OR.B.EQ.D)GO TO 980
      X=V(D,B,C,A)
      V(D,B,C,A)=V(B,A,C,D)
      V(B,A,C,D)=V(A,D,C,B)
      V(A,D,C,B)=X
 980  CONTINUE
      GO TO 1000
 16   CONTINUE
      call mtranrec(v,nu,6)
      call mtra1(v,nu,10)
      GO TO 1000
      DO 990 A=1,NU
      DO 990 B=1,A
      DO 990 C=1,NU
      DO 990 D=1,B
      X=V(B,D,C,A)
      V(B,D,C,A)=V(A,B,C,D)
      V(A,B,C,D)=V(D,A,C,B)
      V(D,A,C,B)=X
      IF (A.EQ.B.OR.B.EQ.D)GO TO 990
      X=V(A,D,C,B)
      V(A,D,C,B)=V(B,A,C,D)
      V(B,A,C,D)=V(D,B,C,A)
      V(D,B,C,A)=X
  990 CONTINUE
      GO TO 1000
 17   CONTINUE
      call mtranrec(v,nu,8)
      call mtranrec(v,nu,2)
      GO TO 1000
      DO 991 A=1,NU
      DO 991 B=1,A
      DO 991 C=1,A
      LC=C
      IF (A.EQ.C)LC=B
      DO 991 D=1,LC
      X=V(A,B,C,D)
      V(A,B,C,D)=V(D,C,A,B)
      V(D,C,A,B)=V(B,A,D,C)
      V(B,A,D,C)=V(C,D,B,A)
      V(C,D,B,A)=X
      IF (C.EQ.D.OR.A.EQ.C.AND.B.EQ.D.OR.A.EQ.B)GO TO 991
      X=V(A,B,D,C)
      V(A,B,D,C)=V(C,D,A,B)
      V(C,D,A,B)=V(B,A,C,D)
      V(B,A,C,D)=V(D,C,B,A)
      V(D,C,B,A)=X
 991  CONTINUE
      GO TO 1000
 18   CONTINUE
	call mtranrec(v,nu,1)
	call mtranrec(v,nu,8)
	call mtranrec(v,nu,2)
      GO TO 1000
      DO 992 A=1,NU
      DO 992 B=1,NU
      DO 992 D=1,A
      LC=NU
      IF (A.EQ.d)LC=B
      DO 992 c=1,LC
      X=V(A,B,C,D)
      V(A,B,C,D)=V(D,C,B,A)
      V(D,C,B,A)=X
 992  CONTINUE
      GO TO 1000
 19   CONTINUE
	call mtra1(v,nu,19)
	goto 1000
      DO 993 A=1,NU
c      write(6,*)'a=',a
      DO 993 B=1,a
      DO 993 C=1,nu
      DO 993 D=1,c
      X=V(A,B,C,D)
      V(A,B,C,D)=V(B,A,D,C)
      V(B,A,D,C)=X
      if(a.eq.b.or.c.eq.d)goto 9993
      x=v(a,b,d,c)
      v(a,b,d,c)=v(b,a,c,d)
      v(b,a,c,d)=x
 9993 continue
c      write(6,909)a,b,c,d,v(a,b,c,d),b,a,d,c,v(b,a,d,c)
 993  CONTINUE
 909  format('a,b,c,d,v:',4i3,i8,5x,4i3,i8)
      GO TO 1000
 20   CONTINUE
      call mtranrec(v,nu,2)
      call mtranrec(v,nu,6)
      GO TO 1000
      DO 994 A=1,NU
      DO 994 B=1,A
      DO 994 C=1,A
      LD=B
      IF (A.EQ.B)Ld=C
      DO 994 D=1,LD
      X=V(A,B,C,D)
      V(A,B,C,D)=v(B,C,D,A)
      V(B,C,D,A)=V(C,D,A,B)
      V(C,D,A,B)=V(D,A,B,C)
      V(D,A,B,C)=X
      IF (B.EQ.D.OR.A.EQ.B.AND.C.EQ.D.OR.A.EQ.C)GO TO 994
      X=V(A,D,C,B)
      V(A,D,C,B)=v(D,C,B,A)
      V(D,C,B,A)=v(C,B,A,D)
      V(C,B,A,D)=v(B,A,D,C)
      V(B,A,D,C)=X
 994  CONTINUE
      GO TO 1000
 21   CONTINUE
      call mtranrec(v,nu,6)
      call mtranrec(v,nu,3)
      GO TO 1000
      DO 995 A=1,NU
      DO 995 B=1,A
      DO 995 C=1,A
      LD=B
      IF (A.EQ.B)LD=C
      DO 995 D=1,LD
      X=V(A,B,C,D)
      V(A,B,C,D)=V(D,A,B,C)
      V(D,A,B,C)=V(C,D,A,B)
      V(C,D,A,B)=V(B,C,D,A)
      V(B,C,D,A)=X
      IF (B.EQ.D.OR.A.EQ.B.AND.C.EQ.D.OR.A.EQ.C)GO TO 995
      X=V(A,D,C,B)
      V(A,D,C,B)=V(B,A,D,C)
      V(B,A,D,C)=V(C,B,A,D)
      V(C,B,A,D)=V(D,C,B,A)
      V(D,C,B,A)=X
 995  CONTINUE
      GO TO 1000
 22   CONTINUE
      call mtranrec(v,nu,3)
      call mtranrec(v,nu,6)
      GO TO 1000
      DO 996 A=1,NU
      DO 996 B=1,A
      DO 996 D=1,A
      LC=B
      IF (A.EQ.B)LC=D
      DO 996 C=1,LC
      X=V(A,B,C,D)
      V(A,B,C,D)=V(B,D,A,C)
      V(B,D,A,C)=V(D,C,B,A)
      V(D,C,B,A)=V(C,A,D,B)
      V(C,A,D,B)=X
      IF (B.EQ.C.OR.A.EQ.B.AND.C.EQ.D.OR.A.EQ.D)GO TO 996
      X=V(A,C,B,D)
      V(A,C,B,D)=V(C,D,A,B)
      V(C,D,A,B)=V(D,B,C,A)
      V(D,B,C,A)=V(B,A,D,C)
      V(B,A,D,C)=X
 996  CONTINUE
      GO TO 1000
 23   CONTINUE
      call mtranrec(v,nu,2)
      call mtranrec(v,nu,8)
      GO TO 1000
      DO 997 A=1,NU
      DO 997 B=1,A
      DO 997 D=1,A
      LC=B
      IF (A.EQ.B)LC=D
      DO 997 C=1,LC
      X=V(A,B,C,D)
      V(A,B,C,D)=V(C,A,D,B)
      V(C,A,D,B)=V(D,C,B,A)
      V(D,C,B,A)=V(B,D,A,C)
      V(B,D,A,C)=X
      IF (B.EQ.C.OR.A.EQ.B.AND.C.EQ.D.OR.A.EQ.D)GO TO 997
      X=V(A,C,B,D)
      V(A,C,B,D)=V(B,A,D,C)
      V(B,A,D,C)=V(D,B,C,A)
      V(D,B,C,A)=V(C,D,A,B)
      V(C,D,A,B)=X
 997  CONTINUE
      GO TO 1000
 1000 continue
      call iexit(19)
      RETURN
      END
      SUBROUTINE MTRANSnew(V,NU,ID)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      DIMENSION V(NU,NU,NU,NU)
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,
     *22,23),ID
    1 CONTINUE
      DO 100 K=1,NU
      DO 100 L=1,K
      DO 100 J=1,NU
      DO 100 I=1,NU
      X=V(I,J,K,L)
      V(I,J,K,L)=V(I,J,L,K)
      V(I,J,L,K)=X
  100 CONTINUE
      GO TO 1000
    2 CONTINUE
      DO 200 B=1,NU
      DO 200 C=1,B
      DO 200 D=1,C
      DO 200 A=1,NU
      X=V(A,D,B,C)
      V(A,D,B,C)=V(A,B,C,D)
      V(A,B,C,D)=V(A,C,D,B)
      V(A,C,D,B)=X
      IF(B.EQ.C.OR.C.EQ.D)GO TO 200
      X=V(A,B,D,C)
      V(A,B,D,C)=V(A,D,C,B)
      V(A,D,C,B)=V(A,C,B,D)
      V(A,C,B,D)=X
  200 CONTINUE
      GO TO 1000
    3 CONTINUE
      DO 300 B=1,NU
      DO 300 C=1,B
      DO 300 D=1,C
      DO 300 A=1,NU
      X=V(A,C,D,B)
      V(A,C,D,B)=V(A,B,C,D)
      V(A,B,C,D)=V(A,D,B,C)
      V(A,D,B,C)=X
      IF (B.EQ.C.OR.C.EQ.D)GO TO 300
      X=V(A,B,D,C)
      V(A,B,D,C)=V(A,C,B,D)
      V(A,C,B,D)=V(A,D,C,B)
      V(A,D,C,B)=X
  300 CONTINUE
      GO TO 1000
    4 CONTINUE
	call mtranrec(v,nu,3)
	call mtranrec(v,nu,8)
	goto 1000
      DO 400 A=1,NU
      DO 400 B=1,A
      DO 400 C=1,A
      LC=C
      IF (A.EQ.C)LC=B
      DO 400 D=1,LC
      X=V(D,C,A,B)
      V(D,C,A,B)=V(A,B,C,D)
      V(A,B,C,D)=V(C,D,B,A)
      V(C,D,B,A)=V(B,A,D,C)
      V(B,A,D,C)=X
      IF (C.EQ.D.OR.A.EQ.C.AND.B.EQ.D.OR.A.EQ.B)GO TO 400
      X=V(C,D,A,B)
      V(C,D,A,B)=V(A,B,D,C)
      V(A,B,D,C)=V(D,C,B,A)
      V(D,C,B,A)=V(B,A,C,D)
      V(B,A,C,D)=X
  400 CONTINUE
      GO TO 1000
    5 CONTINUE
	call mtranrec(v,nu,8)
	call mtranrec(v,nu,1)
	goto 1000
      DO 500 A=1,NU
      DO 500 C=1,A
      DO 500 D=1,C
      DO 500 B=1,NU
      X=V(C,B,D,A)
      V(C,B,D,A)=V(A,B,C,D)
      V(A,B,C,D)=V(D,B,A,C)
      V(D,B,A,C)=X
      IF (A.EQ.C.OR.C.EQ.D)GO TO 500
      X=V(A,B,D,C)
      V(A,B,D,C)=V(C,B,A,D)
      V(C,B,A,D)=V(D,B,C,A)
      V(D,B,C,A)=X
  500 CONTINUE
      GO TO 1000
    6 CONTINUE
      DO 600 D=1,NU
      DO 600 C=1,NU
      DO 600 A=1,NU
      DO 600 B=1,A
      X=V(B,A,C,D)
      V(B,A,C,D)=V(A,B,C,D)
      V(A,B,C,D)=X
  600 CONTINUE
      GO TO 1000
    7 CONTINUE
      DO 700 D=1,NU
      DO 700 B=1,NU
      DO 700 C=1,B
      DO 700 A=1,NU
      X=V(A,B,C,D)
      V(A,B,C,D)=V(A,C,B,D)
      V(A,C,B,D)=X
  700 CONTINUE
      GO TO 1000
    8 CONTINUE
      DO 800 D=1,NU
      DO 800 B=1,NU
      DO 800 A=1,NU
      DO 800 C=1,A
      X=V(A,B,C,D)
      V(A,B,C,D)=V(C,B,A,D)
      V(C,B,A,D)=X
  800 CONTINUE
      GO TO 1000
    9 CONTINUE
	call mtra1(v,nu,10)
	call mtranrec(v,nu,8)
	goto 1000
      DO 900 A=1,NU
      DO 900 B=1,NU
      DO 900 C=1,A
      LIMD=NU
      IF (A.EQ.C)LIMD=B
      DO 900 D=1,LIMD
      X=V(A,B,C,D)
      V(A,B,C,D)=V(C,D,A,B)
      V(C,D,A,B)=X
  900 CONTINUE
      GO TO 1000
 10   CONTINUE
	call mtra1(v,nu,10)
	goto 1000
      DO 950 A=1,NU
      DO 950 B=1,NU
      DO 950 C=1,NU
      DO 950 D=1,B
      X=V(A,B,C,D)
      V(A,B,C,D)=V(A,D,C,B)
      V(A,D,C,B)=X
 950  CONTINUE
      GOTO 1000
 11   CONTINUE
	call mtranrec(v,nu,1)
	call mtranrec(v,nu,8)
	call mtranrec(v,nu,1)
	goto 1000
      DO 960 A=1,NU
      DO 960 B=1,NU
      DO 960 C=1,NU
      DO 960 D=1,A
      X=V(A,B,C,D)
      V(A,B,C,D)=V(D,B,C,A)
      V(D,B,C,A)=X
 960  CONTINUE
      GOTO 1000
 12   CONTINUE
	call mtranrec(v,nu,7)
	call mtranrec(v,nu,6)
	goto 1000
      call tranmd(v,nu,nu,nu,nu,231)      
      GOTO 1000
 13   CONTINUE
	call mtranrec(v,nu,6)
	call mtranrec(v,nu,7)
	goto 1000
      call tranmd(v,nu,nu,nu,nu,312)
      GOTO 1000
 14   CONTINUE
	call mtranrec(v,nu,1)
	call mtranrec(v,nu,8)
      GO TO 1000
      DO 970 A=1,NU
      DO 970 B=1,NU
      DO 970 C=1,A
      DO 970 D=1,C
      X=V(D,B,A,C)
      V(D,B,A,C)=V(A,B,C,D)
      V(A,B,C,D)=V(C,B,D,A)
      V(C,B,D,A)=X
      IF (A.EQ.C.OR.C.EQ.D)GO TO 970
      X=V(D,B,C,A)
      V(D,B,C,A)=V(C,B,A,D)
      V(C,B,A,D)=V(A,B,D,C)
      V(A,B,D,C)=X
  970 CONTINUE
      GO TO 1000
 15   CONTINUE
      call mtra1(v,nu,10)
      call mtranrec(v,nu,6)
      GO TO 1000
      DO 980 A=1,NU
      DO 980 B=1,A
      DO 980 C=1,NU
      DO 980 D=1,B
      X=V(D,A,C,B)
      V(D,A,C,B)=V(A,B,C,D)
      V(A,B,C,D)=V(B,D,C,A)
      V(B,D,C,A)=X
      IF (A.EQ.B.OR.B.EQ.D)GO TO 980
      X=V(D,B,C,A)
      V(D,B,C,A)=V(B,A,C,D)
      V(B,A,C,D)=V(A,D,C,B)
      V(A,D,C,B)=X
 980  CONTINUE
      GO TO 1000
 16   CONTINUE
      call mtranrec(v,nu,6)
      call mtra1(v,nu,10)
      GO TO 1000
      DO 990 A=1,NU
      DO 990 B=1,A
      DO 990 C=1,NU
      DO 990 D=1,B
      X=V(B,D,C,A)
      V(B,D,C,A)=V(A,B,C,D)
      V(A,B,C,D)=V(D,A,C,B)
      V(D,A,C,B)=X
      IF (A.EQ.B.OR.B.EQ.D)GO TO 990
      X=V(A,D,C,B)
      V(A,D,C,B)=V(B,A,C,D)
      V(B,A,C,D)=V(D,B,C,A)
      V(D,B,C,A)=X
  990 CONTINUE
      GO TO 1000
 17   CONTINUE
      call mtranrec(v,nu,8)
      call mtranrec(v,nu,2)
      GO TO 1000
      DO 991 A=1,NU
      DO 991 B=1,A
      DO 991 C=1,A
      LC=C
      IF (A.EQ.C)LC=B
      DO 991 D=1,LC
      X=V(A,B,C,D)
      V(A,B,C,D)=V(D,C,A,B)
      V(D,C,A,B)=V(B,A,D,C)
      V(B,A,D,C)=V(C,D,B,A)
      V(C,D,B,A)=X
      IF (C.EQ.D.OR.A.EQ.C.AND.B.EQ.D.OR.A.EQ.B)GO TO 991
      X=V(A,B,D,C)
      V(A,B,D,C)=V(C,D,A,B)
      V(C,D,A,B)=V(B,A,C,D)
      V(B,A,C,D)=V(D,C,B,A)
      V(D,C,B,A)=X
 991  CONTINUE
      GO TO 1000
 18   CONTINUE
	call mtranrec(v,nu,1)
	call mtranrec(v,nu,8)
	call mtranrec(v,nu,2)
      GO TO 1000
      DO 992 A=1,NU
      DO 992 B=1,NU
      DO 992 D=1,A
      LC=NU
      IF (A.EQ.d)LC=B
      DO 992 c=1,LC
      X=V(A,B,C,D)
      V(A,B,C,D)=V(D,C,B,A)
      V(D,C,B,A)=X
 992  CONTINUE
      GO TO 1000
 19   CONTINUE
	call mtra1(v,nu,19)
	goto 1000
      DO 993 A=1,NU
c      write(6,*)'a=',a
      DO 993 B=1,a
      DO 993 C=1,nu
      DO 993 D=1,c
      X=V(A,B,C,D)
      V(A,B,C,D)=V(B,A,D,C)
      V(B,A,D,C)=X
      if(a.eq.b.or.c.eq.d)goto 9993
      x=v(a,b,d,c)
      v(a,b,d,c)=v(b,a,c,d)
      v(b,a,c,d)=x
 9993 continue
c      write(6,909)a,b,c,d,v(a,b,c,d),b,a,d,c,v(b,a,d,c)
 993  CONTINUE
 909  format('a,b,c,d,v:',4i3,i8,5x,4i3,i8)
      GO TO 1000
 20   CONTINUE
      call mtranrec(v,nu,2)
      call mtranrec(v,nu,6)
      GO TO 1000
      DO 994 A=1,NU
      DO 994 B=1,A
      DO 994 C=1,A
      LD=B
      IF (A.EQ.B)Ld=C
      DO 994 D=1,LD
      X=V(A,B,C,D)
      V(A,B,C,D)=v(B,C,D,A)
      V(B,C,D,A)=V(C,D,A,B)
      V(C,D,A,B)=V(D,A,B,C)
      V(D,A,B,C)=X
      IF (B.EQ.D.OR.A.EQ.B.AND.C.EQ.D.OR.A.EQ.C)GO TO 994
      X=V(A,D,C,B)
      V(A,D,C,B)=v(D,C,B,A)
      V(D,C,B,A)=v(C,B,A,D)
      V(C,B,A,D)=v(B,A,D,C)
      V(B,A,D,C)=X
 994  CONTINUE
      GO TO 1000
 21   CONTINUE
      call mtranrec(v,nu,6)
      call mtranrec(v,nu,3)
      GO TO 1000
      DO 995 A=1,NU
      DO 995 B=1,A
      DO 995 C=1,A
      LD=B
      IF (A.EQ.B)LD=C
      DO 995 D=1,LD
      X=V(A,B,C,D)
      V(A,B,C,D)=V(D,A,B,C)
      V(D,A,B,C)=V(C,D,A,B)
      V(C,D,A,B)=V(B,C,D,A)
      V(B,C,D,A)=X
      IF (B.EQ.D.OR.A.EQ.B.AND.C.EQ.D.OR.A.EQ.C)GO TO 995
      X=V(A,D,C,B)
      V(A,D,C,B)=V(B,A,D,C)
      V(B,A,D,C)=V(C,B,A,D)
      V(C,B,A,D)=V(D,C,B,A)
      V(D,C,B,A)=X
 995  CONTINUE
      GO TO 1000
 22   CONTINUE
      call mtranrec(v,nu,3)
      call mtranrec(v,nu,6)
      GO TO 1000
      DO 996 A=1,NU
      DO 996 B=1,A
      DO 996 D=1,A
      LC=B
      IF (A.EQ.B)LC=D
      DO 996 C=1,LC
      X=V(A,B,C,D)
      V(A,B,C,D)=V(B,D,A,C)
      V(B,D,A,C)=V(D,C,B,A)
      V(D,C,B,A)=V(C,A,D,B)
      V(C,A,D,B)=X
      IF (B.EQ.C.OR.A.EQ.B.AND.C.EQ.D.OR.A.EQ.D)GO TO 996
      X=V(A,C,B,D)
      V(A,C,B,D)=V(C,D,A,B)
      V(C,D,A,B)=V(D,B,C,A)
      V(D,B,C,A)=V(B,A,D,C)
      V(B,A,D,C)=X
 996  CONTINUE
      GO TO 1000
 23   CONTINUE
      call mtranrec(v,nu,2)
      call mtranrec(v,nu,8)
      GO TO 1000
      DO 997 A=1,NU
      DO 997 B=1,A
      DO 997 D=1,A
      LC=B
      IF (A.EQ.B)LC=D
      DO 997 C=1,LC
      X=V(A,B,C,D)
      V(A,B,C,D)=V(C,A,D,B)
      V(C,A,D,B)=V(D,C,B,A)
      V(D,C,B,A)=V(B,D,A,C)
      V(B,D,A,C)=X
      IF (B.EQ.C.OR.A.EQ.B.AND.C.EQ.D.OR.A.EQ.D)GO TO 997
      X=V(A,C,B,D)
      V(A,C,B,D)=V(B,A,D,C)
      V(B,A,D,C)=V(D,B,C,A)
      V(D,B,C,A)=V(C,D,A,B)
      V(C,D,A,B)=X
 997  CONTINUE
      GO TO 1000
 1000 continue
c      call iexit(2)
      RETURN
      END
      SUBROUTINE MTRSM(V,NO,NU,ID)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      DIMENSION V(NU,NU,NU,NU,NO)
      call ienter(20)
      DO 1 I=1,NO
      CALL MTRANS(V(1,1,1,1,I),NU,ID)
    1 CONTINUE
      call iexit(20)
      RETURN
      END
      SUBROUTINE MTRSMN(V,NO,NU,ID)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      DIMENSION V(NU,NU,NU,NU,NO,NO)
      call ienter(21)
      DO 1 I=1,NO
      DO 1 J=1,NO
      CALL MTRANS(V(1,1,1,1,I,J),NU,ID)
    1 CONTINUE
      call iexit(21)
      RETURN
      END
      
      SUBROUTINE SYMT3(V,NU,ID)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      DIMENSION V(NU,NU,NU)
      DATA TWO/2.0D+0/
      call ienter(42)
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
      call iexit(42)
      RETURN
      END
      SUBROUTINE SYMT3inew(V,no,NU)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      DIMENSION V(NU,NU,NU,no)
      DATA TWO/2.0D+0/
      do 100 i=1,NO
      DO 100 A=1,NU
      DO 100 B=1,A
      DO 100 C=1,B
      XAB=V(A,B,C,i)
      XAC=V(A,C,B,i)
      XBA=V(B,A,C,i)
      XBC=V(B,C,A,i)
      XCA=V(C,A,B,i)
      XCB=V(C,B,A,i)
      V(A,B,C,i)=TWO*XAB-XAC-XBA
      V(A,C,B,i)=TWO*XAC-XAB-XCA
      V(B,A,C,i)=TWO*XBA-XBC-XAB
      V(B,C,A,i)=TWO*XBC-XBA-XCB
      V(C,A,B,i)=TWO*XCA-XCB-XAC
      V(C,B,A,i)=TWO*XCB-XCA-XBC
  100 CONTINUE
      RETURN
      END
      SUBROUTINE SYMT3I(V,no,nu)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      DIMENSION V(NU,NU,NU,no)
      DATA TWO/2.0D+0/
    1 CONTINUE
      do 100 i=1,no
      DO 100 A=1,NU
      DO 100 B=1,A
      DO 100 C=1,B
         if (a.eq.b.and.b.eq.c)goto 100
      XAB=V(A,B,C,i)
      XAC=V(A,C,B,i)
      XBA=V(B,A,C,i)
      XBC=V(B,C,A,i)
      XCA=V(C,A,B,i)
      XCB=V(C,B,A,i)
      V(A,B,C,i)=TWO*XAB-XAC-XCB
      V(A,C,B,i)=TWO*XAC-XAB-XBC
      V(B,A,C,i)=TWO*XBA-XBC-XCA
      V(B,C,A,i)=TWO*XBC-XBA-XAC
      V(C,A,B,i)=TWO*XCA-XCB-XBA
      V(C,B,A,i)=TWO*XCB-XCA-XAB
  100 CONTINUE
      RETURN
      END
      SUBROUTINE DESYMT3I(V,no,nu)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      DIMENSION V(NU,NU,NU,no)
      DATA TWO/2.0D+0/,zero/0.0d+0/,half/0.5d+0/
    1 CONTINUE
      do 100 i=1,no
      DO 100 A=1,NU
      DO 100 B=1,A
      DO 100 C=1,B
      a1=V(A,B,C,i)
      a2=V(A,C,B,i)
      a3=V(B,A,C,i)
      a4=V(B,C,A,i)
      a5=V(C,A,B,i)
      a6=V(C,B,A,i)
      x1=zero
      x4=half*(a4-a1)
      x2=half*(a2+x4)
      x6=two*x4-x2-a4
      x5=two*x6-a6
      x3=two*x5-x6-a5
      V(A,B,C,i)=x1
      V(A,C,B,i)=x2
      V(B,A,C,i)=x3
      V(B,C,A,i)=x4
      V(C,A,B,i)=x5
      V(C,B,A,i)=x6
  100 CONTINUE
      RETURN
      END
      subroutine drt2sec(no,nu,t,oeh,oep)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical print
      common/flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/totmem/mem
      common/itrat/icycle,mx,icn
      dimension t(*),oeh(no),oep(nu)
      print=iflags(1).gt.10
      i1=1         !ti
      i2=i1+nu3    !t2pp
      i3=i2+no2u2  !t2vo
      i4=i3+no2u2  !voe
      i5=i4+no2u2  !vo
      i6=i5+no2u2  !pz
      i7=i6+no4    !v
      it=i7+nu4    
      if(icycle.eq.1.and.print)write(6,99)mem,it
 99   format('Space usage in t2sec: available - ',i8,'   used - ',i8)
      call t2sec(no,nu,t,t(i2),t(i3),t(i4),t(i5),t(i6),t(i7),oeh,oep)
      return
      end
      SUBROUTINE T2SEC(NO,NU,TI,T2PP,T2VO,VOE,VO,PZ,V,OEH,OEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical print
      common/flags/iflags(100)
      COMMON/NEWTAPE/NT2T4,NT4INT
      COMMON/ITRAT/ITER,MAXIT,ICNV
      COMMON/UNLIN/UN(15)
      INTEGER A,B,E,F
      DIMENSION TI(NU,NU,NO),T2PP(NU,NU,NO,NO),T2VO(NO,NU,NU,NO),
     *VOE(NO,NU,NU,NO),VO(NU,NU,NO,NO),PZ(NO,NO,NO,NO),
     *V(NU,NU,NU,NU),OEH(NO),OEP(NU)
      DATA TWO/2.0D+00/,ZERO/0.0D+0/
      call ienter(54)
      print=iflags(1).gt.12
      NO2=NO*NO
      no4=no2*no2
      NU2=NU*NU
      nu4=nu2*nu2
      NOU2=NO*NU2
      NO2U2=NO2*NU2
      CALL RO2HPP(1,NO,NU,TI,VOE)
      CALL RO2PPH(2,NO,NU,TI,VO)
      CALL RDOV4(1,NU,NO,TI,PZ)
      CALL MTRANS(PZ,NO,9)
      CALL RDOV4(0,NO,NU,TI,V)
      e3h=zero
      e3p=zero
      e3r=zero
      e3rt=zero
      e3=zero
      e2t=zero
      e2tt=zero
      e2e=zero
      e3t=zero
      DO 5 I=1,NO
      DO 5 J=1,NO
      DO 5 A=1,NU
      DO 5 B=1,NU
      DEN=OEH(I)+OEH(J)-OEP(A)-OEP(B)
      T2PP(A,B,J,I)=VOE(I,A,B,J)/DEN
      x2t=voe(i,a,b,j)*(two*voe(i,a,b,j)-voe(i,b,a,j))
      e2t=e2t+x2t/(den*den)
      e2tt=e2tt+x2t/(den*den*den)
      e2e=e2e+x2t/den
 5    CONTINUE
      UN(1)=E2E
      UN(2)=E2T
      UN(3)=E2TT
      E4D=ZERO
      O4D1=ZERO
      O4D2=ZERO
      DO 19 II=1,NO
      DO 9 JJ=1,NO
      DO 9 A=1,NU
      DO 9 B=1,NU
      x1=zero
      x2=zero
      x3=zero
      X1T=ZERO
      X2T=ZERO
      X3T=ZERO
      x3tst=zero
      DEN=OEH(II)+OEH(JJ)-OEP(A)-OEP(B)
      X=ZERO
      DO 8 E=1,NU
      DO 8 F=1,NU
      DENP=OEH(II)+OEH(JJ)-OEP(E)-OEP(F)
      XX=T2PP(E,F,JJ,II)*V(A,E,B,F)
      X1=X1+XX
      X1T=X1T+XX/DENP
 8    CONTINUE
      DO 7 M=1,NO
      DO 7 N=1,NO
      DENH=OEH(M)+OEH(N)-OEP(A)-OEP(B)
      XX=T2PP(A,B,N,M)*PZ(II,JJ,M,N)
      X2=X2+XX
      X2T=X2T+XX/DENH
 7    CONTINUE
      DO 6 M=1,NO
      DO 6 E=1,NU
      DENAI=OEH(II)+OEH(M)-OEP(A)-OEP(E)
      DENAJ=OEH(JJ)+OEH(M)-OEP(A)-OEP(E)
      DENBI=OEH(II)+OEH(M)-OEP(B)-OEP(E)
      DENBJ=OEH(JJ)+OEH(M)-OEP(B)-OEP(E)
      XAI=T2PP(A,E,M,II)*(TWO*VOE(M,E,B,JJ)-VO(B,E,JJ,M))
     *   -T2PP(E,A,M,II)*     VOE(M,E,B,JJ)
      XAJ=T2PP(E,A,M,JJ)*VO(B,E,II,M)
      XBI=T2PP(E,B,M,II)*VO(A,E,JJ,M)
      XBJ=T2PP(B,E,M,JJ)*(TWO*VOE(M,E,A,II)-VO(A,E,II,M))
     *   -T2PP(E,B,M,JJ)   *VOE(M,E,A,II)
C     
      X3TST=X3TST+XAI-XAJ-XBI+XBJ
      X3T=X3T+XAI/DENAI-XAJ/DENAJ-XBI/DENBI+XBJ/DENBJ
      X3=X3+(T2PP(A,E,M,II)*VOE(M,E,B,JJ)+T2PP(B,E,M,JJ)*VOE(M,E,A,II))*
     *TWO  - T2PP(E,A,M,II)*VOE(M,E,B,JJ)-T2PP(E,B,M,JJ)*VOE(M,E,A,II)
     *     - T2PP(A,E,M,II)* VO(B,E,JJ,M)-T2PP(B,E,M,JJ)* VO(A,E,II,M)
     *     - T2PP(E,B,M,II)* VO(A,E,JJ,M)-T2PP(E,A,M,JJ)* VO(B,E,II,M)
 6    CONTINUE
      TI(A,B,JJ)=(X1T+X2T+X3T)/DEN
      T2VO(II,A,B,JJ)=(X1+x2+x3)/DEN
      e3p=e3p+x1*(two*voe(ii,a,b,jj)-voe(ii,b,a,jj))/den
      e3h=e3h+x2*(two*voe(ii,a,b,jj)-voe(ii,b,a,jj))/den
      e3r=e3r+x3*(two*voe(ii,a,b,jj)-voe(ii,b,a,jj))/den
      e3rt=e3rt+x3tst*(two*voe(ii,a,b,jj)-voe(ii,b,a,jj))/den
      e3 =e3 +t2vo(ii,a,b,jj)*(two*voe(ii,a,b,jj)-voe(ii,b,a,jj))
      e3t=e3t+t2vo(ii,a,b,jj)*(two*voe(ii,a,b,jj)-voe(ii,b,a,jj))/den
 9    CONTINUE
      DO 18 JJ=1,NO
      DO 18 A=1,NU
      DO 18 B=1,NU
      DEN=OEH(II)+OEH(JJ)-OEP(A)-OEP(B)
      XX1=(TWO*T2VO(II,A,B,JJ)-T2VO(II,B,A,JJ))
      XX=T2VO(II,A,B,JJ)*XX1
      E4D=E4D+XX*DEN
      O4D2=O4D2+XX
      O4D1=O4D1+TI(A,B,JJ)*XX1*DEN
 18   continue
 19   CONTINUE
      UN(9)=E4D
      UN(10)=O4D1
      UN(11)=O4D2
c      WRITE(6,*)'TEST OF E4D', E4D
      if (iter.eq.1.and.print)write(6,120)e2e,e3p,e3h,e3r,e3rt,e3
 120  format('E2,E3P,E3H,E3R,e3rt,E3:',6F14.10)
      if (iter.eq.1)then
      UN(4)=E3
      UN(5)=E3T
      WRITE(6,124)E3
      if(print)then
      write(6,*)'     UNLINKED DIAGRAMS'
      write(6,121)(e2t*e2e)
      write(6,122)(e2t*e3),(two*e3t*e2e)
      write(6,123)(two*e3t*e3),(e2tt*e2e*e2e)
      endif
 121  format(' O2*E2:',F15.10)
 122  format(' O2*E3,O3*E2:',2F15.10)
 123  format(' O3*E3,O2P*E2**2:',2F15.10)
 124  format('THIRD  ORDER ENERGY:',F15.10)
      WRITE(6,*)'*****************************************************'
      endif
      CALL ZEROMA(TI,1,NOU2)
      DO 161 I=1,NO
      DO 160 A=1,NU
      DO 160 B=1,NU
      DO 160 J=1,NO
      TI(A,B,J)=T2VO(I,A,B,J)
  160 CONTINUE
      IAS=I+NO
      WRITE(NT2T4,REC=IAS)TI
  161 CONTINUE
      call iexit(54)
      RETURN
      END
      SUBROUTINE T4SQUA(I,J,K,L,LS,NO,NU,EQ,EQS,EQD,EQ2,EQ3,EQQ,
     *TI,T2,VOE,VO,T4,V,OEH,OEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SDTQ123
      INTEGER A,B,C,D,E
      COMMON/ITRAT/ITER,MAXIT,ICCNV
      COMMON /NEWOPT/NOPT(6)
      COMMON/NEWTAPE/NT2T4,NT4INT
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      COMMON /RESLTS/CMP(30)
      DIMENSION T4(NU,NU,NU,NU),T2(NO,NU,NU,NO),VOE(NO,NU,NU,NO),
     *TI(NU,NU,NO),VO(NU,NU,NO,NO),V(1),OEH(NO),OEP(NU)
      DATA ZERO/0.0D+0/,TWO/2.0D+0/,THREE/3.0D+0/,FOUR/4.0D+0/,
     *     HALF/0.5D+0/,EIGHT/8.0D+0/,SIXTN/16.0D+0/,ONE/1.0D+0/
      call ienter(72)
      RIT=ZERO
      TCC=ZERO
      SDTQ123=NOPT(1).GE.3.AND.NOPT(1).LE.5
      no2u2=no*no*nu*nu
      nu4=nu*nu*nu*nu
      CALL RO2HPP(1,NO,NU,TI,VOE)
      DIJKL=OEH(I)+OEH(J)+OEH(K)+OEH(L)
      RITQ=RIT
      nou2=no*no*Nu
      DO 361 I1=1,NO
      IAS =I1+NO
      READ (NT2T4,REC=IAS)TI
      DO 360 J1=1,NO
      DO 360 A =1,NU
      DO 360 B= 1,NU
      VO(A,B,J1,I1)=TI(A,B,J1)
 360  CONTINUE
 361  CONTINUE
      IF(LS.EQ.0)GOTO 100
      X1=ZERO
      X2=ZERO
      X3=ZERO
      X4=ZERO
      X5=ZERO
      XDS=ZERO
      DO 50 A=1,NU
      DO 60 B=1,NU
      DO 70 C=1,NU
      IF (A.EQ.B.AND.A.EQ.C)GOTO 70
      DABC=OEP(A)+OEP(B)+OEP(C)
      DO 80 D=1,NU
      IF(A.EQ.B.AND.B.EQ.D.OR.A.EQ.C.AND.C.EQ.D.OR.
     *B.EQ.C.AND.C.EQ.D)GO TO 80
      DENOM=DIJKL-DABC-OEP(D)
      D4ABCD=(T2(I,A,B,J)*T2(K,C,D,L)+T2(I,A,C,K)*T2(J,B,D,L)
     *      + T2(I,A,D,L)*T2(J,B,C,K))*DENOM
      T4ABCD=T4(A,B,C,D)*DENOM
      X1=X1+T4ABCD*T4(A,B,C,D)*TWO/THREE
      X2=X2-T4ABCD*T4(A,B,D,C)*TWO
      X3=X3+T4ABCD*T4(A,D,B,C)*FOUR/THREE
      X4=X4+T4ABCD*T4(B,A,D,C)/TWO
      X5=X5-T4ABCD*T4(D,A,B,C)/TWO
      XDS=XDS+D4ABCD*(T4(A,B,C,D)*TWO/THREE
     *       -T4(A,B,D,C)*TWO+T4(A,D,B,C)*FOUR/THREE
     *       +T4(B,A,D,C)/TWO-T4(D,A,B,C)/TWO)
 80   CONTINUE
 70   CONTINUE
 60   CONTINUE
 50   CONTINUE
      XSUM=X1+X2+X3+X4+X5
      ET6=ET6+XSUM
      ET5=ET5+XDS
 900  FORMAT('E4SQ  FOR AL IJKL:E5Q,E6Q:',2F15.10)
      WRITE(6,900)ET5,ET6
      IF (I.GE.J.AND.J.GE.K.AND.K.GE.L) GO TO 100
      GO TO 200
 100  CONTINUE
      XDS0=ZERO
      XDS1=ZERO
      XDS2=ZERO
      XDS3=ZERO
      XC=ZERO
      XQQ=ZERO
      DIJ=OEH(I)+OEH(J)
      DIK=OEH(I)+OEH(K)
      DIL=OEH(I)+OEH(L)
      DJK=OEH(J)+OEH(K)
      DJL=OEH(J)+OEH(L)
      DKL=OEH(K)+OEH(L)
      IF(SDTQ123.and.nopt(6).ne.8) THEN
      IF(ITER.EQ.1.OR.(ITER.GT.1.AND.I.EQ.2.AND.J.EQ.2)) THEN
      iiii=nt4
      nt4=not4
      not4=iiii
      ENDIF
      IF(NOPT(4).EQ.ITER)THEN
         KK=2*(NOPT(4)/2)
      IF (KK.NE.NOPT(4).AND.I.EQ.2) THEN
      IIII=NT4
      NT4=NOT4
      NOT4=IIII
      ENDIF 
      ENDIF
      IAS =IT4(I,J,K,L)
      IF(ITER.EQ.1.AND.(L.LT.(NO-1)))THEN
      IIII=NT4
      NT4=NOT4
      NOT4=IIII
      ENDIF
      ENDIF
      DO 150 A=1,NU
      OEPA=OEP(A)
      DO 160 B=1,NU
      OEPB=OEP(B)
      DO 170 C=1,NU
      OEPC=OEP(C)
      IF (A.EQ.B.AND.B.EQ.C)GOTO 170
      DABC=OEPA+OEPB+OEPC
      DO 180 D=1,NU      
      OEPD=OEP(D)
      IF(A.EQ.B.AND.B.EQ.D.OR.A.EQ.C.AND.C.EQ.D.OR.
     *B.EQ.C.AND.C.EQ.D)GO TO 180
      DENOM=DIJKL-DABC-OEPD
      T4ABCD=T4(A,B,C,D)
      T4ABCD=T4ABCD*DENOM
      D1=T4(A,B,C,D)
      D2=T4(A,B,D,C)+T4(A,D,C,B)+T4(A,C,B,D)
     *  +T4(D,B,C,A)+T4(C,B,A,D)+T4(B,A,C,D)
      D3=T4(A,D,B,C)+T4(A,C,D,B)+T4(D,B,A,C)+T4(C,B,D,A)
     *  +T4(D,A,C,B)+T4(B,D,C,A)+T4(C,A,B,D)+T4(B,C,A,D)
      D4=T4(B,A,D,C)+T4(C,D,A,B)+T4(D,C,B,A)
      D5=T4(B,C,D,A)+T4(B,D,A,C)+T4(C,A,D,B)
     *  +T4(C,D,B,A)+T4(D,A,B,C)+T4(D,C,A,B)
C
      VIABJ=VOE(I,A,B,J)
      VIACK=VOE(I,A,C,K)
      VIADL=VOE(I,A,D,L)
      VJBCK=VOE(J,B,C,K)
      VJBDL=VOE(J,B,D,L)
      VKCDL=VOE(K,C,D,L)
C
      VIABJ=VIABJ/(DIJ-OEPA-OEPB)
      VIACK=VIACK/(DIK-OEPA-OEPC)
      VIADL=VIADL/(DIL-OEPA-OEPD)
      VJBCK=VJBCK/(DJK-OEPB-OEPC)
      VJBDL=VJBDL/(DJL-OEPB-OEPD)
      VKCDL=VKCDL/(DKL-OEPC-OEPD)
C
      D0ABCD=(VIABJ*VKCDL+VIACK*VJBDL+VIADL*VJBCK)*DENOM
C
      D1ABCD=(T2(I,A,B,J)*VKCDL+T2(I,A,C,K)*VJBDL
     *       +T2(I,A,D,L)*VJBCK+T2(J,B,C,K)*VIADL
     *       +T2(J,B,D,L)*VIACK+T2(K,C,D,L)*VIABJ)*DENOM*HALF
C
      D2ABCD=(T2(I,A,B,J)*T2(K,C,D,L)+T2(I,A,C,K)*T2(J,B,D,L)
     *      + T2(I,A,D,L)*T2(J,B,C,K))*DENOM
C
      D3ABCD=(VO(A,B,J,I)*VKCDL+VO(A,C,K,I)*VJBDL
     *       +VO(A,D,L,I)*VJBCK+VO(B,C,K,J)*VIADL
     *       +VO(B,D,L,J)*VIACK+VO(C,D,L,K)*VIABJ)*DENOM
C
      F=D1*SIXTN-EIGHT*D2+FOUR*(D3+D4)-TWO*D5
      XDS0=XDS0+F*D0ABCD
      XDS1=XDS1+F*D1ABCD
      XDS2=XDS2+F*D2ABCD
      XDS3=XDS3+F*D3ABCD
      XC  =XC  +F*T4ABCD
 180  CONTINUE
 170  CONTINUE
 160  CONTINUE
 150  CONTINUE
      CF=ONE
      IF(I.EQ.J.OR.J.EQ.K.OR.K.EQ.L) CF=HALF
      IF(I.EQ.J.AND.K.EQ.L)CF=HALF*HALF
      EQ =EQ +CF*XDS0
      EQS=EQS+CF*XDS1
      EQD=EQD+CF*XDS2
      EQ2=EQ2+CF*XC
      EQ3=EQ3+CF*XDS3
      EQQ=EQQ+CF*XQQ
      IF (J.EQ.NO.AND.L.EQ.(NO-1)) THEN
      CMP(6)=EQD
      WRITE(6,1010)EQ,(EQ+RITQ)
      WRITE(6,1020)EQS,(EQS+RITQ)
      WRITE(6,1030)EQD,(EQD+RITQ)
      WRITE(6,1040)EQ2,(RITQ+EQ+EQ2)
      WRITE(6,1050)(RITQ+EQS+EQ2),(RITQ+EQD+EQ2)
      WRITE(6,1060)EQ3
      WRITE(6,1070)EQQ
      CMP(7)=EQ2+EQ3
      ENDIF
 1010 FORMAT(1X,'     <WW|T4>  =',F15.10,4X,'E+    <WW|T4> =',F15.10)
 1020 FORMAT(1X,'    <WT2|T4>  =',F15.10,4X,'E+   <WT2|T4> =',F15.10)
 1030 FORMAT(1X,'   <T2T2|T4>  =',F15.10,4X,'E+  <T2T2|T4> =',F15.10)
 1040 FORMAT(1X,'      <T4|T4> =',F15.10,4X,'E+ <WW+T4|T4> =',F15.10)
 1050 FORMAT(1X,'E+<WT2+T4|T4> =',F15.10,4X,'E+<T2T2+T4|T4>=',F15.10)
 200  CONTINUE
 1060 FORMAT(1X,'  <WT2(2)|T4> =',F15.10)
 1070 FORMAT(1X,'   <T4(3)|T4> =',F15.10)
      call iexit(72)
      RETURN
      END
      SUBROUTINE TRANT3(V,NU,ID)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      DIMENSION V(NU,NU,NU)
      call ienter(45)
      GO TO (1,2,3,4,5),ID
    1 CONTINUE
      DO 100 B=1,NU
      DO 100 C=1,B
      DO 100 A=1,NU
      X=V(A,B,C)
      V(A,B,C)=V(A,C,B)
      V(A,C,B)=X
  100 CONTINUE
      GO TO 1000
    2 CONTINUE
      DO 200 C=1,NU
      DO 200 A=1,NU
      DO 200 B=1,A
      X=V(A,B,C)
      V(A,B,C)=V(B,A,C)
      V(B,A,C)=X
  200 CONTINUE
      GO TO 1000
    3 CONTINUE
      DO 300 B=1,NU
      DO 300 A=1,NU
      DO 300 C=1,A
      X=V(A,B,C)
      V(A,B,C)=V(C,B,A)
      V(C,B,A)=X
  300 CONTINUE
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
      call iexit(45)
      RETURN
      END
      SUBROUTINE TRANVT6(VT,NO,NU,ID)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      DIMENSION VT(NU,NU,NO,NO,NO,NO)
      GOTO (1,2,3,4,5,6)ID
 1    CONTINUE
      DO 10 A=1,NU
      DO 10 B=1,NU
      DO 10 I=1,NO
      DO 10 K=1,NO
      DO 10 L=1,K
      DO 10 J=1,L
      X=VT(A,B,I,J,K,L)
      VT(A,B,I,J,K,L)=VT(A,B,I,K,L,J)
      VT(A,B,I,K,L,J)=VT(A,B,I,L,J,K)
      VT(A,B,I,L,J,K)=X
      IF (K.EQ.L.OR.L.EQ.J)GOTO 10
      X=VT(A,B,I,K,J,L)
      VT(A,B,I,K,J,L)=VT(A,B,I,J,L,K)
      VT(A,B,I,J,L,K)=VT(A,B,I,L,K,J)
      VT(A,B,I,L,K,J)=X
 10   CONTINUE
      GOTO 1000
 2    CONTINUE
      DO 20 A=1,NU
      DO 20 B=1,NU
      DO 20 I=1,NO
      DO 20 J=1,NO
      DO 20 K=1,NO
      DO 20 L=1,J
      X=VT(A,B,I,J,K,L)
      VT(A,B,I,J,K,L)=VT(A,B,I,L,K,J)
      VT(A,B,I,L,K,J)=X
 20   CONTINUE
      GOTO 1000
 3    CONTINUE
      DO 30 A=1,NU
      DO 30 B=1,A
      DO 30 I=1,NO
      DO 30 J=1,NO
      DO 30 K=1,NO
      DO 30 L=1,NO
      X=VT(A,B,I,J,K,L)
      VT(A,B,I,J,K,L)=VT(B,A,I,J,K,L)
      VT(B,A,I,J,K,L)=X
 30   CONTINUE
      GOTO 1000
 4    CONTINUE
      DO 40 A=1,NU
      DO 40 B=1,NU
      DO 40 I=1,NO
      DO 40 J=1,NO
      DO 40 K=1,J
      DO 40 L=1,NO
      X=VT(A,B,I,J,K,L)
      VT(A,B,I,J,K,L)=VT(A,B,I,K,J,L)
      VT(A,B,I,K,J,L)=X
 40   CONTINUE
      GOTO 1000
 5    CONTINUE
      DO 50 A=1,NU
      DO 50 B=1,NU
      DO 50 I=1,NO
      DO 50 J=1,NO
      DO 50 K=1,NO
      DO 50 L=1,K
      X=VT(A,B,I,J,K,L)
      VT(A,B,I,J,K,L)=VT(A,B,I,J,L,K)
      VT(A,B,I,J,L,K)=X
 50   CONTINUE
      GOTO 1000
 6    CONTINUE
      DO 60 A=1,NU
      DO 60 B=1,NU
      DO 60 I=1,NO
      DO 60 J=1,NO
      DO 60 K=1,J
      DO 60 L=1,K
      X=VT(A,B,I,K,L,J)
      VT(A,B,I,K,L,J)=VT(A,B,I,J,K,L)
      VT(A,B,I,J,K,L)=VT(A,B,I,L,J,K)
      VT(A,B,I,L,J,K)=X
      IF (J.EQ.K.OR.K.EQ.L)GOTO 60
      X=VT(A,B,I,J,L,K)
      VT(A,B,I,J,L,K)=VT(A,B,I,K,J,L)
      VT(A,B,I,K,J,L)=VT(A,B,I,L,K,J)
      VT(A,B,I,L,K,J)=X
 60   CONTINUE
      GOTO 1000
 1000 CONTINUE
      RETURN
      END
      SUBROUTINE TRT3ALL(NO,NU,T3,IS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C
      DIMENSION T3(NO,NU,NU,NU,NO,NO)
      call ienter(51)
      DO 10 K=1,NO
      DO 10 J=1,NO
      DO 10 C=1,NU
      DO 10 A=1,NU
      DO 10 B=1,A
      DO 10 I=1,NO
      X=T3(I,A,B,C,J,K)
      T3(I,A,B,C,J,K)=T3(I,B,A,C,J,K)
      T3(I,B,A,C,J,K)=X
 10   CONTINUE
      call iexit(51)
      RETURN
      END
      SUBROUTINE MTRANREC(V,NU,ID)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      DIMENSION V(NU,NU,NU,NU)
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,
     *22,23),ID
    1 CONTINUE
      DO 100 K=1,NU
      DO 100 L=1,K
      DO 100 J=1,NU
      DO 100 I=1,NU
      X=V(I,J,K,L)
      V(I,J,K,L)=V(I,J,L,K)
      V(I,J,L,K)=X
  100 CONTINUE
      GO TO 1000
    2 CONTINUE
      DO 200 B=1,NU
      DO 200 C=1,B
      DO 200 D=1,C
      DO 200 A=1,NU
      X=V(A,D,B,C)
      V(A,D,B,C)=V(A,B,C,D)
      V(A,B,C,D)=V(A,C,D,B)
      V(A,C,D,B)=X
      IF(B.EQ.C.OR.C.EQ.D)GO TO 200
      X=V(A,B,D,C)
      V(A,B,D,C)=V(A,D,C,B)
      V(A,D,C,B)=V(A,C,B,D)
      V(A,C,B,D)=X
  200 CONTINUE
      GO TO 1000
    3 CONTINUE
      DO 300 B=1,NU
      DO 300 C=1,B
      DO 300 D=1,C
      DO 300 A=1,NU
      X=V(A,C,D,B)
      V(A,C,D,B)=V(A,B,C,D)
      V(A,B,C,D)=V(A,D,B,C)
      V(A,D,B,C)=X
      IF (B.EQ.C.OR.C.EQ.D)GO TO 300
      X=V(A,B,D,C)
      V(A,B,D,C)=V(A,C,B,D)
      V(A,C,B,D)=V(A,D,C,B)
      V(A,D,C,B)=X
  300 CONTINUE
      GO TO 1000
    4 CONTINUE
	call mtransnew(v,nu,3)
	call mtransnew(v,nu,8)
	goto 1000
      DO 400 A=1,NU
      DO 400 B=1,A
      DO 400 C=1,A
      LC=C
      IF (A.EQ.C)LC=B
      DO 400 D=1,LC
      X=V(D,C,A,B)
      V(D,C,A,B)=V(A,B,C,D)
      V(A,B,C,D)=V(C,D,B,A)
      V(C,D,B,A)=V(B,A,D,C)
      V(B,A,D,C)=X
      IF (C.EQ.D.OR.A.EQ.C.AND.B.EQ.D.OR.A.EQ.B)GO TO 400
      X=V(C,D,A,B)
      V(C,D,A,B)=V(A,B,D,C)
      V(A,B,D,C)=V(D,C,B,A)
      V(D,C,B,A)=V(B,A,C,D)
      V(B,A,C,D)=X
  400 CONTINUE
      GO TO 1000
    5 CONTINUE
      DO 500 A=1,NU
      DO 500 C=1,A
      DO 500 D=1,C
      DO 500 B=1,NU
      X=V(C,B,D,A)
      V(C,B,D,A)=V(A,B,C,D)
      V(A,B,C,D)=V(D,B,A,C)
      V(D,B,A,C)=X
      IF (A.EQ.C.OR.C.EQ.D)GO TO 500
      X=V(A,B,D,C)
      V(A,B,D,C)=V(C,B,A,D)
      V(C,B,A,D)=V(D,B,C,A)
      V(D,B,C,A)=X
  500 CONTINUE
      GO TO 1000
    6 CONTINUE
      DO 600 D=1,NU
      DO 600 C=1,NU
      DO 600 A=1,NU
      DO 600 B=1,A
      X=V(B,A,C,D)
      V(B,A,C,D)=V(A,B,C,D)
      V(A,B,C,D)=X
  600 CONTINUE
c      do 600 d=1,nu
c      do 600 c=1,nu
c      call transq(v(1,1,c,d),nu)
c 600  continue
      GO TO 1000
    7 CONTINUE
      DO 700 D=1,NU
      DO 700 B=1,NU
      DO 700 C=1,B
      DO 700 A=1,NU
      X=V(A,B,C,D)
      V(A,B,C,D)=V(A,C,B,D)
      V(A,C,B,D)=X
  700 CONTINUE
      GO TO 1000
    8 CONTINUE
      DO 800 D=1,NU
      DO 800 B=1,NU
      DO 800 A=1,NU
      DO 800 C=1,A
      X=V(A,B,C,D)
      V(A,B,C,D)=V(C,B,A,D)
      V(C,B,A,D)=X
  800 CONTINUE
      GO TO 1000
    9 CONTINUE
c      DO 900 A=1,NU
cv      DO 900 B=1,NU
c      DO 900 C=1,A
c      LIMD=NU
c      IF (A.EQ.C)LIMD=B
c      DO 900 D=1,LIMD
c      X=V(A,B,C,D)
c      V(A,B,C,D)=V(C,D,A,B)
c      V(C,D,A,B)=X
c  900 CONTINUE
      nu2=nu*nu
c      call transq1(v,nu2)
      GO TO 1000
 10   CONTINUE
      DO 950 A=1,NU
      DO 950 B=1,NU
      DO 950 C=1,NU
      DO 950 D=1,B
      X=V(A,B,C,D)
      V(A,B,C,D)=V(A,D,C,B)
      V(A,D,C,B)=X
 950  CONTINUE
      GOTO 1000
 11   CONTINUE
      DO 960 A=1,NU
      DO 960 B=1,NU
      DO 960 C=1,NU
      DO 960 D=1,A
      X=V(A,B,C,D)
      V(A,B,C,D)=V(D,B,C,A)
      V(D,B,C,A)=X
 960  CONTINUE
      GOTO 1000
 12   CONTINUE
      call tranmd(v,nu,nu,nu,nu,231)      
      GOTO 1000
 13   CONTINUE
      call tranmd(v,nu,nu,nu,nu,312)
      GOTO 1000
 14   CONTINUE
      DO 970 A=1,NU
      DO 970 B=1,NU
      DO 970 C=1,A
      DO 970 D=1,C
      X=V(D,B,A,C)
      V(D,B,A,C)=V(A,B,C,D)
      V(A,B,C,D)=V(C,B,D,A)
      V(C,B,D,A)=X
      IF (A.EQ.C.OR.C.EQ.D)GO TO 970
      X=V(D,B,C,A)
      V(D,B,C,A)=V(C,B,A,D)
      V(C,B,A,D)=V(A,B,D,C)
      V(A,B,D,C)=X
  970 CONTINUE
      GO TO 1000
 15   CONTINUE
      DO 980 A=1,NU
      DO 980 B=1,A
      DO 980 C=1,NU
      DO 980 D=1,B
      X=V(D,A,C,B)
      V(D,A,C,B)=V(A,B,C,D)
      V(A,B,C,D)=V(B,D,C,A)
      V(B,D,C,A)=X
      IF (A.EQ.B.OR.B.EQ.D)GO TO 980
      X=V(D,B,C,A)
      V(D,B,C,A)=V(B,A,C,D)
      V(B,A,C,D)=V(A,D,C,B)
      V(A,D,C,B)=X
 980  CONTINUE
      GO TO 1000
 16   CONTINUE
      DO 990 A=1,NU
      DO 990 B=1,A
      DO 990 C=1,NU
      DO 990 D=1,B
      X=V(B,D,C,A)
      V(B,D,C,A)=V(A,B,C,D)
      V(A,B,C,D)=V(D,A,C,B)
      V(D,A,C,B)=X
      IF (A.EQ.B.OR.B.EQ.D)GO TO 990
      X=V(A,D,C,B)
      V(A,D,C,B)=V(B,A,C,D)
      V(B,A,C,D)=V(D,B,C,A)
      V(D,B,C,A)=X
  990 CONTINUE
      GO TO 1000
 17   CONTINUE
      DO 991 A=1,NU
      DO 991 B=1,A
      DO 991 C=1,A
      LC=C
      IF (A.EQ.C)LC=B
      DO 991 D=1,LC
      X=V(A,B,C,D)
      V(A,B,C,D)=V(D,C,A,B)
      V(D,C,A,B)=V(B,A,D,C)
      V(B,A,D,C)=V(C,D,B,A)
      V(C,D,B,A)=X
      IF (C.EQ.D.OR.A.EQ.C.AND.B.EQ.D.OR.A.EQ.B)GO TO 991
      X=V(A,B,D,C)
      V(A,B,D,C)=V(C,D,A,B)
      V(C,D,A,B)=V(B,A,C,D)
      V(B,A,C,D)=V(D,C,B,A)
      V(D,C,B,A)=X
 991  CONTINUE
      GO TO 1000
 18   CONTINUE
      DO 992 A=1,NU
      DO 992 B=1,NU
      DO 992 D=1,A
      LC=NU
      IF (A.EQ.d)LC=B
      DO 992 c=1,LC
      X=V(A,B,C,D)
      V(A,B,C,D)=V(D,C,B,A)
      V(D,C,B,A)=X
 992  CONTINUE
      GO TO 1000
 19   CONTINUE
      DO 993 A=1,NU
c      write(6,*)'a=',a
      DO 993 B=1,a
      DO 993 C=1,nu
      DO 993 D=1,c
      X=V(A,B,C,D)
      V(A,B,C,D)=V(B,A,D,C)
      V(B,A,D,C)=X
      if(a.eq.b.or.c.eq.d)goto 9993
      x=v(a,b,d,c)
      v(a,b,d,c)=v(b,a,c,d)
      v(b,a,c,d)=x
 9993 continue
c      write(6,909)a,b,c,d,v(a,b,c,d),b,a,d,c,v(b,a,d,c)
 993  CONTINUE
 909  format('a,b,c,d,v:',4i3,i8,5x,4i3,i8)
      GO TO 1000
 20   CONTINUE
      DO 994 A=1,NU
      DO 994 B=1,A
      DO 994 C=1,A
      LD=B
      IF (A.EQ.B)Ld=C
      DO 994 D=1,LD
      X=V(A,B,C,D)
      V(A,B,C,D)=v(B,C,D,A)
      V(B,C,D,A)=V(C,D,A,B)
      V(C,D,A,B)=V(D,A,B,C)
      V(D,A,B,C)=X
      IF (B.EQ.D.OR.A.EQ.B.AND.C.EQ.D.OR.A.EQ.C)GO TO 994
      X=V(A,D,C,B)
      V(A,D,C,B)=v(D,C,B,A)
      V(D,C,B,A)=v(C,B,A,D)
      V(C,B,A,D)=v(B,A,D,C)
      V(B,A,D,C)=X
 994  CONTINUE
      GO TO 1000
 21   CONTINUE
      DO 995 A=1,NU
      DO 995 B=1,A
      DO 995 C=1,A
      LD=B
      IF (A.EQ.B)LD=C
      DO 995 D=1,LD
      X=V(A,B,C,D)
      V(A,B,C,D)=V(D,A,B,C)
      V(D,A,B,C)=V(C,D,A,B)
      V(C,D,A,B)=V(B,C,D,A)
      V(B,C,D,A)=X
      IF (B.EQ.D.OR.A.EQ.B.AND.C.EQ.D.OR.A.EQ.C)GO TO 995
      X=V(A,D,C,B)
      V(A,D,C,B)=V(B,A,D,C)
      V(B,A,D,C)=V(C,B,A,D)
      V(C,B,A,D)=V(D,C,B,A)
      V(D,C,B,A)=X
 995  CONTINUE
      GO TO 1000
 22   CONTINUE
      DO 996 A=1,NU
      DO 996 B=1,A
      DO 996 D=1,A
      LC=B
      IF (A.EQ.B)LC=D
      DO 996 C=1,LC
      X=V(A,B,C,D)
      V(A,B,C,D)=V(B,D,A,C)
      V(B,D,A,C)=V(D,C,B,A)
      V(D,C,B,A)=V(C,A,D,B)
      V(C,A,D,B)=X
      IF (B.EQ.C.OR.A.EQ.B.AND.C.EQ.D.OR.A.EQ.D)GO TO 996
      X=V(A,C,B,D)
      V(A,C,B,D)=V(C,D,A,B)
      V(C,D,A,B)=V(D,B,C,A)
      V(D,B,C,A)=V(B,A,D,C)
      V(B,A,D,C)=X
 996  CONTINUE
      GO TO 1000
 23   CONTINUE
      DO 997 A=1,NU
      DO 997 B=1,A
      DO 997 D=1,A
      LC=B
      IF (A.EQ.B)LC=D
      DO 997 C=1,LC
      X=V(A,B,C,D)
      V(A,B,C,D)=V(C,A,D,B)
      V(C,A,D,B)=V(D,C,B,A)
      V(D,C,B,A)=V(B,D,A,C)
      V(B,D,A,C)=X
      IF (B.EQ.C.OR.A.EQ.B.AND.C.EQ.D.OR.A.EQ.D)GO TO 997
      X=V(A,C,B,D)
      V(A,C,B,D)=V(B,A,D,C)
      V(B,A,D,C)=V(D,B,C,A)
      V(D,B,C,A)=V(C,D,A,B)
      V(C,D,A,B)=X
 997  CONTINUE
      GO TO 1000
 1000 continue
c      call iexit(2)
      RETURN
      END
      SUBROUTINE MTRA1(V,NU,IDENT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      LOGICAL AEB,AEC,AED,BEC,BED,CED,ACD,ABD
      DIMENSION V(1)
      NU2=NU*NU
      NU3=NU2*NU
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,
     *22,23),IDENT
    1 CONTINUE
      DO 100 J=1,NU
         j1=(j-1)*nu
         DO 101 K=1,NU
            k0=k-1
            k2=k0*nu2
            k3=k0*nu3
            DO 102 L=1,K
               l0=l-1
               jkl=j1+k2+l0*nu3
               jlk=j1+k3+l0*nu2
               DO 103 I=1,NU
                  ijkl=jkl+i
                  ijlk=jlk+i
                  X=V(ijkl)
                  V(ijkl)=V(ijlk)
                  V(ijlk)=X
 103           CONTINUE
 102        CONTINUE
 101     CONTINUE
 100  CONTINUE
      GO TO 1000
 2    CONTINUE
      DO 120 B=1,NU
         ib=b-1
         ib1=ib*nu
         ib2=ib*nu2
         ib3=ib*nu3
         DO 121 C=1,B
            BEC=B.EQ.C
            ic=c-1
            ic1=ic*nu
            ic2=ic*nu2
            ic3=ic*nu3
            DO 122 D=1,C
               CED=C.EQ.D
               id=d-1
               id1=id*nu
               id2=id*nu2
               id3=id*nu3
               idbc=id1+ib2+ic3
               ibcd=ib1+ic2+id3
               icdb=ic1+id2+ib3
               if(bec.OR.ced)GOTO 124
               ibdc=ib1+id2+ic3
               idcb=id1+ic2+ib3
               icbd=ic1+ib2+id3
 124           CONTINUE
               DO 123 A=1,NU
                  iadbc=idbc+a
                  iabcd=ibcd+a
                  iacdb=icdb+a
                  X=V(iadbc)
                  V(iadbc)=V(iabcd)
                  V(iabcd)=V(iacdb)
                  V(iacdb)=X
                  IF(BEC.OR.CED)GO TO 123
                  iabdc=a+ibdc
                  iadcb=a+idcb
                  iacbd=a+icbd
                  X=V(iabdc)
                  V(iabdc)=V(iadcb)
                  V(iadcb)=V(iacbd)
                  V(iacbd)=X
 123           CONTINUE
 122        CONTINUE
 121     CONTINUE
 120  CONTINUE
      GO TO 1000
    3 CONTINUE
      DO 300 B=1,NU
         ib=b-1
         ib1=ib*nu
         ib2=ib*nu2
         ib3=ib*nu3
         DO 301 C=1,B
            ic=c-1
            ic1=ic*nu
            ic2=ic*nu2
            ic3=ic*nu3
            DO 302 D=1,C
               id=d-1
               id1=id*nu
               id2=id*nu2
               id3=id*nu3
               icdb=ic1+id2+ib3
               ibcd=ib1+ic2+id3
               idbc=id1+ib2+ic3
               if(b.ne.c.and.c.ne.d)then
                  ibdc=ib1+id2+ic3
                  icbd=ic1+ib2+id3
                  idcb=id1+ic2+ib3
               endif
               DO 303 A=1,NU
                  iacdb=a+icdb
                  iabcd=a+ibcd
                  iadbc=a+idbc
                  X=V(iacdb)
                  V(IACDB)=V(IABCD)
                  V(IABCD)=V(IADBC)
                  V(IADBC)=X
                  IF (B.EQ.C.OR.C.EQ.D)GO TO 303
                  iabdc=a+ibdc
                  iacbd=a+icbd
                  iadcb=a+idcb
                  X=V(IABDC)
                  V(IABDC)=V(IACBD)
                  V(IACBD)=V(IADCB)
                  V(IADCB)=X
 303           CONTINUE
 302        CONTINUE
 301     CONTINUE
 300  CONTINUE
      GO TO 1000
    4 CONTINUE
      DO 400 A=1,NU
         ia=a-1
         ia1=ia*nu
         ia2=ia*nu2
         ia3=ia*nu3
         DO 401 B=1,A
            aeb=a.eq.b
            ib=b-1
            ib1=ib*nu
            ib2=ib*nu2
            ib3=ib*nu3
            DO 402 C=1,A
               aec=a.eq.c
               ic=c-1
               ic1=ic*nu
               ic2=ic*nu2
               ic3=ic*nu3
               icab=ic1+ia2+ib3
               iabc=a  +ib1+ic2
               icba=c  +ib2+ia3
               ibac=b  +ia1+ic3
               LC=C
               IF (AEC)LC=B
               id1=-nu
               id2=-nu2
               id3=-nu3
               DO 403 D=1,LC
                  ced=c.eq.d
                  bed=b.eq.d
c                  id=d-1
c                  id1=id*nu
c                  id2=id*nu2
c                  id3=id*nu3
                  id1=id1+nu
                  id2=id2+nu2
                  id3=id3+nu3
                  idcab=d   +icab
                  iabcd=iabc+id3
                  icdba=icba+id1
                  ibadc=ibac+id2
                  X=V(IDCAB)
                  V(IDCAB)=V(IABCD)
                  V(IABCD)=V(ICDBA)
                  V(ICDBA)=V(IBADC)
                  V(IBADC)=X
                  IF (CED.OR.AEC.AND.BED.OR.AEB)GO TO 403
                  icdab=c+id1+ia2+ib3
                  iabdc=a+ib1+id2+ic3
                  idcba=d+ic1+ib2+ia3
                  ibacd=b+ia1+ic2+id3
                  X=V(ICDAB)
                  V(ICDAB)=V(IABDC)
                  V(IABDC)=V(IDCBA)
                  V(IDCBA)=V(IBACD)
                  V(IBACD)=X
 403           CONTINUE
 402        CONTINUE
 401     CONTINUE
 400  CONTINUE
      GO TO 1000
    5 CONTINUE
      DO 500 A=1,NU
         ia=a-1
         ia2=ia*nu2
         ia3=ia*nu3
         DO 501 C=1,A
            AEC=A.EQ.C
            ic=c-1
            ic2=ic*nu2
            ic3=ic*nu3
            DO 502 D=1,C
               ACD=AEC.OR.C.EQ.D
               id=d-1
               id2=id*nu2
               id3=id*nu3
               icda=c+id2+ia3
               iacd=a+ic2+id3
               idac=d+ia2+ic3
               IF(ACD)GOTO 504
               iadc=a+id2+ic3
               icad=c+ia2+id3
               idca=d+ic2+ia3
 504           CONTINUE
               ib1=-nu
              DO 503 B=1,NU
                  ib1=ib1+nu
c                  ib=b-1
c                  ib1=ib*nu
                  icbda=ib1+icda
                  iabcd=ib1+iacd
                  idbac=ib1+idac
                  X=V(ICBDA)
                  V(ICBDA)=V(IABCD)
                  V(IABCD)=V(IDBAC)
                  V(IDBAC)=X
                  IF (ACD)GO TO 503
                  iabdc=ib1+iadc
                  icbad=ib1+icad
                  idbca=ib1+idca
                  X=V(IABDC)
                  V(IABDC)=V(ICBAD)
                  V(ICBAD)=V(IDBCA)
                  V(IDBCA)=X
 503           CONTINUE
 502        CONTINUE
 501     CONTINUE
 500  CONTINUE
      GO TO 1000
 6    CONTINUE
      DO 600 D=1,NU
         id=d-1
         id3=id*nu3
         DO 601 C=1,NU
            ic=c-1
            ic2=ic*nu2
            DO 602 A=1,NU
               ia=a-1
               ia1=ia*nu
               iacd =ia1+ic2+id3
               iacd_=a+ic2+id3
               ib1=-nu
               ib2=-nu2
               ib3=-nu3
               DO 603 B=1,A
                  ib1=ib1+nu
                  ib2=ib2+nu2
                  ib3=ib3+nu3
c                  ib=b-1
c                  ib1=ib*nu
c                  ib2=ib*nu2
c                  ib3=ib*nu3
                  ibacd=b  +iacd
                  iabcd=ib1+iacd_
                  X=V(IBACD)
                  V(IBACD)=V(IABCD)
                  V(IABCD)=X
 603           CONTINUE
 602        CONTINUE
 601     CONTINUE
 600  CONTINUE
      GO TO 1000
    7 CONTINUE
      DO 700 D=1,NU
         id=d-1
         id3=id*nu3
         DO 701 B=1,NU
            ib=b-1
            ib1=ib*nu
            ib2=ib*nu2
            DO 702 C=1,B
               ic=c-1
               ic1=ic*nu
               ic2=ic*nu2
               ibcd=ib1+ic2+id3
               icbd=ic1+ib2+id3
               DO 703 A=1,NU
                  iabcd=a+ibcd
                  iacbd=a+icbd
                  X=V(IABCD)
                  V(IABCD)=V(IACBD)
                  V(IACBD)=X
 703           CONTINUE
 702        CONTINUE
 701     CONTINUE
 700  CONTINUE
      GO TO 1000
 8    CONTINUE
      DO 800 D=1,NU
         id=d-1
         id3=id*nu3
         DO 801 B=1,NU
            ib=b-1
            ib1=ib*nu
            DO 802 A=1,NU
               ia=a-1
               ia2=ia*nu2
               iabd=a  +ib1+id3
               ibad=ib1+ia2+id3
               ic2=-nu2
               DO 803 C=1,A
c                  ic=c-1
c                  ic2=ic*nu2
                  ic2=ic2+nu2
                  iabcd=iabd+ic2
                  icbad=c   +ibad
                  X=V(IABCD)
                  V(IABCD)=V(ICBAD)
                  V(ICBAD)=X
 803           CONTINUE
 802        CONTINUE
 801     CONTINUE
 800  CONTINUE
      GO TO 1000
    9 CONTINUE
      DO 900 A=1,NU2
         ia=a-1
         ia2=ia*nu2
         ib2=-nu2
         DO 901 B=1,A
            ib2=ib2+nu2
c            ib=b-1
c            ib2=ib*nu2
            iab=a+ib2
            iba=b+ia2
            X=V(IAB)
            V(IAB)=V(IBA)
            V(IBA)=X
 901     CONTINUE
 900  CONTINUE
      GO TO 1000
 10   CONTINUE
      DO 950 B=1,NU
         ib=b-1
         ib1=ib*nu
         ib3=ib*nu3
         DO 951 D=1,B
            id=d-1
            id1=id*nu
            id3=id*nu3
            DO 952 C=1,NU
               ic=c-1
               ic2=ic*nu2
               ibcd=ib1+ic2+id3
               idcb=id1+ic2+ib3
               DO 953 A=1,NU
                  iabcd=a+ibcd
                  iadcb=a+idcb
                  X=V(IABCD)
                  V(IABCD)=V(IADCB)
                  V(IADCB)=X
 953           CONTINUE
 952        CONTINUE
 951     CONTINUE
 950  CONTINUE
      GOTO 1000
 11   CONTINUE
      DO 960 A=1,NU
         ia=a-1
         ia3=ia*nu3
         DO 961 D=1,A
            id=d-1
            id3=id*nu3
            ic2=-nu2
            iad=a+id3
            ida=d+ia3
            DO 962 C=1,NU
c               ic2=n2(c)
               ic2=ic2+nu2
c               ic=c-1
c               ic2=ic*nu2
               iacd=ic2+iad
               idca=ic2+ida
               ib1=-nu
               DO 963 B=1,NU
                  ib1=ib1+nu
c                  ib1=n1(b)
                  iabcd=ib1+iacd
                  idbca=ib1+idca
                  X=V(IABCD)
                  V(IABCD)=V(IDBCA)
                  V(IDBCA)=X
 963           CONTINUE
 962        CONTINUE
 961     CONTINUE
 960  CONTINUE
      GOTO 1000
 12   CONTINUE
      call tranmd(v,nu,nu,nu,nu,231)      
      GOTO 1000
 13   CONTINUE
      call tranmd(v,nu,nu,nu,nu,312)
      GOTO 1000
 14   CONTINUE
      DO 970 A=1,NU
         ia=a-1
         ia2=ia*nu2
         ia3=ia*nu3
         DO 971 C=1,A
            AEC=A.EQ.C
            ic=c-1
            ic2=ic*nu2
            ic3=ic*nu3
            DO 972 D=1,C
               id=d-1
               id2=id*nu2
               id3=id*nu3
               ACD=AEC.OR.C.EQ.D
               idac=d+ia2+ic3
               iacd=a+ic2+id3
               icda=c+id2+ia3
               IF (ACD)GOTO 974
               idca=d+ic2+ia3
               icad=c+ia2+id3
               iadc=a+id2+ic3
 974           CONTINUE
               ib1=-nu
               DO 973 B=1,NU
                  ib1=ib1+nu
c                  ib=b-1
c                  ib1=ib*nu
                  idbac=ib1+idac
                  iabcd=ib1+iacd
                  icbda=ib1+icda
                  X=V(IDBAC)
                  V(IDBAC)=V(IABCD)
                  V(IABCD)=V(ICBDA)
                  V(ICBDA)=X
                  IF (ACD)GO TO 973
                  idbca=ib1+idca
                  icbad=ib1+icad
                  iabdc=ib1+iadc
                  X=V(IDBCA)
                  V(IDBCA)=V(ICBAD)
                  V(ICBAD)=V(IABDC)
                  V(IABDC)=X
 973           CONTINUE
 972        CONTINUE
 971     CONTINUE
 970  CONTINUE
      GO TO 1000
 15   CONTINUE
      DO 980 A=1,NU
         ia=a-1
         ia1=ia*nu
         ia2=ia*nu2
         ia3=ia*nu3
         DO 981 B=1,A
            ib=b-1
            ib1=ib*nu
            ib2=ib*nu2
            ib3=ib*nu3
            AEB=A.EQ.B
            DO 982 D=1,B
               id=d-1
               id1=id*nu
               id2=id*nu2
               id3=id*nu3
               idab=d+ia1+ib3
               iabd=a+ib1+id3
               ibda=b+id1+ia3
               ABD=AEB.OR.B.EQ.D
               IF(ABD)GOTO 984
               idba=d+ib1+ia3
               ibad=b+ia1+id3
               iadb=a+id1+ib3
 984           CONTINUE
               ic2=-nu2
               DO 983 C=1,NU
                  ic2=ic2+nu2
c                  ic=c-1
c                  ic2=ic*nu2
                  idacb=ic2+idab
                  iabcd=ic2+iabd
                  ibdca=ic2+ibda
                  X=V(IDACB)
                  V(IDACB)=V(IABCD)
                  V(IABCD)=V(IBDCA)
                  V(IBDCA)=X
                  IF (ABD)GO TO 983
                  idbca=ic2+idba
                  ibacd=ic2+ibad
                  iadcb=ic2+iadb
                  X=V(IDBCA)
                  V(IDBCA)=V(IBACD)
                  V(IBACD)=V(IADCB)
                  V(IADCB)=X
 983           CONTINUE
 982        CONTINUE
 981     CONTINUE
 980  CONTINUE
      GO TO 1000
 16   CONTINUE
      DO 990 A=1,NU
         ia=a-1
         ia1=ia*nu
         ia2=ia*nu2
         ia3=ia*nu3
         DO 991 B=1,A
            AEB=A.EQ.B
            ib=b-1
            ib1=ib*nu
            ib2=ib*nu2
            ib3=ib*nu3
            DO 992 D=1,B
               ABD=AEB.OR.B.EQ.D
               id=d-1
               id1=id*nu
               id2=id*nu2
               id3=id*nu3
               ibda=b+id1+ia3
               iabd=a+ib1+id3
               idab=d+ia1+ib3
               IF (ABD)GO TO 994
               iadb=a+id1+ib3
               ibad=b+ia1+id3
               idba=d+ib1+ia3
 994           CONTINUE
               ic2=-nu2
               DO 993 C=1,NU
                  ic2=ic2+nu2
c                  ic=c-1
c                  ic2=ic*nu2
                  ibdca=ic2+ibda
                  iabcd=ic2+iabd
                  idacb=ic2+idab
                  X=V(IBDCA)
                  V(IBDCA)=V(IABCD)
                  V(IABCD)=V(IDACB)
                  V(IDACB)=X
                  IF (ABD)GO TO 993
                  iadcb=ic2+iadb
                  ibacd=ic2+ibad
                  idbca=ic2+idba
                  X=V(IADCB)
                  V(IADCB)=V(IBACD)
                  V(IBACD)=V(IDBCA)
                  V(IDBCA)=X
 993           CONTINUE
 992        CONTINUE
 991     CONTINUE
 990  CONTINUE
      GO TO 1000
 17   CONTINUE
      DO 170 A=1,NU
         ia=a-1
         ia1=ia*nu
         ia2=ia*nu2
         ia3=ia*nu3
         DO 171 B=1,A
            AEB=A.EQ.B
            ib=b-1
            ib1=ib*nu
            ib2=ib*nu2
            ib3=ib*nu3
            DO 172 C=1,A
               AEC=A.EQ.C
               ic=c-1
               ic1=ic*nu
               ic2=ic*nu2
               ic3=ic*nu3
               LD=C
               IF (AEC)LD=B
               iabc=a  +ib1+ic2
               icab=ic1+ia2+ib3
               ibac=b  +ia1+ic3
               icba=c  +ib2+ia3
               IF (AEB)GO TO 174
               iabc_=a  +ib1+ic3
               icab_=c  +ia2+ib3
               ibac_=b  +ia1+ic2
               icba_=ic1+ib2+ia3
 174           CONTINUE
               id1=-nu
               id2=-nu2
               id3=-nu3
               DO 173 D=1,LD
                  CED=C.EQ.D
                  BED=B.EQ.D
                  id1=id1+nu
                  id2=id2+nu2
                  id3=id3+nu3
c                  id=d-1
c                  id1=id*nu
c                  id2=id*nu2
c                  id3=id*nu3
                  iabcd=id3+iabc
                  idcab=d  +icab
                  ibadc=id2+ibac
                  icdba=id1+icba
                  X=V(IABCD)
                  V(IABCD)=V(IDCAB)
                  V(IDCAB)=V(IBADC)
                  V(IBADC)=V(ICDBA)
                  V(ICDBA)=X
                  IF (CED.OR.AEC.AND.BED.OR.AEB)GO TO 173
                  iabdc=id2+iabc_
                  icdab=id1+icab_
                  ibacd=id3+ibac_
                  idcba=d  +icba_
                  X=V(IABDC)
                  V(IABDC)=V(ICDAB)
                  V(ICDAB)=V(IBACD)
                  V(IBACD)=V(IDCBA)
                  V(IDCBA)=X
 173           CONTINUE
 172        CONTINUE
 171     CONTINUE
 170  CONTINUE
      GO TO 1000
 18   CONTINUE
      DO 180 A=1,NU
         ia=a-1
         ia3=ia*nu3
         DO 181 B=1,NU
            ib=b-1
            ib1=ib*nu
            ib2=ib*nu2
            DO 182 D=1,A
               id=d-1
               id3=id*nu3
               iabd=a+ib1+id3
               idba=d+ib2+ia3
               LC=NU
               IF (A.EQ.D)LC=B
               ic1=-nu
               ic2=-nu2
               DO 183 c=1,LC
                  ic1=ic1+nu
                  ic2=ic2+nu2
c                  ic=c-1
c                  ic1=ic*nu
c                  ic2=ic*nu2
                  iabcd=ic2+iabd
                  idcba=ic1+idba
                  X=V(IABCD)
                  V(IABCD)=V(IDCBA)
                  V(IDCBA)=X
 183           CONTINUE
 182        CONTINUE
 181     CONTINUE
 180  CONTINUE
      GO TO 1000
 19   CONTINUE
      DO 190 C=1,nu
         ic=c-1
         ic2=ic*nu2
         ic3=ic*nu3
         DO 191 D=1,c
            CED=C.EQ.D
            id=d-1
            id2=id*nu2
            id3=id*nu3
            DO 192 A=1,NU
               ia=a-1
               ia1=ia*nu
               iacd=a  +ic2+id3
               iadc=ia1+id2+ic3
               if(CED)goto 194
               iadc_=a  +id2+ic3
               iacd_=ia1+ic2+id3
 194           CONTINUE
               ib1=-nu
               DO 193 B=1,a
                  AEB=A.EQ.B
c                  ib=b-1
                  ib1=ib1+nu
c                  ib1=ib*nu
                  iabcd=ib1+iacd
                  ibadc=b  +iadc
                  X=V(IABCD)
                  V(IABCD)=V(IBADC)
                  V(IBADC)=X
                  if(AEB.OR.CED)goto 193
                  iabdc=ib1+iadc_
                  ibacd=b  +iacd_
                  x=v(iabdc)
                  v(iabdc)=v(ibacd)
                  v(ibacd)=x
 193           CONTINUE
 192        CONTINUE
 191     CONTINUE
 190  CONTINUE
      GO TO 1000
 20   CONTINUE
      DO 200 A=1,NU
         ia=a-1
         ia1=ia*nu
         ia2=ia*nu2
         ia3=ia*nu3
         DO 201 B=1,A
            AEB=A.EQ.B
            ib=b-1
            ib1=ib*nu
            ib2=ib*nu2
            ib3=ib*nu3
            DO 202 C=1,A
               AEC=A.EQ.C
               ic=c-1
               ic1=ic*nu
               ic2=ic*nu2
               ic3=ic*nu3
               iabc =a  +ib1+ic2
               ibca =b  +ic1+ia3
               icab =c  +ia2+ib3
               iabc_=ia1+ib2+ic3
               LD=B
               IF (A.EQ.B)Ld=C
               IF(AEC)GOTO 204
               iacb=a  +ic2+ib3
               icba=ic1+ib2+ia3
               icba_=c +ib1+ia2
               ibac =b +ia1+ic3
 204           CONTINUE
               id1=-nu
               id2=-nu2
               id3=-nu3
               DO 203 D=1,LD
                  BED=B.EQ.D
                  CED=C.EQ.D
                  id1=id1+nu
                  id2=id2+nu2
                  id3=id3+nu3
c                  id=d-1
c                  id1=id*nu
c                  id2=id*nu2
c                  id3=id*nu3
                  iabcd=id3+iabc
                  ibcda=id2+ibca
                  icdab=id1+icab
                  idabc=d  +iabc_
                  X=V(IABCD)
                  V(IABCD)=v(IBCDA)
                  V(IBCDA)=V(ICDAB)
                  V(ICDAB)=V(IDABC)
                  V(IDABC)=X
                  IF (BED.OR.AEB.AND.CED.OR.AEC)GO TO 203
                  iadcb=id1+iacb
                  idcba=d  +icba
                  icbad=id3+icba_
                  ibadc=id2+ibac
                  X=V(IADCB)
                  V(IADCB)=v(IDCBA)
                  V(IDCBA)=v(ICBAD)
                  V(ICBAD)=v(IBADC)
                  V(IBADC)=X
 203           CONTINUE
 202        CONTINUE
 201     CONTINUE
 200  CONTINUE
      GO TO 1000
 21   CONTINUE
      DO 210 A=1,NU
         ia=a-1
         ia1=ia*nu
         ia2=ia*nu2
         ia3=ia*nu3
         DO 211 B=1,A
            AEB=A.EQ.B
            ib=b-1
            ib1=ib*nu
            ib2=ib*nu2
            ib3=ib*nu3
            DO 212 C=1,A
               AEC=A.EQ.C
               ic=c-1
               ic1=ic*nu
               ic2=ic*nu2
               ic3=ic*nu3
               iabc =a  +ib1+ic2
               iabc_=ia1+ib2+ic3
               icab =c  +ia2+ib3
               ibca =b  +ic1+ia3
               IF (AEC)GO TO 214
               iacb =a  +ic2+ib3
               ibac =b  +ia1+ic3
               icba =c  +ib1+ia2
               icba_=ic1+ib2+ia3
 214           CONTINUE
               LD=B
               IF (AEB)LD=C
               id1=-nu
               id2=-nu2
               id3=-nu3
               DO 213 D=1,LD
                  BED=B.EQ.D
                  CED=C.EQ.D
                  id1=id1+nu
                  id2=id2+nu2
                  id3=id3+nu3
c                  id=d-1
c                  id1=id*nu
c                  id2=id*nu2
c                  id3=id*nu3
                  iabcd=id3+iabc
                  idabc=d  +iabc_
                  icdab=id1+icab
                  ibcda=id2+ibca
                  X=V(IABCD)
                  V(IABCD)=V(IDABC)
                  V(IDABC)=V(ICDAB)
                  V(ICDAB)=V(IBCDA)
                  V(IBCDA)=X
                  IF (BED.OR.AEB.AND.CED.OR.AEC)GO TO 213
                  iadcb=id1+iacb
                  ibadc=id2+ibac
                  icbad=id3+icba
                  idcba=d  +icba_
                  X=V(IADCB)
                  V(IADCB)=V(IBADC)
                  V(IBADC)=V(ICBAD)
                  V(ICBAD)=V(IDCBA)
                  V(IDCBA)=X
 213           CONTINUE
 212        CONTINUE
 211     CONTINUE
 210  CONTINUE
      GO TO 1000
 22   CONTINUE
      DO 220 A=1,NU
         ia=a-1
         ia1=ia*nu
         ia2=ia*nu2
         ia3=ia*nu3
         DO 221 B=1,A
            AEB=A.EQ.B
            ib=b-1
            ib1=ib*nu
            ib2=ib*nu2
            ib3=ib*nu3
            DO 222 D=1,A
               AED=A.EQ.D
               id=d-1
               id1=id*nu
               id2=id*nu2
               id3=id*nu3
               iabd=a  +ib1+id3
               ibda=b  +id1+ia2
               idba=d  +ib2+ia3
               iadb=ia1+id2+ib3
               IF (AED)GO TO 224
               iabd_=a  +ib2+id3
               idab =id1+ia2+ib3
               idba_=d  +ib1+ia3
               ibad =b  +ia1+id2
 224           CONTINUE
               LC=B
               IF (AEB)LC=D
               ic1=-nu
               ic2=-nu2
               ic3=-nu3
               DO 223 C=1,LC
                  BEC=B.EQ.C
                  CED=C.EQ.D
                  ic1=ic1+nu
                  ic2=ic2+nu2
                  ic3=ic3+nu3
c                  ic=c-1
c                  ic1=ic*nu
c                  ic2=ic*nu2
c                  ic3=ic*nu3
                  iabcd=ic2+iabd
                  ibdac=ic3+ibda
                  idcba=ic1+idba
                  icadb=c  +iadb
                  X=V(IABCD)
                  V(IABCD)=V(IBDAC)
                  V(IBDAC)=V(IDCBA)
                  V(IDCBA)=V(ICADB)
                  V(ICADB)=X
                  IF (BEC.OR.AEB.AND.CED.OR.AED)GO TO 223
                  iacbd=ic1+iabd_
                  icdab=c  +idab
                  idbca=ic2+idba_
                  ibadc=ic3+ibad
                  X=V(IACBD)
                  V(IACBD)=V(ICDAB)
                  V(ICDAB)=V(IDBCA)
                  V(IDBCA)=V(IBADC)
                  V(IBADC)=X
 223           CONTINUE
 222        CONTINUE
 221     CONTINUE
 220  CONTINUE
      GO TO 1000
 23   CONTINUE
      DO 230 A=1,NU
         ia=a-1
         ia1=ia*nu
         ia2=ia*nu2
         ia3=ia*nu3
         DO 231 B=1,A
            AEB=A.EQ.B
            ib=b-1
            ib1=ib*nu
            ib2=ib*nu2
            ib3=ib*nu3
            DO 232 D=1,A
               AED=A.EQ.D
               id=d-1
               id1=id*nu
               id2=id*nu2
               id3=id*nu3
               iabd=a  +ib1+id3
               iadb=ia1+id2+ib3
               idba=d  +ib2+ia3
               ibda=b  +id1+ia2
               IF (AED)GO TO 234
               iabd_=a  +ib2+id3
               ibad =b  +ia1+id2
               idba_=d  +ib1+ia3
               idab =id1+ia2+ib3
 234           CONTINUE
               LC=B
               IF (AEB)LC=D
               ic1=-nu
               ic2=-nu2
               ic3=-nu3
               DO 233 C=1,LC
                  BEC=B.EQ.C
                  CED=C.EQ.D
                  ic1=ic1+nu
                  ic2=ic2+nu2
                  ic3=ic3+nu3
c                  ic=c-1
c                  ic1=ic*nu
c                  ic2=ic*nu2
c                  ic3=ic*nu3
                  iabcd=ic2+iabd
                  icadb=c  +iadb
                  idcba=ic1+idba
                  ibdac=ic3+ibda
                  X=V(IABCD)
                  V(IABCD)=V(ICADB)
                  V(ICADB)=V(IDCBA)
                  V(IDCBA)=V(IBDAC)
                  V(IBDAC)=X
                  IF (BEC.OR.AEB.AND.CED.OR.AED)GO TO 233
                  iacbd=ic1+iabd_
                  ibadc=ic3+ibad
                  idbca=ic2+idba_
                  icdab=c  +idab
                  X=V(IACBD)
                  V(IACBD)=V(IBADC)
                  V(IBADC)=V(IDBCA)
                  V(IDBCA)=V(ICDAB)
                  V(ICDAB)=X
 233           CONTINUE
 232        CONTINUE
 231     CONTINUE
 230  CONTINUE
      GO TO 1000
 1000 continue
c      call iexit(2)
      RETURN
      END
      SUBROUTINE TRANMD(A,N1,N2,N3,N4,IJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1)
      call ienter(43)
      n12=n1*n2
      n123=n12*n3
      IF(IJ.EQ.12) GO TO 12
      IF(IJ.EQ.13) GO TO 13
      IF(IJ.EQ.14) GO TO 14
      IF(IJ.EQ.23) GO TO 23
      IF(IJ.EQ.24) GO TO 24
      IF(IJ.EQ.34) GO TO 34
      IF(IJ.EQ.231) GO TO 231
      IF(IJ.EQ.312) GO TO 312
      IF(IJ.EQ.341) GO TO 341
      IF(IJ.EQ.413) GO TO 413
      IF(IJ.EQ.1234) GO TO 1234
      GOTO 100
   12 CONTINUE
      DO 10 L=1,N4     
         n123l=n123*(l-1)
      DO 10 K=1,N3      
         n12k=(k-1)*n12+n123l
      DO 10 I=1,N1      
         n1i=(i-1)*n1
      DO 10 J=1,I      
         ijkl=n12k+(j-1)*n1+i
         jikl=n12k+n1i+j
      X=A(ijkl)
      A(IJKL)=A(JIKL)
      A(jIKL)=X
   10 CONTINUE
      GO TO 100
   13 CONTINUE
      DO 20 L=1,N4      
         n123l=n123*(l-1)
      DO 20 I=1,N1      
         n12i=n123l+n12*(i-1)
      DO 20 J=1,N2      
         n1j=n1*(j-1)
      DO 20 K=1,I      
         ijkl=n123l+(k-1)*n12+n1j+i
         kjil=n12i+n1j+k
      X=A(IJKL)
      A(IJKL)=A(KJIL)
      A(KJIL)=X
   20 CONTINUE
      GO TO 100
   14 CONTINUE
      DO 25 I=1,N1
         n123i=n123*(i-1)
      DO 25 J=1,N2
         n1j=n1*(j-1)
      DO 25 K=1,N3
         n12k=n12*(k-1)+n1j
      DO 25 L=1,I
         ijkl=(l-1)*n123+n12k+i
         ljki=n123i     +n12k+l
      X=A(IJKL)
      A(IJKL)=A(LJKI)
      A(LJKI)= X
   25 CONTINUE
      GO TO 100
   23 CONTINUE
      DO 30 J=1,N2  
         j1=j-1
         n1j=n1*j1
         n12j=n12*j1
      DO 30 K=1,J   
         k1=k-1
         n1k=n1*k1
         n12k=n12*k1
      DO 30 L=1,N4
         n123l=n123*(l-1)
      DO 30 I=1,N1      
         n123li=n123l+i
         ijkl=n123li+n12k+n1j
         ikjl=n123li+n12j+n1k
      X=A(IJKL)
      A(IJKL)=A(IKJL)
      A(IKJL)=X
   30 CONTINUE
      GO TO 100
   24 CONTINUE
      DO 40 J=1,N2      
         j1=j-1
         n123j=n123*j1
         n1j  =n1*j1
      DO 40 L=1,J
         l1=l-1
         n123l=n123*l1
         n1l=n1*l1
      DO 40 K=1,N3
         n12k=n12*(k-1)
      DO 40 I=1,N1     
         n12ki=n12k+i
         ijkl=n123l+n12ki+n1j
         ilkj=n123j+n12ki+n1l
      X=A(IJKL)
      A(IJKL)=A(ILKJ)
      A(ILKJ)=X
   40 CONTINUE
      GO TO 100
   34 CONTINUE
      DO 50 K=1,N3      
         k1=k-1
         n12k=n12*k1
         n123k=n123*k1
      DO 50 L=1,K      
         l1=l-1
         n12l=n12*l1
         n123l=n123*l1
      DO 50 J=1,N2
         j1=j-1
         n1j=n1*(j-1)
      DO 50 I=1,N1    
         n1ji=n1j+i
         ijkl=n123l+n12k+n1ji
         ijlk=n123k+n12l+n1ji
      X=A(IJKL)
      A(IJKL)=A(IJLK)
      A(IJLK)=X
   50 CONTINUE
      GO TO 100
 231  CONTINUE
      DO 60 L=1,N4
         n123l=n123*(l-1)
      DO 60 J=1,N1
         j1=j-1
         n12j=n12*j1
         n1j=n1*j1
      DO 60 K=1,J
         n12jk=n12j+k
         k1=k-1
         n12k=n12*k1
         n1k=n1*k1
      DO 60 I=1,K
         i1=i-1
         ijkl=n123l+n12k+n1j+i
         jkil=n123l+i1*n12+n1k+j
         kijl=n123l+n12jk +i1*n1
      X=A(IJKL)
      A(IJKL)=A(JKIL)
      A(JKIL)=A(KIJL)
      A(KIJL)=X
      IF(J.EQ.K.OR.K.EQ.I) GOTO 60
         jikl=n123l+n12k  +n1*i1+j
         ikjl=n123l+n12j  +n1k  +i
         kjil=n123l+n12*i1+n1j  +k
      X=A(JIKL)
      A(JIKL)=A(IKJL)
      A(IKJL)=A(KJIL)
      A(KJIL)=X
 60   CONTINUE
      GOTO 100
 312  continue
      DO 70 L=1,N4
         n123l=n123*(l-1)
      DO 70 I=1,N1
         i1=i-1
         n1i=n1*i1
         n12i=n12*i1
      DO 70 J=1,I
         j1=j-1
         n1j =n1*j1
         n12j=n12*j1
      DO 70 K=1,J
         k1=k-1
         n1k=n1*k1
         n12k=n12*k1
         ijkl=n123l+n12k+n1j+i
         jkil=n123l+n12i+n1k+j
         kijl=n123l+n12j+n1i+k
      X=A(JKIL)
      A(JKIL)=A(IJKL)
      A(IJKL)=A(KIJL)
      A(KIJL)=X
      IF (I.EQ.J.OR.J.EQ.K)GOTO 70
         ikjl=n123l+n12j+n1k+i
         jikl=n123l+n12k+n1i+j
         kjil=n123l+n12i+n1j+k
      X=A(IKJL)
      A(IKJL)=A(JIKL)
      A(JIKL)=A(KJIL)
      A(KJIL)=X
 70   continue
      GOTO 100
 341  CONTINUE
      DO 80 L=1,N2
      DO 80 J=1,N1
      DO 80 K=1,J
      DO 80 I=1,K
         iljk=(k-1)*n123+(j-1)*n12+(l-1)*n1+i
         jlki=(i-1)*n123+(k-1)*n12+(l-1)*n1+j
         klij=(j-1)*n123+(i-1)*n12+(l-1)*n1+k
      X=A(ILJK)
      A(ILJK)=A(JLKI)
      A(JLKI)=A(KLIJ)
      A(KLIJ)=X
      IF(J.EQ.K.OR.K.EQ.I) GOTO 80
         ilkj=(j-1)*n123+(k-1)*n12+(l-1)*n1+i
         jlik=(k-1)*n123+(i-1)*n12+(l-1)*n1+j
         klji=(i-1)*n123+(j-1)*n12+(l-1)*n1+k
      X=A(JLIK)
      A(JLIK)=A(ILKJ)
      A(ILKJ)=A(KLJI)
      A(KLJI)=X
 80   CONTINUE
      GOTO 100
 413  CONTINUE
      DO 90 L=1,N2
      DO 90 I=1,N1
      DO 90 J=1,I
      DO 90 K=1,J
         jlki=(i-1)*n123+(k-1)*n12+(l-1)*n1+j
         iljk=(k-1)*n123+(j-1)*n12+(l-1)*n1+i
         klij=(j-1)*n123+(i-1)*n12+(l-1)*n1+k
      X=A(JLKI)
      A(JLKI)=A(ILJK)
      A(ILJK)=A(KLIJ)
      A(KLIJ)=X
      IF (I.EQ.J.OR.J.EQ.K)GOTO 90
         ilkj=(j-1)*n123+(k-1)*n12+(l-1)*n1+i
         jlik=(k-1)*n123+(i-1)*n12+(l-1)*n1+j
         klji=(i-1)*n123+(j-1)*n12+(l-1)*n1+k
      X=A(ILKJ)
      A(ILKJ)=A(JLIK)
      A(JLIK)=A(KLJI)
      A(KLJI)=X
 90   continue
      GO TO 100
 1234 CONTINUE
C      write(6,76)a
      DO 95 I=1,N1
      DO 95 J=1,N2
      DO 95 K=1,J
      DO 95 L=1,I
         ijkl=(l-1)*n123+(k-1)*n12+(j-1)*n1+i
         lkji=(i-1)*n123+(j-1)*n12+(k-1)*n1+l
      X=A(IJKL)
      A(IJKL)=A(LKJI)
      A(LKJI)=X
      if (i.eq.l.or.k.eq.j) goto 95
         ljki=(i-1)*n123+(k-1)*n12+(j-1)*n1+l
         ikjl=(l-1)*n123+(j-1)*n12+(k-1)*n1+i
      X=A(LJKI)
      A(LJKI)=A(IKJL)
      A(IKJL)=x
 95   CONTINUE
 100  CONTINUE
 77   FORMAT('IJKL:',4I3)
 76   format(4f15.10)
      call iexit(43)
      RETURN
      END
