      SUBROUTINE ADT2T4D(NO,NU,TI,T2,O2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER A,B
      COMMON /NEWTAPE/NT2T4,NT4INT
      COMMON/NEWCCSD/NTT2
      COMMON/ITRAT/ITER,MAXIT,ICCNV
      DIMENSION TI(NU,NU,NO),T2(NO,NU,NU,NO),O2(NO,NU,NU,NO)
      DATA ZERO/0.0D+0/,TRES/1.0D-10/,DAMP/0.35D+0/,ONE/1.0D+0/
      call ro2pph(0,no,nu,ti,o2)
      NLAST=NO+2
      ioff=1
      DO 1 I=1,NO
      call ro2hpp(0,no,nu,ti,o2)
      IASV=I+NLAST
      READ(NTT2,REC=IASV)TI
      DO 10 J=1,NO
      DO 10 A=1,NU
      DO 10 B=1,NU
      XN=TI(A,B,J)+T2(I,A,B,J)
      TI(A,B,J)=XN
 10   CONTINUE
      WRITE(NTT2,REC=IASV)TI
      call tranmd(ti,nu,nu,no,1,12)
      ioff=ioff+no
 1    CONTINUE
      RETURN
      END
      SUBROUTINE DENMT1(T1HP,O1HP,EH,EP,NH,NP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A
      COMMON/ITRAT/ITER,MAXIT,ICCNV
      COMMON/DMP/IDAMP,DAMP
      COMMON/NEWCCSD/NTT2
      DIMENSION T1HP(NH,NP),EH(NH),EP(NP),O1HP(NH,NP)
      DATA TRES/1.0D-10/,ZERO/0.0D+0/,ONE/1.0D+0/
      CALL RDT1HP(NH,NP,O1HP)
      DO 10 I=1,NH
      DO 10 A=1,NP
      XT1=T1HP(I,A)
      XDN=XT1/(EH(I)-EP(A))
      T1HP(I,A)=XDN
   10 CONTINUE
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
      SUBROUTINE REMD2(NH,NP,TI,T2HPPH,EH,EP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      COMMON/ITRAT/ITER,MAXIT,ICCNV
      DIMENSION T2HPPH(NH,NP,NP,NH),EH(NH),EP(NP),TI(NP,NP,NH)
      DATA ZERO/0.0D+0/,TRES/1.0D-10/,DAMP/0.0D+0/,ONE/1.0D+0/
      DO 10 I=1,NH
      DO 10 A=1,NP
      DO 10 B=1,NP
      DO 10 J=1,NH
      DEN=EH(I)+EH(J)-EP(A)-EP(B)
      XT2=T2HPPH(I,A,B,J)
      XDN=XT2*DEN
      T2HPPH(I,A,B,J)=XDN
   10 CONTINUE
      RETURN
      END
      subroutine drt3den(no,nu,t,eh,ep)
      implicit double precision (a-h,o-z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension t(*),eh(no),ep(nu)
      i1=1
      i2=i1+nu3
      call t3den(no,nu,t,t(i2),eh,ep)
      return
      end
      SUBROUTINE T3DEN(NO,NU,T3,O3,OEH,OEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C
      LOGICAL GIJK,IEJ,GABC,CGAB,JEK,BEC,AEB,AEC
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION T3(NU,NU,NU),OEH(NO),OEP(NU),O3(NU,NU,NU)
      COMMON/RESLTS/CMP(30)
      common/newopt/nopt(6)
      common/activ/noa,nua
      COMMON/DMP/IDAMP,DAMP
      COMMON/ITRAT/ITER,MAXIT,ICCNV
      DATA ZERO/0.0D+0/,TRES/1.0D-10/,ONE/1.0D+0/,half/0.5d+0/
      call ienter(57)
      iocc=no-noa
      kkias=0
      E3O=ZERO
      E3E=ZERO
      DO 10 I=1,NO
      DO 10 J=1,I
      DO 10 K=1,J
      if (i.ne.k)then
         kkias=kkias+1
      else 
         goto 10
      endif
      DIJK=OEH(I)+OEH(J)+OEH(K)
      IF(I.EQ.J.AND.J.EQ.K)THEN
      CALL ZEROMA(T3,1,NU3)
      CALL ZEROMA(O3,1,NU3)
      GO TO 9
      ENDIF
      if(i.le.iocc.and.nopt(6).eq.-1)then
         call zeroma(t3,1,nu3)
      goto 353
      else
         CALL RDVT3N(kkIAS,NU,T3)
      endif
      IF(ITER.GE.idamp)CALL RDVT3O(kkIAS,NU,O3)
 9    CONTINUE
      DO 111 C=1,NU
      DO 111 B=1,NU
      DO 111 A=1,NU
      IF (I.EQ.K)GOTO 111
      DENOM=DIJK-OEP(A)-OEP(B)-OEP(C)
      XDN=T3(A,B,C)/DENOM
      IF(ITER.GE.idamp)THEN
      T3(A,B,C)=DAMP*O3(A,B,C)+(ONE-DAMP)*XDN
       ELSE
      T3(A,B,C)=XDN
      ENDIF
      if(a.gt.nua.and.b.gt.nua.and.c.gt.nua.and.nopt(6).eq.-1)
     *t3(a,b,c)=zero
      IF (A.EQ.B.AND.B.EQ.C)T3(A,B,C)=ZERO
 111  CONTINUE
      IEJ = I.EQ.J
      JEK = J.EQ.K
      GIJK=I.GT.J.AND.J.GT.K
      DO 11 A=1,NU
      DO 11 B=1,A
      DO 11 C=1,B
      IF(I.EQ.K.OR.A.EQ.C)GOTO 11
      GABC=A.GT.B.AND.B.GT.C
      AEB=A.EQ.B
      BEC=B.EQ.C
      CGAB=C.GT.A.AND.C.GT.B
      IF (GIJK)THEN
      IF (GABC) GOTO 130
      IF (AEB)  GOTO 131
      IF (BEC)  GOTO 132
 130  CONTINUE
      X=T3(A,B,C)
      T3(A,B,C)=ZERO
      T3(A,C,B)=T3(A,C,B)-X
      T3(B,A,C)=T3(B,A,C)-X
      T3(B,C,A)=T3(B,C,A)-X
      T3(C,A,B)=T3(C,A,B)-X
      T3(C,B,A)=T3(C,B,A)-X
      GOTO 135
 131  CONTINUE
      X=T3(A,B,C)
      T3(A,B,C)=ZERO
      T3(A,C,B)=T3(A,C,B)-X
      T3(C,A,B)=T3(C,A,B)-X
      GOTO 135
 132  CONTINUE
      X=T3(A,B,C)
      T3(A,B,C)=ZERO
      T3(B,A,C)=T3(B,A,C)-X
      T3(B,C,A)=T3(B,C,A)-X
      GOTO 135
 135  CONTINUE
      ELSE
      IF (IEJ)THEN
      IF (GABC)GOTO 140
      IF (AEB) GOTO 141
      IF (BEC) GOTO 142
 140  CONTINUE
      X=T3(A,B,C)
      T3(A,B,C)=ZERO
      T3(A,C,B)=T3(A,C,B)-X
      T3(B,A,C)=T3(B,A,C)-X
      T3(B,C,A)=T3(B,C,A)-X
      T3(C,A,B)=T3(C,A,B)-X
      T3(C,B,A)=T3(C,B,A)-X
      GOTO 145
 141  CONTINUE
      X=T3(A,B,C)
      T3(A,B,C)=ZERO
      T3(A,C,B)=T3(A,C,B)-X
      T3(C,A,B)=T3(C,A,B)-X
      GOTO 145
 142  CONTINUE
      X=T3(A,B,C)
      T3(A,B,C)=ZERO
      T3(B,A,C)=T3(B,A,C)-X
      T3(B,C,A)=T3(B,C,A)-X
      GOTO 145
 145  CONTINUE
      ELSE
      IF(GABC)GOTO 150
      IF (AEB)GOTO 151
      IF (BEC)GOTO 152
 150  CONTINUE
      X=T3(A,B,C)
      T3(A,B,C)=ZERO
      T3(A,C,B)=T3(A,C,B)-X
      T3(B,A,C)=T3(B,A,C)-X
      T3(B,C,A)=T3(B,C,A)-X
      T3(C,A,B)=T3(C,A,B)-X
      T3(C,B,A)=T3(C,B,A)-X
      GOTO 155
 151  CONTINUE
      X=T3(A,B,C)
      T3(A,B,C)=ZERO
      T3(A,C,B)=T3(A,C,B)-X
      T3(C,A,B)=T3(C,A,B)-X
      GOTO 155
 152  CONTINUE
      X=T3(A,B,C)
      T3(A,B,C)=ZERO
      T3(B,A,C)=T3(B,A,C)-X
      T3(B,C,A)=T3(B,C,A)-X
      GOTO 155
 155  CONTINUE
      ENDIF
      ENDIF
 11   CONTINUE
353   continue
      iasv=kkias
      CALL WRVT3N(IASV,NU,T3)
 10   CONTINUE
 32   FORMAT('WIGNER RULE TRIPLES:ODD,EVN:',2F15.10)
      call iexit(57)
      RETURN
      END
      SUBROUTINE T4DEN(i,j,k,l,NU,DIJKL,T4,OEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      LOGICAL IEJKEL, AEBCED
      COMMON/ITRAT/ITER,MAXIT,ICCNV
      DIMENSION T4(NU,NU,NU,NU),OEP(NU)
      common/newopt/nopt(6)
      common/activ/noa,nua	
      DATA ZERO/0.0D+0/,TRES/1.0D-10/,DAMP/0.35D+0/,ONE/1.0D+0/,
     *TENTH/0.1D+0/,FIVE/5.0D+0/,THIRD/0.3333333D+0/,FIFTY/50.0d+0/
      IEJKEL=I.EQ.J.AND.K.EQ.L
      call ienter(64)
      DO 20 A=1,NU
      DO 20 B=1,NU
      DO 10 C=1,NU
      DO 10 D=1,NU
      iabcd=a+b+c+d
      DEN=DIJKL-OEP(A)-OEP(B)-OEP(C)-OEP(D)
      XT4=T4(D,C,B,A)
      XABCD=XT4/DEN
      if(nopt(6).lt.0)then
      if(a.gt.nua.and.b.gt.nua.and.c.gt.nua.or.
     *   a.gt.nua.and.b.gt.nua.and.d.gt.nua.or.
     *   a.gt.nua.and.c.gt.nua.and.d.gt.nua.or.
     *   b.gt.nua.and.c.gt.nua.and.d.gt.nua)xabcd=zero
      endif
      T4(D,C,B,A)=XABCD
 10   CONTINUE
      T4(A,A,A,B)=ZERO
      T4(A,A,B,A)=ZERO
      T4(A,B,A,A)=ZERO
      T4(B,A,A,A)=ZERO
      T4(A,B,B,B)=ZERO
      T4(B,A,B,B)=ZERO
      T4(B,B,A,B)=ZERO
      T4(B,B,B,A)=ZERO
 20   CONTINUE
 1000 format('i,j,k,l,a,b,c,d:',8i3,f15.10)
      DO 30 A=1,NU
      DO 30 B=1,A
      DO 30 C=1,B
      DO 30 D=1,C
      AEBCED=A.EQ.B.AND.C.EQ.D
      IF (IEJKEL.AND.AEBCED)THEN
      X=T4(A,A,C,C)
      T4(A,A,C,C)=ZERO
      T4(A,C,A,C)=T4(A,C,A,C)-X
      T4(A,C,C,A)=T4(A,C,C,A)-X
      T4(C,A,A,C)=T4(C,A,A,C)-X
      T4(C,A,C,A)=T4(C,A,C,A)-X
      T4(C,C,A,A)=T4(C,C,A,A)-X
      ENDIF
 30   CONTINUE
      call iexit(64)
      RETURN
      END
