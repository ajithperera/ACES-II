      subroutine drf1int(no,nu,t)
      implicit double precision(a-h,o-z)
      logical print
      common /flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/itrat/icycle,mx,icn
      common/totmem/nmem
      dimension t(*)
      print=iflags(1).gt.10
      i1=1            !o1
      i2=i1+nou       !fpp
      i3=i2+nu2       !ti
      ni3=nu3
      if(nu3.lt.(no2u2/2+1))ni3=no2u2/2+1
      i4=i3+ni3       !voe
      i5=i4+no2u2     !vdppph
      itot=i5+nou3 
      if(icycle.eq.1.and.print)write(6,99)nmem,itot
 99   format('Space usage in f1int: available - ',i8,'   used - ',i8)
      call f1int(no,nu,t,t(i2),t(i3),t(i4),t(i5))
      return
      end
      SUBROUTINE F1INT(NO,NU,O1,FPP,TI,VOE,VDPPPH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION O1(1),FPP(1),TI(1),VOE(1),VDPPPH(1)
      call ienter(7)
      NLAST=5*NO+2*NU
      call ro1(nou,1,o1)
c
      call rf(nou,nlast+2,fpp)
      CALL RO2HPH(1,NO,NU,TI,VOE)
      CALL SYMT21(VOE,NO,NU,NO,NU,13)
      CALL MATMUL(VOE,   O1,FPP,NOU,1,NOU, 0,0)
      call wr(nall4,nlast+4,nou,fpp)
      call rf(no2,nlast+1,fpp)
      CALL RDVEM4(1,NU,NO,TI,VDPPPH)
      CALL SYMT21(VDPPPH,NO,NO,NO,NU,23)
      CALL MATMUL(VDPPPH,O1,FPP,NO2,1,NOU, 0,0)
      CALL TRT1(NO,NU,TI,O1)
      CALL TRANSQ(FPP,NO)
      CALL RO2PHP(0,NO,NU,TI,VDPPPH)
      CALL TRANMD(VDPPPH,NU,NO,NU,NO,13)
      CALL MATMUL(VOE,  VDPPPH,  FPP,NO,NO,NOU2,0,0)
      call rf(nou,nlast+4,ti)
      CALL MATMUL(ti,   O1,FPP,NO,NO,NU,  0,0)
      call wr(nall4,nlast+5,no2,fpp)
      call rf(nu2,nlast+3,fpp)
      CALL MATMUL(VDPPPH,    VOE, FPP,NU, NU,NO2U,0,1)
      CALL RDVEM4(0,NO,NU,TI,VDPPPH)
      CALL SYMT21(VDPPPH,NU,NU,NU,NO,23)
      call ro1(nou,1,o1)
      CALL TRT1(NO,NU,TI,O1)
      CALL MATMUL(VDPPPH,O1,FPP,NU2,1, NOU,0,0)
      call rf(nou,nlast+4,ti)
      CALL MATMUL(O1,  ti, FPP,NU, NU,NO, 0,1)
      CALL TRANSQ(FPP,NU)
      call wr(nall4,nlast+6,nu2,fpp)
      call iexit(7)
      RETURN
      END
      subroutine drintqua(no,nu,t)
      implicit double precision(a-h,o-z)
      logical print
      common /flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/itrat/icycle,mx,icn
      common/totmem/mem
      dimension t(*)
      print=iflags(1).gt.10
      i1=1            !ti
      i2=i1+nu3       !t2pp
      i3=i2+no2u2     !t2
      i4=i3+no2u2     !voe
      i5=i4+no2u2     !v
      i6=i5+nu4       !ve
      i7=i6+nou3
      itot=i7+nou
      if(icycle.eq.1.and.print) write(6,99)mem,itot
 99   format('Space usage in intqua: available - ',i8,'   used - ',i8)
      call intqua(0,no,nu,t,t(i2),t(i3),t(i4),t(i5),t(i6),t(i7))
      return
      end
      SUBROUTINE INTQUA(IC,NO,NU,TI,T2PP,T2,VOE,V,VE,O1HP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SDTQ2,SDTQ3,SDT
      COMMON/NEWOPT/NOPT(6)
      common/wpak/nfr(15),nsz(15)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION TI(1),T2PP(1),T2(1),VOE(1),V(1),VE(1),O1HP(1)
      DATA TWO/2.0D+00/,HALF/0.5D+00/,ONE/1.0D+0/
      call ienter(11)
      SDT  =NOPT(1).GE.2
      SDTQ2=NOPT(1).EQ.4
      SDTQ3=NOPT(1).EQ.5
      CALL RO2PPH(0,NO,NU,TI,T2PP)
      CALL TRANMD(T2PP,NU,NU,NO,NO,34)
      CALL RDT1HP(NO,NU,O1HP)
      CALL RO2HHP(1,NO,NU,TI,VOE)
      CALL TRANMD(VOE,NO,NO,NU,NU,34)
      call ienter(121)
      if(.not.sdt)goto 3000
      CALL RDOV4(0,NO,NU,TI,V)
      CALL RDVEM4(0,NO,NU,TI,VE)
      IF (SDTQ3.AND.IC.EQ.1) GO TO 10
      CALL RDVEM4(-10,NO,NU,TI,VE)
      CALL VECMUL(VE,NOU3,HALF)
 10   continue
C
      IF(SDT)THEN
      call mtrans(v,nu,2)
      CALL MATMUL(T2PP,VOE,V ,NU2,NU2,NO2,0,0)
      IF(SDTQ2)THEN
      CALL WROV4A(0,NO,NU,TI,V)
      ENDIF
      CALL MTRANS(V,NU,3)
      ENDIF
      CALL VECMUL(V,NU4,HALF)
      CALL MATMUL(VE,O1HP,V,NU3,NU,NO,0,1)
      call symv(v,nu)
      if(sdt)then
      CALL MTRANS(V,NU,1)
      CALL MTRANS(V,NU,7)
      endif
      IF (SDTQ3.AND.IC.EQ.1)THEN
      CALL WROV4A(0,NO,NU,TI,V)
      ELSE
      CALL WROV4(0,NO,NU,TI,V)
      ENDIF
 3000 continue
      call iexit(121)
      call ienter(122)
      CALL RDOV4(1,NU,NO,TI,V)
      CALL MTRANS(V,NO,7)
      CALL RDVEM4(1,NU,NO,TI,VE)
      IF (SDTQ3.AND.IC.EQ.1) GO TO 11
      CALL RDVEM4(-11,NU,NO,TI,VE)
      CALL VECMUL(VE,NO3U,HALF)
 11   CONTINUE
      call tranmd(v,no,no,no,no,312)
      CALL MATMUL(VOE,T2PP,V,NO2,NO2,NU2,0,0)
      IF(SDTQ2)CALL WROV4A(1,NU,NO,TI,V)
      CALL VECMUL(V,NO4,HALF)
      call tranmd(v,no,no,no,no,231)
      CALL TRT1(NO,NU,TI,O1HP)
      CALL MATMUL(VE,O1HP,V,NO3,NO,NU,0,0)
      call symv(v,no)
      call tranmd(v,no,no,no,no,312)
      IF (SDTQ3.AND.IC.EQ.1) THEN
      CALL WROV4A(1,NU,NO,TI,V)
      ELSE
       CALL WROV4(1,NU,NO,TI,V)
      ENDIF
      call iexit(122)
      call ienter(123)
      CALL RO2HPP(1,NO,NU,TI,V)
      IR10=10
      IF(SDTQ3.AND.IC.EQ.1)IR10=0
      CALL RDVEM4(IR10,NO,NU,TI,VE)
      CALL TRANMD(VE,NU,NU,NU,NO,12)
      CALL TRT1(NU,NO,TI,O1HP)
      CALL MATMUL(O1HP,VE,V,NO,NOU2,NU,0,0)
      CALL TRANMD(V,NO,NU,NU,NO,1234)
      CALL TRANSQ(V,NOU)
      CALL RDVEM4(1,NU,NO,TI,VE)
      CALL TRT1(NO,NU,TI,O1HP)
      CALL MATMUL(O1HP,VE,V,NU,NO2U,NO,0,1)
      CALL TRANSQ(V,NOU)
      CALL TRANMD(V,NO,NU,NU,NO,1234)
      CALL WRVOEP(1,NO,NU,TI,V)
c
      CALL RO2HPH(0,NO,NU,TI,T2PP)
      call reapak(3,voe,ve)
      CALL VECMUL(T2PP,NO2U2,HALF)
      CALL SYMT21(t2pp,no,nu,no,nu,13)
      call zeroma(t2,1,no2u2)
      call mmf(nou,nou,nou,nsz(3),t2pp,voe,ve,t2,1,0)
      CALL desm21(t2pp,no,nu,no,nu,13)
      CALL TRANMDpak(ve,no,nu,nu,no,nsz(3),23)
      call mmf(nou,nou,nou,nsz(3),t2pp,voe,ve,t2,0,1)
      call vecadd(v,T2,no2u2)
      call wrvoet2(no,nu,ti,v)
      IF (.NOT.SDT)GOTO 300
      call ro2hpp(1,no,nu,ti,voe)
      call vecadd(voe,T2,no2u2)
      call vecadd(voe,T2,no2u2)
      IF(SDTQ2)CALL WRVOEA(NO,NU,TI,VOE)
      call vecadd(v,T2,no2u2)
      IF(SDTQ3.AND.IC.EQ.1) THEN
      CALL WRVOEA(NO,NU,TI,V)
      ELSE
       CALL WRVOE(NO,NU,TI,V)
      ENDIF
c
 300  CONTINUE
      call iexit(123)
      call ienter(124)
      CALL RO2HPP(2,NO,NU,TI,T2)
      CALL TRANSQ(T2,NOU)
      CALL RDVEM4(1,NU,NO,TI,V)
      CALL TRANMD(V,NO,NO,NO,NU,13)
      CALL MATMUL(O1HP,V,T2,NU,NO2U,NO,0,1)
      CALL TRANSQ(T2,NOU)
      CALL TRANMD(T2,NO,NU,NU,NO,1234)
      CALL WRVOEP(2,NO,NU,TI,T2)
c
      CALL TRANMD(T2PP,NO,NU,No,Nu,13)
      CALL VECMUL(T2PP,NO2U2,TWO)
      call ro2hpp(1,no,nu,ti,voe)
      call tranmd(voe,no,nu,nu,no,23)
      CALL MATMUL(T2PP,VOE,T2,NOU,NOU,NOU,1,1)
      IF(SDTQ2) THEN
      CALL RO2HPP(2,NO,NU,TI,VOE)
      CALL VECADD(VOE,T2,NO2U2)
      CALL WRVOA(NO,NU,TI,VOE)
      ENDIF
      call rovoep(2,no,nu,ti,voe)
      CALL VECMUL(T2,NO2U2,half)
      call vecadd(voe,T2,no2u2)
      CALL TRANMD(voe,NO,NU,NU,NO,23)
      CALL RDVEM4(IR10,NO,NU,TI,VE)
      CALL TRANMD(VE,NU,NU,NU,NO,13)
      CALL TRT1(NU,NO,TI,O1HP)
      CALL MATMUL(O1HP,VE,voe,NO,NOU2,NU,0,0)
      CALL TRANMD(voe,NO,NU,NU,NO,23)
      call wrvoot2(no,nu,ti,voe)
      IF (.NOT.SDT)goto 9999
      call vecadd(voe,T2,no2u2)
      IF(SDTQ3.AND.IC.EQ.1) THEN
      CALL WRVOA(NO,NU,TI,voe)
      ELSE
       CALL WRVO(NO,NU,TI,voe)
      ENDIF
      CALL RO2HPP(2,NO,NU,TI,T2)
      CALL RDVEM4(0,NO,NU,TI,VE)
      CALL TRANMD(VE,NU,NU,NU,NO,13)
      CALL MATMUL(O1HP,VE,T2,NO,NOU2,NU,0,0)
      CALL WRVOEP(3,NO,NU,TI,T2)
9999  continue
      call iexit(124)
      call iexit(11)
      RETURN
      END
      subroutine drt4fht4(i,j,k,l,no,nu,o1)
      implicit double precision (a-h,o-z)
      logical print
      common /flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension o1(1)
      common/itrat/icycle,mx,icn
      common/totmem/mem
      print=iflags(1).gt.10
      i1=1          !t2
      i2=i1+no2u2   !t4
      i3=i2+nu4     !o4
      i4=i3+no*nu4  !t4n
      i5=i4+nu4     !fhh
      i6=i5+no2 
      if(icycle.eq.2.and.i.eq.2.and.print)write(6,98)mem,i6
 98   format('Space usage in t4fht4: available - ',i10,'   used - ',i10)
      call t4fht4(i,j,k,l,no,nu,o1(i3),o1(i2),o1(i4),o1(i5))
      return
      end
      SUBROUTINE T4FHT4(I,J,K,L,NO,NU,O4,T4,T4N,FHH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D,E,F
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION O4(NU,NU,NU,NU,NO),T4(NU,NU,NU,NU),T4N(NU,NU,NU,NU),
     *FHH(NO,NO)
      DATA ZERO/0.0D+0/,TWO/2.0D+0/,half/0.5d+0/
      call ienter(65)
      NLAST=5*NO+2*NU
      IASV=NLAST+5
      READ(NALL4,REC=IASV)FHH
      DO 10 M=1,NO
      IF (I.EQ.J.AND.J.EQ.M.OR.I.EQ.K.AND.K.EQ.M
     *.OR.J.EQ.K.AND.K.EQ.M)THEN
      CALL ZEROMA(O4(1,1,1,1,M),1,NU4)
      GOTO 10
      ENDIF
      CALL RDIJKM(I,J,K,M,NO,NU,T4N,O4(1,1,1,1,M))
 10   CONTINUE
      if(k.eq.l)call vecmul(t4,nu4,half)
      CALL MATMUL(O4,FHH(1,L),T4,NU4,1,NO,0,1)
      if(k.eq.l)then
      call symt411(t4,nu,1)
      goto 111
      endif
      DO 11 M=1,NO
      IF (I.EQ.J.AND.J.EQ.M.OR.I.EQ.M.AND.M.EQ.L
     *.OR.J.EQ.M.AND.M.EQ.L)THEN
      CALL ZEROMA(O4(1,1,1,1,M),1,NU4)
      GOTO 11
      ENDIF
      CALL RDIJMK(I,J,L,M,NO,NU,T4N,O4(1,1,1,1,M))
 11   CONTINUE
      if(j.eq.k)call vecmul(t4,nu4,half)
      CALL MATMUL(O4,FHH(1,K),T4,NU4,1,NO,0,1)
 111  continue
      if(j.eq.k)then
         call symt411(t4,nu,7)
         goto 112
      endif
      DO 12 M=1,NO
      IF (I.EQ.M.AND.M.EQ.K.OR.I.EQ.M.AND.M.EQ.L
     *.OR.M.EQ.K.AND.K.EQ.L)THEN
      CALL ZEROMA(O4(1,1,1,1,M),1,NU4)
      GOTO 12
      ENDIF
      CALL RDIMJK(I,K,L,M,NO,NU,T4N,O4(1,1,1,1,M))
 12   CONTINUE
      if(i.eq.j)call vecmul(t4,nu4,half)
      CALL MATMUL(O4,FHH(1,J),T4,NU4,1,NO,0,1)
 112  continue
      if(i.eq.j)then
         call symt411(t4,nu,6)
         goto 113
      endif
      DO 13 M=1,NO
      IF (M.EQ.J.AND.J.EQ.K.OR.M.EQ.J.AND.J.EQ.L
     *.OR.M.EQ.K.AND.K.EQ.L)THEN
      CALL ZEROMA(O4(1,1,1,1,M),1,NU4)
      GOTO 13
      ENDIF
      CALL RDMIJK(J,K,L,M,NO,NU,T4N,O4(1,1,1,1,M))
 13   CONTINUE
      CALL MATMUL(O4,FHH(1,I),T4,NU4,1,NO,0,1)
 113  continue
      call iexit(65)
      RETURN
      END
      subroutine drt4fpt4(i,j,k,l,no,nu,o1)
      implicit double precision (a-h,o-z)
      logical print
      common /flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension o1(1)
      common/itrat/icycle,mx,icn
      common/totmem/mem
      print=iflags(1).gt.10
      i1=1          !t2
      i2=i1+no2u2   !t4
      i3=i2+nu4     !o4
      i4=i3+nu4     !fpp
      it=i4+nu2
      if(icycle.eq.2.and.i.eq.2.and.print)write(6,98)mem,it
 98   format('Space usage in t4fpt4: available - ',i10,'   used - ',i10)
      call t4fpt4(i,j,k,l,no,nu,o1(i3),o1(i2),o1(i4))
      return
      end
      SUBROUTINE T4FPT4(I,J,K,L,NO,NU,O4,T4,FPP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D,E,F
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION O4(NU,NU,NU,NU),T4(NU,NU,NU,NU),FPP(NU,NU)
      DATA ZERO/0.0D+0/,HALF/0.5D+0/
      call ienter(66)
      NLAST=5*NO+2*NU
      IASV=NLAST+6
      READ(NALL4,REC=IASV)FPP
      IAS=IT4(I,J,K,L)
      CALL RDT4(IAS,NU,O4)
      CALL MATMUL(O4,FPP,T4,NU3,NU,NU,0,0)
      CALL MTRANS(O4,NU,1)
      CALL MTRANS(T4,NU,1)
      CALL MATMUL(O4,FPP,T4,NU3,NU,NU,0,0)
      CALL MTRANS(O4,NU,4)
      CALL MTRANS(T4,NU,4)
      CALL MATMUL(O4,FPP,T4,NU3,NU,NU,0,0)
      CALL MTRANS(O4,NU,23)
      CALL MTRANS(T4,NU,23)
      CALL MATMUL(O4,FPP,T4,NU3,NU,NU,0,0)
      CALL MTRANS(T4,NU,11)
      call iexit(66)
      RETURN
      END
      subroutine drt4hht4(i,j,k,l,no,nu,o1)
      implicit double precision (a-h,o-z)
      logical print
      common /flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension o1(1)
      common/itrat/icycle,mx,icn
      common/totmem/mem
      print=iflags(1).gt.10
      i1=1          !t2
      i2=i1+no2u2   !t4
      i3=i2+nu4     !o4
      i4=i3+nu4*no  !t4n
      i5=i4+nu4     !pz
      it=i5+no4
      if(icycle.eq.2.and.i.eq.2.and.print)write(6,98)mem,it
 98   format('Space usage in t4hht4: available - ',i10,'   used - ',i10)
      call t4hht4(i,j,k,l,no,nu,o1(i3),o1(i2),o1(i4),o1(i5))
      return
      end
      SUBROUTINE T4HHT4(I,J,K,L,NO,NU,O4,T4,T4N,PZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL NE1,NE2,IEJ,JEK,KEL
      COMMON/NEWOPT/NOPT(6)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION O4(NU,NU,NU,NU,NO),T4(NU,NU,NU,NU),PZ(NO,NO,NO,NO),
     *     T4N(NU,NU,NU,NU)
      DATA ZERO/0.0D+0/,HALF/0.5D+0/
      call ienter(67)
      IEJ=I.EQ.J
      JEK=J.EQ.K
      KEL=K.EQ.L
      I11=1
      IF(NOPT(1).GT.5)I11=11
      CALL RDOV4(I11,NU,NO,T4N,PZ)
      if(i.eq.j)goto 500
      if(j.eq.k)call vecmul(t4,nu4,half)
      DO 101 N=1,NO
      NE1=N.EQ.K
      NE2=N.EQ.L
      DO 10 M=1,NO
      IF (M.EQ.N.AND.(NE1.OR.NE2).OR.KEL.AND.(M.EQ.K.OR.NE1))THEN
      CALL ZEROMA(O4(1,1,1,1,M),1,NU4)
      GO TO 10
      ENDIF
      CALL RDMNKL(K,L,M,N,NO,NU,T4N,O4(1,1,1,1,M))
 10   CONTINUE
      CALL MATMUL(O4,PZ(1,n,I,J),T4,NU4,1,NO,0,0)
 101  continue
      if(j.eq.k)then
      call symt411(t4,nu,7)
      goto1111
      endif
      if(k.eq.l)call vecmul(t4,nu4,half)
      DO 111 N=1,NO
      NE1=N.EQ.J
      NE2=N.EQ.L
      DO 11 M=1,NO
      IF (M.EQ.N.AND.(NE1.OR.NE2))THEN
      CALL ZEROMA(O4(1,1,1,1,M),1,NU4)
      GO TO 11
      ENDIF
      CALL RDMJNL(J,L,M,N,NO,NU,T4N,O4(1,1,1,1,M))
 11   CONTINUE
      CALL MATMUL(O4,PZ(1,N,I,K),T4,NU4,1,NO,0,0)
 111  CONTINUE
 1111 continue
      if(k.eq.l)then
      call symt411(t4,nu,1)
      goto 1211
      endif
      DO 121 N=1,NO
      NE1=N.EQ.J
      NE2=N.EQ.K
      DO 12 M=1,NO
      IF (M.EQ.N.AND.(NE1.OR.NE2).OR.JEK.AND.(NE2.OR.M.EQ.K))THEN
      CALL ZEROMA(O4(1,1,1,1,M),1,NU4)
      GO TO 12
      ENDIF
      CALL RDMJKN(J,K,M,N,NO,NU,T4N,O4(1,1,1,1,M))
 12   CONTINUE
      CALL MATMUL(O4,PZ(1,N,I,L),T4,NU4,1,NO,0,0)
 121  CONTINUE
 1211 continue
      if(k.eq.l)call vecmul(t4,nu4,half)
      DO 131 N=1,NO
      NE1=N.EQ.I
      NE2=N.EQ.L
      DO 13 M=1,NO
      IF(M.EQ.N.AND.(NE1.OR.NE2))THEN
      CALL ZEROMA(O4(1,1,1,1,M),1,NU4)
      GO TO 13
      ENDIF
      CALL RDIMNL(I,L,M,N,NO,NU,T4N,O4(1,1,1,1,M))
 13   CONTINUE
      CALL MATMUL(O4,PZ(1,N,J,K),T4,NU4,1,NO,0,0)
 131  CONTINUE
      if(k.eq.l)then 
      call symt411(t4,nu,1)
      goto 1411
      endif
      if(j.eq.k)call vecmul(t4,nu4,half)
      DO 141 N=1,NO
      NE1=N.EQ.I
      NE2=N.EQ.K
      DO 14 M=1,NO
      IF (M.EQ.N.AND.(NE1.OR.NE2))THEN
      CALL ZEROMA(O4(1,1,1,1,M),1,NU4)
      GO TO 14
      ENDIF
      CALL RDIMKN(I,K,M,N,NO,NU,T4N,O4(1,1,1,1,M))
 14   CONTINUE
      CALL MATMUL(O4,PZ(1,N,J,L),T4,NU4,1,NO,0,0)
 141  CONTINUE
 1411 CONTINUE
      if(j.eq.k)then 
      call symt411(t4,nu,7)
      goto 1511
      endif
      DO 151 N=1,NO
      NE1=N.EQ.I
      NE2=N.EQ.J
      DO 15 M=1,NO
      IF (M.EQ.N.AND.(NE1.OR.NE2).OR.IEJ.AND.(NE1.OR.M.EQ.I))THEN
      CALL ZEROMA(O4(1,1,1,1,M),1,NU4)
      GO TO 15
      ENDIF
      CALL RDIJMN(I,J,M,N,NO,NU,T4N,O4(1,1,1,1,M))
 15   CONTINUE
      CALL MATMUL(O4,PZ(1,N,K,L),T4,NU4,1,NO,0,0)
 151  CONTINUE
 1511 continue
      goto 1000
 500  continue
      DO 1019 N=1,NO
      NE1=N.EQ.K
      NE2=N.EQ.L
      DO 109 M=1,NO
      IF (M.EQ.N.AND.(NE1.OR.NE2).OR.KEL.AND.(M.EQ.K.OR.NE1))THEN
      CALL ZEROMA(O4(1,1,1,1,M),1,NU4)
      GO TO 109
      ENDIF
      CALL RDMNKL(K,L,M,N,NO,NU,T4N,O4(1,1,1,1,M))
 109  CONTINUE
      CALL MATMUL(O4,PZ(1,n,I,J),T4,NU4,1,NO,0,0)
 1019 continue
      call vecmul(t4,nu4,half)
      DO 1119 N=1,NO
      NE1=N.EQ.J
      NE2=N.EQ.L
      DO 119 M=1,NO
      IF (M.EQ.N.AND.(NE1.OR.NE2))THEN
      CALL ZEROMA(O4(1,1,1,1,M),1,NU4)
      GO TO 119
      ENDIF
      CALL RDMJNL(J,L,M,N,NO,NU,T4N,O4(1,1,1,1,M))
 119  CONTINUE
      CALL MATMUL(O4,PZ(1,N,I,K),T4,NU4,1,NO,0,0)
 1119 CONTINUE
      call symt411(t4,nu,6)
      call vecmul(t4,nu4,half)
      DO 1219 N=1,NO
      NE1=N.EQ.J
      NE2=N.EQ.K
      DO 129 M=1,NO
      IF (M.EQ.N.AND.(NE1.OR.NE2).OR.JEK.AND.(NE2.OR.M.EQ.K))THEN
      CALL ZEROMA(O4(1,1,1,1,M),1,NU4)
      GO TO 129
      ENDIF
      CALL RDMJKN(J,K,M,N,NO,NU,T4N,O4(1,1,1,1,M))
 129  CONTINUE
      CALL MATMUL(O4,PZ(1,N,I,L),T4,NU4,1,NO,0,0)
 1219 CONTINUE
      call symt411(t4,nu,6)
      DO 1519 N=1,NO
      NE1=N.EQ.I
      NE2=N.EQ.J
      DO 159 M=1,NO
      IF (M.EQ.N.AND.(NE1.OR.NE2).OR.IEJ.AND.(NE1.OR.M.EQ.I))THEN
      CALL ZEROMA(O4(1,1,1,1,M),1,NU4)
      GO TO 159
      ENDIF
      CALL RDIJMN(I,J,M,N,NO,NU,T4N,O4(1,1,1,1,M))
 159  CONTINUE
      CALL MATMUL(O4,PZ(1,N,K,L),T4,NU4,1,NO,0,0)
 1519 CONTINUE
 1000 continue
      call iexit(67)
      RETURN
      END
      subroutine drt4hpt4(i,j,k,l,no,nu,o1)
      implicit double precision (a-h,o-z)
      logical print
      common /flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension o1(1)
      common/itrat/icycle,mx,icn
      common/totmem/mem
      print=iflags(1).gt.10
      i1=1          !t2
      i2=i1+no2u2   !t4
      i3=i2+nu4     !o4
      i4=i3+nu4*no  !voe
      i5=i4+no2u2   !vo
      i6=i5+no2u2   !t4n
      it=i6+nu4
      if(icycle.eq.2.and.i.eq.2.and.print)write(6,98)mem,it
 98   format('Space usage in t4hp4: available - ',i10,'   used - ',i10)
      call t4hpt4(i,j,k,l,no,nu,o1(i3),o1(i2),o1(i4),o1(i5),o1(i6))
      return
      end
      SUBROUTINE T4HPT4(I,J,K,L,NO,NU,O4,T4,VOE,VO,T4N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D,E,F
      LOGICAL IJ,IK,IL,JK,JL,KL
      COMMON/NEWOPT/NOPT(6)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION O4(NU,NU,NU,NU,NO),T4(NU,NU,NU,NU),VOE(NU,NO,NU,NO),
     *VO(NU,NO,NU,NO),T4N(NU,NU,NU,NU)
      DATA ZERO/0.0D+0/,TWO/2.0D+0/,HALF/0.5D+00/
      call ienter(68)
      IJ=I.EQ.J
      IK=I.EQ.K
      IL=I.EQ.L
      JK=J.EQ.K
      JL=J.EQ.L
      KL=K.EQ.L
      I11=1
      IF (NOPT(1).GT.5)I11=11
      I12=I11+1
      CALL RO2PHP(I11,NO,NU,T4N,VOE)
      CALL RO2PHP(I12,NO,NU,T4N,VO)
      CALL TRANSQ(VOE,NOU)
      CALL TRANSQ(VO,NOU)
      if(k.eq.l)call vecmul(t4,nu4,half)
      DO 713 M=1,NO
      IF (IJ.AND.(JK.OR.J.EQ.M).OR.(IK.OR.JK).AND.K.EQ.M)then
      call zeroma(o4(1,1,1,1,m),1,nu4)
      GOTO 713
      endif
      CALL RDIJKM(I,J,K,M,NO,NU,T4N,O4(1,1,1,1,M))
 713  CONTINUE
      CALL MATMUL(O4,VO(1,1,1,L),T4,NU3,NU,NOU,0,1)
      CALL MTRANS(T4,NU,1)
      CALL MTRSM(O4,NO,NU,1)
      if(j.eq.k)call vecmul(t4,nu4,half)
      CALL MATMUL(O4,VO(1,1,1,L),T4,NU3,NU,NOU,0,1)
c
      if(j.eq.k)then
      CALL MTRANS(T4,NU,1)
      CALL MTRSM(O4,NO,NU,5)
      call symt411(t4,nu,7)
      else
      CALL MTRANS(T4,NU,3)
      CALL MTRSM(O4,NO,NU,3)
      if(i.eq.j)call vecmul(t4,nu4,half)
      CALL MATMUL(O4,VO(1,1,1,L),T4,NU3,NU,NOU,0,1)
      CALL MTRANS(T4,NU,10)
      CALL MTRSM(O4,NO,NU,16)
      endif
      if(i.eq.j)then
      call symt411(t4,nu,6)
      CALL MTRANS(T4,NU,1)
      else
      CALL MTRANS(T4,NU,11)
      CALL MATMUL(O4,VO(1,1,1,L),T4,NU3,NU,NOU,0,1)
      CALL MTRANS(T4,NU,14)
      endif
      CALL MTRSM(O4,NO,NU,11)
      call smijkm(i,j,k,no,nu,t4n,o4)
      CALL MTRANS(T4,NU,1)
      CALL MATMUL(O4,VOE(1,1,1,L),T4,NU3,NU,NOU,0,0)
      if(k.eq.l)then
      call symt411(t4,nu,1)
      call mtrans(t4,nu,15)
      goto 1312
      endif
      if(j.eq.k)call vecmul(t4,nu4,half)
      DO 712 M=1,NO
      IF (IJ.AND.(JL.OR.J.EQ.M).OR.L.EQ.M.AND.(IL.OR.JL))THEN
      call zeroma(o4(1,1,1,1,m),1,nu4)
      GOTO 712
      ENDIF
      CALL RDIJMK(I,J,L,M,NO,NU,T4N,O4(1,1,1,1,M))
 712   CONTINUE
      CALL MATMUL(O4,VO(1,1,1,K),T4,NU3,NU,NOU,0,1)
      CALL MTRSM(O4,NO,NU,1)
      CALL MTRANS(T4,NU,1)
      CALL MATMUL(O4,VO(1,1,1,K),T4,NU3,NU,NOU,0,1)
      CALL MTRSM(O4,NO,NU,5)
      CALL MTRANS(T4,NU,5)
      if(i.eq.j)call vecmul(t4,nu4,half)
      CALL MATMUL(O4,VO(1,1,1,K),T4,NU3,NU,NOU,0,1)
      if(i.eq.j)then
      CALL MTRANS(T4,NU,11)
      call symt411(t4,nu,6)
      CALL MTRANS(T4,NU,1)
      CALL MTRSM(O4,NO,NU,14)
      else
      CALL MTRSM(O4,NO,NU,15)
      CALL MTRANS(T4,NU,15)
      CALL MATMUL(O4,VO(1,1,1,K),T4,NU3,NU,NOU,0,1)
      CALL MTRSM(O4,NO,NU,2)
      CALL MTRANS(T4,NU,2)
      endif
      call smijmk(i,j,l,no,nu,t4n,o4)
      CALL MATMUL(O4,VOE(1,1,1,K),T4,NU3,NU,NOU,0,0)
      CALL MTRANS(T4,NU,22)
 1312 continue
      if(j.eq.k)then
      call symt411(t4,nu,1)
      call mtrans(t4,nu,16)
      goto 1711
      endif
      if(i.eq.j)call vecmul(t4,nu4,half)
      DO 711 M=1,NO
      IF(IK.AND.(KL.OR.K.EQ.M).OR.L.EQ.M.AND.(IL.OR.KL))THEN
      CALL ZEROMA(O4(1,1,1,1,M),1,NU4)
      GOTO 711
      ENDIF
      CALL RDIMJKN(I,K,L,M,NO,NU,T4N,O4(1,1,1,1,M))
 711   CONTINUE
      CALL MTRANS(t4,nu,10)
      CALL MTRSM(O4,NO,NU,11)
      CALL MATMUL(O4,VO(1,1,1,J),T4,NU3,NU,NOU,0,1)
      CALL MTRANS(t4,nu,15)
      CALL MTRSM(O4,NO,NU,15)
      CALL MATMUL(O4,VO(1,1,1,J),T4,NU3,NU,NOU,0,1)
      CALL MTRSM(O4,NO,NU,2)
      CALL MTRANS(T4,NU,2)
      if(k.eq.l)call vecmul(t4,nu4,half)
      CALL MATMUL(O4,VO(1,1,1,J),T4,NU3,NU,NOU,0,1)
      CALL MTRSM(O4,NO,NU,1)
      CALL MTRANS(T4,NU,1)
      if(k.eq.l)then
         call symt411(t4,nu,1)
      else
      CALL MATMUL(O4,VO(1,1,1,J),T4,NU3,NU,NOU,0,1)
      endif
      call mtrsm(o4,no,nu,3)
      call smimjk(i,k,l,no,nu,t4n,o4)
      CALL MTRANS(T4,NU,3)
      CALL MATMUL(O4,VOE(1,1,1,J),T4,NU3,NU,NOU,0,0)
      CALL MTRANS(T4,NU,2)
 1711 continue
      if(i.eq.j)then
         call symt411(t4,nu,6)
         goto 1710
      endif
      DO 710 M=1,NO
      IF(JK.AND.(KL.OR.K.EQ.M).OR.L.EQ.M.AND.(JL.OR.KL)) THEN
      CALL ZEROMA(O4(1,1,1,1,M),1,NU4)
      GOTO 710
      ENDIF
      CALL RDMIJKS(J,K,L,M,NO,NU,T4N,O4(1,1,1,1,M))
 710  CONTINUE
      CALL MTRANS(T4,NU,21)
      CALL MATMUL(O4,VO(1,1,1,I),T4,NU3,NU,NOU,0,1)
      CALL MTRANS(T4,NU,18)
      CALL MTRSM(O4,NO,NU,18)
      if(j.eq.k)call vecmul(t4,nu4,half)
      CALL MATMUL(O4,VO(1,1,1,I),T4,NU3,NU,NOU,0,1)
      if(j.eq.k)then
      CALL MTRSM(O4,NO,NU,10)
      CALL MTRANS(T4,NU,10)
         call symt411(t4,nu,7)
      else
      CALL MTRANS(T4,NU,2)
      CALL MTRSM(O4,NO,NU,2)
      if(k.eq.l)call vecmul(t4,nu4,half)
      CALL MATMUL(O4,VO(1,1,1,I),T4,NU3,NU,NOU,0,1)
      CALL MTRANS(T4,NU,1)
      CALL MTRSM(O4,NO,NU,1)
      endif
      if(l.eq.k)then
         call symt411(t4,nu,1)
      else
      CALL MATMUL(O4,VO(1,1,1,I),T4,NU3,NU,NOU,0,1)
      endif
      CALL MTRSM(O4,NO,NU,21)
      call smmijk(j,k,l,no,nu,t4n,o4)
      CALL MTRANS(T4,NU,21)
      CALL MATMUL(O4,VOE(1,1,1,I),T4,NU3,NU,NOU,0,0)
      CALL MTRANS(T4,NU,20)
 1710 continue
      call iexit(68)
      RETURN
      END
      subroutine smijkm(i,j,k,no,nu,t4,o4)
      implicit double precision (a-h,o-z)
      integer a,b,c,d
      dimension t4(nu,nu,nu,nu),o4(nu,nu,nu,nu,no)
      data two/2.0d+0/
      nu4=nu*nu*nu*nu
      do 10 m=1,no
      call veccop(nu4,t4,o4(1,1,1,1,m))
      call iperm4(i,j,k,m,iprm)
      goto (1,2,3,4)iprm
 1    continue
      do 11 d=1,nu
      do 11 c=1,nu
      do 11 b=1,nu
      do 11 a=1,nu
      o4(a,b,c,d,m)=t4(a,b,c,d)*two-t4(a,b,d,c)-t4(a,d,c,b)-t4(d,b,c,a)
 11   continue
      goto 1000
 2    continue
      call mtrans(t4,nu,1)
      do 12 d=1,nu
      do 12 c=1,nu
      do 12 a=1,nu
      do 12 b=1,nu
      o4(a,b,d,c,m)=t4(a,b,c,d)*two-t4(a,b,d,c)-t4(a,c,b,d)-t4(c,b,a,d)
 12   continue
      goto 1000
 3    continue
      call mtrans(t4,nu,2)
      do 13 a=1,nu
      do 13 b=1,nu
      do 13 c=1,nu
      do 13 d=1,nu
      o4(a,c,d,b,m)=t4(a,b,c,d)*two-t4(a,d,c,b)-t4(a,c,b,d)-t4(b,a,c,d)
 13   continue
      goto 1000
 4    continue
      call mtrans(t4,nu,20)
      do 14 a=1,nu
      do 14 b=1,nu
      do 14 c=1,nu
      do 14 d=1,nu
      o4(b,c,d,a,m)=t4(a,b,c,d)*two-t4(d,b,c,a)-t4(c,b,a,d)-t4(b,a,c,d)
 14   continue
 1000 continue
 10   continue
      return
      end
      subroutine smijmk(i,j,k,no,nu,t4,o4)
      implicit double precision (a-h,o-z)
      integer a,b,c,d
      dimension t4(nu,nu,nu,nu),o4(nu,nu,nu,nu,no)
      data two/2.0d+0/
      nu4=nu*nu*nu*nu
      do 10 m=1,no
      call veccop(nu4,t4,o4(1,1,1,1,m))
      call iperm4(i,j,k,m,iprm)
      goto (1,2,3,4)iprm
 1    continue
      do 11 a=1,nu
      do 11 b=1,nu
      do 11 c=1,nu
      do 11 d=1,nu
      o4(a,b,c,d,m)=t4(a,b,c,d)*two-t4(a,b,d,c)-t4(d,b,c,a)-t4(a,d,c,b)
 11   continue
      goto 1000
 2    continue
      call mtrans(t4,nu,1)
      do 12 a=1,nu
      do 12 b=1,nu
      do 12 c=1,nu
      do 12 d=1,nu
      o4(a,b,c,d,m)=t4(a,b,d,c)*two-t4(a,b,c,d)-t4(d,b,a,c)-t4(a,d,b,c)
 12   continue
      goto 1000
 3    continue
      call mtrans(t4,nu,2)
      do 13 a=1,nu
      do 13 b=1,nu
      do 13 c=1,nu
      do 13 d=1,nu
      o4(a,b,c,d,m)=t4(a,d,b,c)*two-t4(a,c,b,d)-t4(d,a,b,c)-t4(a,b,d,c)
 13   continue
      goto 1000
 4    continue
      call mtrans(t4,nu,20)
      do 14 a=1,nu
      do 14 b=1,nu
      do 14 c=1,nu
      do 14 d=1,nu
      o4(a,b,c,d,m)=t4(d,a,b,c)*two-t4(c,a,b,d)-t4(a,d,b,c)-t4(b,a,d,c)
 14   continue
 1000 continue
 10   continue
      return
      end
      subroutine smimjk(i,j,k,no,nu,t4,o4)
      implicit double precision (a-h,o-z)
      integer a,b,c,d
      dimension t4(nu,nu,nu,nu),o4(nu,nu,nu,nu,no)
      data two/2.0d+0/
      nu4=nu*nu*nu*nu
      do 10 m=1,no
      call veccop(nu4,t4,o4(1,1,1,1,m))
      call iperm4(i,j,k,m,iprm)
      goto (1,2,3,4)iprm
 1    continue
      do 11 a=1,nu
      do 11 b=1,nu
      do 11 c=1,nu
      do 11 d=1,nu
      o4(a,b,c,d,m)=t4(a,b,c,d)*two-t4(a,d,c,b)-t4(a,b,d,c)-t4(d,b,c,a)
 11   continue
      goto 1000
 2    continue
      call mtrans(t4,nu,1)
      do 12 a=1,nu
      do 12 b=1,nu
      do 12 c=1,nu
      do 12 d=1,nu
      o4(a,b,c,d,m)=t4(a,b,d,c)*two-t4(a,d,b,c)-t4(a,b,c,d)-t4(d,b,a,c)
 12   continue
      goto 1000
 3    continue
      call mtrans(t4,nu,2)
      do 13 a=1,nu
      do 13 b=1,nu
      do 13 c=1,nu
      do 13 d=1,nu
      o4(a,b,c,d,m)=t4(a,d,b,c)*two-t4(a,b,d,c)-t4(a,c,b,d)-t4(d,a,b,c)
 13   continue
      goto 1000
 4    continue
      call mtrans(t4,nu,20)
      do 14 a=1,nu
      do 14 b=1,nu
      do 14 c=1,nu
      do 14 d=1,nu
      o4(a,b,c,d,m)=t4(d,a,b,c)*two-t4(b,a,d,c)-t4(c,a,b,d)-t4(a,d,b,c)
 14   continue
 1000 continue
 10   continue
      return
      end
      subroutine smmijk(i,j,k,no,nu,t4,o4)
      implicit double precision (a-h,o-z)
      integer a,b,c,d
      dimension t4(nu,nu,nu,nu),o4(nu,nu,nu,nu,no)
      data two/2.0d+0/
      nu4=nu*nu*nu*nu
      do 10 m=1,no
      call veccop(nu4,t4,o4(1,1,1,1,m))
      call iperm4(i,j,k,m,iprm)
      goto (1,2,3,4)iprm
 1    continue
      do 11 a=1,nu
      do 11 b=1,nu
      do 11 c=1,nu
      do 11 d=1,nu
      o4(a,b,c,d,m)=t4(a,b,c,d)*two-t4(d,b,c,a)-t4(a,d,c,b)-t4(a,b,d,c)
 11   continue
      goto 1000
 2    continue
      call mtrans(t4,nu,1)
      do 12 a=1,nu
      do 12 b=1,nu
      do 12 c=1,nu
      do 12 d=1,nu
      o4(a,b,c,d,m)=t4(a,b,d,c)*two-t4(d,b,a,c)-t4(a,d,b,c)-t4(a,b,c,d)
 12   continue
      goto 1000
 3    continue
      call mtrans(t4,nu,2)
      do 13 a=1,nu
      do 13 b=1,nu
      do 13 c=1,nu
      do 13 d=1,nu
      o4(a,b,c,d,m)=t4(a,d,b,c)*two-t4(d,a,b,c)-t4(a,b,d,c)-t4(a,c,b,d)
 13   continue
      goto 1000
 4    continue
      call mtrans(t4,nu,20)
      do 14 a=1,nu
      do 14 b=1,nu
      do 14 c=1,nu
      do 14 d=1,nu
      o4(a,b,c,d,m)=t4(d,a,b,c)*two-t4(a,d,b,c)-t4(b,a,d,c)-t4(c,a,b,d)
 14   continue
 1000 continue
 10   continue
      return
      end
      subroutine drt4ppt2(i,j,k,l,no,nu,o1)
      implicit double precision (a-h,o-z)
      logical print
      common /flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension o1(1)
      common/itrat/icycle,mx,icn
      common/totmem/mem
      print=iflags(1).gt.10
      i1=1          !t2new
      i2=i1+no2u2   !t4
      i3=i2+nu4     !v
      i4=i3+nu3     !t2
      it=i4+no2u2 
      if(icycle.eq.1.and.i.eq.2.and.print)write(6,98)mem,it
 98   format('Space usage in t4ppt2 available - ',i10,'   used - ',i10)
      call t4ppt2(i,j,k,l,no,nu,o1(i3),o1(i4),o1(i2))
      return
      end
      SUBROUTINE T4PPT2(I,J,K,L,NO,NU,V,T2,T4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer a,b,c,d,e
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION T2(NU,NU,NO,NO),V(NU,NU,NU),T4(NU,NU,NU,NU)
      data zero/0.0d+0/
      call ienter(69)
      CALL RO2PPH(0,NO,NU,V,T2)
c
      kk=it2(k,l)
      IAS=(kk-1)*nu+1+no3
c      IAS=(K-1)*NOU+(L-1)*NU+NO3+1
      DO 190 A=1,NU
      IASV=IAS+A
      CALL RDVT3(IASV,NU,V)
      CALL MATMUL(T2(1,1,i,j),V,T4(1,1,1,A),NU,NU2,NU,0,0)
c      CALL TRANMD(T4(1,1,1,A),NU,NU,NU,1,12)
      CALL TRANT3(T4(1,1,1,A),NU,2)
      CALL MATMUL(T2(1,1,J,I),V,T4(1,1,1,A),NU,NU2,NU,0,0)
 190  CONTINUE
      call mtrans(t4,nu,18)
c
      kk=it2(i,j)
      IAS=(kk-1)*nu+1+no3
      DO 191 A=1,NU
      IASV=IAS+A
      CALL RDVT3(IASV,NU,V)
      CALL MATMUL(T2(1,1,l,k),V,T4(1,1,1,A),NU,NU2,NU,0,0)
      CALL TRANT3(T4(1,1,1,A),NU,2)
      CALL MATMUL(T2(1,1,K,L),V,T4(1,1,1,A),NU,NU2,NU,0,0)
 191  CONTINUE
      CALL MTRANS(T4,NU,16)
c
      kk=it2(i,k)
      IAS=(kk-1)*nu+1+no3
      DO 192 A=1,NU
      IASV=IAS+A
      CALL RDVT3(IASV,NU,V)
      CALL MATMUL(T2(1,1,l,j),V,T4(1,1,1,A),NU,NU2,NU,0,0)
      CALL TRANT3(T4(1,1,1,A),NU,2)
      CALL MATMUL(T2(1,1,J,L),V,T4(1,1,1,A),NU,NU2,NU,0,0)
 192  CONTINUE
      CALL MTRANS(T4,NU,15)
c
      kk=it2(i,l)
      IAS=(kk-1)*nu+1+no3
      DO 193 A=1,NU
      IASV=IAS+A
      CALL RDVT3(IASV,NU,V)
      CALL MATMUL(T2(1,1,k,j),V,T4(1,1,1,A),NU,NU2,NU,0,0)
      CALL TRANT3(T4(1,1,1,A),NU,2)
      CALL MATMUL(T2(1,1,J,K),V,T4(1,1,1,A),NU,NU2,NU,0,0)
 193  continue
c
      kk=it2(j,k)
      IAS=(kk-1)*nu+1+no3
      call mtrans(t4,nu,4)
      DO 194 A=1,NU
      IASV=IAS+A
      CALL RDVT3(IASV,NU,V)
      CALL MATMUL(T2(1,1,l,i),V,T4(1,1,1,A),NU,NU2,NU,0,0)
      CALL TRANT3(T4(1,1,1,A),NU,2)
      CALL MATMUL(T2(1,1,I,L),V,T4(1,1,1,A),NU,NU2,NU,0,0)
 194  CONTINUE
      CALL MTRANS(T4,NU,15)
c
      kk=it2(j,l)
      IAS=(kk-1)*nu+1+no3
      DO 195 A=1,NU
      IASV=IAS+A
      CALL RDVT3(IASV,NU,V)
      CALL MATMUL(T2(1,1,K,i),V,T4(1,1,1,A),NU,NU2,NU,0,0)
      CALL TRANT3(T4(1,1,1,A),NU,2)
      CALL MATMUL(T2(1,1,I,K),V,T4(1,1,1,A),NU,NU2,NU,0,0)
 195  CONTINUE
      call mtrans(t4,nu,7)
      call iexit(69)
      RETURN
      END
      subroutine drt4ppt3(i,j,k,l,no,nu,o1)
      implicit double precision (a-h,o-z)
      logical print
      common /flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension o1(1)
      common/itrat/icycle,mx,icn
      common/totmem/mem
      print=iflags(1).gt.10
      i1=1          !t2new
      i2=i1+no2u2   !t4
      i3=i2+nu4     !vt3m
      i4=i3+nou3    !v
      i5=i4+nou3    !vm
      it=i5+no3u
      if(icycle.eq.2.and.i.eq.2.and.print)write(6,98)mem,it
 98   format('Space usage in t4ppt3 available - ',i10,'   used - ',i10)
      call t4ppt3(i,j,k,l,no,nu,o1(i3),o1(i4),o1(i5),o1(i2))
      return
      end
      SUBROUTINE T4PPT3(I,J,K,L,NO,NU,VT3M,V,VM,T4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SDTQ2,SDTQ4
      COMMON/NEWOPT/NOPT(6)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION VT3M(NU,NU,NU,NO),V(NU,NU,NU,NO),VM(NO,NU,NO,NO),
     *T4(NU,NU,NU,NU)
      call ienter(70)
      SDTQ2=NOPT(1).GT.3
      SDTQ4=NOPT(1).GT.5
      IF (.NOT.SDTQ2)THEN
      CALL RDVEM4(0,NO,NU,VT3M,V)
      CALL RDVEM2(1,NU,NO,VT3M,VM)
      ELSE
      IF(.NOT.SDTQ4)THEN
      CALL RDVEM4A(0,NO,NU,VT3M,V)
      CALL RDVEM2A(1,NU,NO,VT3M,VM)
      ELSE
      CALL RDVEM4A(10,NO,NU,VT3M,V)
      CALL RDVEM2A(11,NU,NO,VT3M,VM)
      ENDIF 
      ENDIF
      CALL TRANMD(VM,NO,NU,NO,NO,13)
      CALL RDVT3ONW(I,J,K,NU,VT3M)
      CALL MATMUL(VT3M,V(1,1,1,L),T4,NU2,NU2,NU,0,0)
      CALL MTRANS(T4,NU,1)
      CALL RDVT3ONW(I,J,L,NU,VT3M)
      CALL MATMUL(VT3M,V(1,1,1,K),T4,NU2,NU2,NU,0,0)
      CALL MTRANS(T4,NU,2)
      CALL RDVT3ONW(I,K,J,NU,VT3M)
      CALL MATMUL(VT3M,V(1,1,1,L),T4,NU2,NU2,NU,0,0)
      CALL MTRANS(T4,NU,1)
      CALL RDVT3ONW(I,K,L,NU,VT3M)
      CALL MATMUL(VT3M,V(1,1,1,J),T4,NU2,NU2,NU,0,0)
      CALL MTRANS(T4,NU,3)
      CALL RDVT3ONW(I,L,J,NU,VT3M)
      CALL MATMUL(VT3M,V(1,1,1,K),T4,NU2,NU2,NU,0,0)
      CALL MTRANS(T4,NU,1)
      CALL RDVT3ONW(I,L,K,NU,VT3M)
      CALL MATMUL(VT3M,V(1,1,1,J),T4,NU2,NU2,NU,0,0)
      CALL MTRANS(T4,NU,4)
      CALL RDVT3ONW(J,K,I,NU,VT3M)
      CALL MATMUL(VT3M,V(1,1,1,L),T4,NU2,NU2,NU,0,0)
      CALL MTRANS(T4,NU,1)
      CALL RDVT3ONW(J,K,L,NU,VT3M)
      CALL MATMUL(VT3M,V(1,1,1,I),T4,NU2,NU2,NU,0,0)
      CALL MTRANS(T4,NU,3)
      CALL RDVT3ONW(J,L,I,NU,VT3M)
      CALL MATMUL(VT3M,V(1,1,1,K),T4,NU2,NU2,NU,0,0)
      CALL MTRANS(T4,NU,1)
      CALL RDVT3ONW(J,L,K,NU,VT3M)
      CALL MATMUL(VT3M,V(1,1,1,I),T4,NU2,NU2,NU,0,0)
      CALL MTRANS(T4,NU,5)
      CALL RDVT3ONW(K,L,I,NU,VT3M)
      CALL MATMUL(VT3M,V(1,1,1,J),T4,NU2,NU2,NU,0,0)
      CALL MTRANS(T4,NU,1)
      CALL RDVT3ONW(K,L,J,NU,VT3M)
      CALL MATMUL(VT3M,V(1,1,1,I),T4,NU2,NU2,NU,0,0)
      CALL MTRANS(T4,NU,4)
      DO 65 M=1,NO
      IF (M.EQ.I.AND.M.EQ.J)THEN
      CALL ZEROMA(V(1,1,1,M),1,NU3)
      GO TO 65
      ENDIF
      CALL RDVT3ONW(I,J,M,NU,V(1,1,1,M))
   65 CONTINUE
      CALL MATMUL(V,VM(1,1,L,K),T4,NU3,NU,NO,0,1)
      CALL MTRANS(T4,NU,1)
      CALL MATMUL(V,VM(1,1,K,L),T4,NU3,NU,NO,0,1)
      CALL MTRANS(T4,NU,2)
      DO 66 M=1,NO
      IF (M.EQ.I.AND.M.EQ.K)THEN
      CALL ZEROMA(V(1,1,1,M),1,NU3)
      GO TO 66
      ENDIF

      CALL RDVT3ONW(I,K,M,NU,V(1,1,1,M))
   66 CONTINUE
      CALL MATMUL(V,VM(1,1,L,J),T4,NU3,NU,NO,0,1)
      CALL MTRANS(T4,NU,1)
      CALL MATMUL(V,VM(1,1,J,L),T4,NU3,NU,NO,0,1)
      CALL MTRANS(T4,NU,3)
      DO 67 M=1,NO
      IF (M.EQ.I.AND.M.EQ.L)THEN
      CALL ZEROMA(V(1,1,1,M),1,NU3)
      GO TO 67
      ENDIF
      CALL RDVT3ONW(I,L,M,NU,V(1,1,1,M))
   67 CONTINUE
      CALL MATMUL(V,VM(1,1,K,J),T4,NU3,NU,NO,0,1)
      CALL MTRANS(T4,NU,1)
      CALL MATMUL(V,VM(1,1,J,K),T4,NU3,NU,NO,0,1)
      CALL MTRANS(T4,NU,4)
      DO 68 M=1,NO
      IF (M.EQ.J.AND.M.EQ.K)THEN
      CALL ZEROMA(V(1,1,1,M),1,NU3)
      GO TO 68
      ENDIF
      CALL RDVT3ONW(J,K,M,NU,V(1,1,1,M))
   68 CONTINUE
      CALL MATMUL(V,VM(1,1,L,I),T4,NU3,NU,NO,0,1)
      CALL MTRANS(T4,NU,1)
      CALL MATMUL(V,VM(1,1,I,L),T4,NU3,NU,NO,0,1)
      CALL MTRANS(T4,NU,3)
      DO 69 M=1,NO
      IF (M.EQ.L.AND.M.EQ.J)THEN
      CALL ZEROMA(V(1,1,1,M),1,NU3)
      GO TO 69
      ENDIF
      CALL RDVT3ONW(J,L,M,NU,V(1,1,1,M))
   69 CONTINUE
      CALL MATMUL(V,VM(1,1,K,I),T4,NU3,NU,NO,0,1)
      CALL MTRANS(T4,NU,1)
      CALL MATMUL(V,VM(1,1,I,K),T4,NU3,NU,NO,0,1)
      CALL MTRANS(T4,NU,5)
      DO 70 M=1,NO
      IF (M.EQ.K.AND.M.EQ.L)THEN
      CALL ZEROMA(V(1,1,1,M),1,NU3)
      GO TO 70
      ENDIF
      CALL RDVT3ONW(K,L,M,NU,V(1,1,1,M))
   70 CONTINUE
      CALL MATMUL(V,VM(1,1,J,I),T4,NU3,NU,NO,0,1)
      CALL MTRANS(T4,NU,1)
      CALL MATMUL(V,VM(1,1,I,J),T4,NU3,NU,NO,0,1)
      CALL MTRANS(T4,NU,4)
      call iexit(70)
      RETURN
      END
      subroutine drt4ppt4(i,j,k,l,no,nu,o1)
      implicit double precision (a-h,o-z)
      logical print
      common /flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension o1(1)
      common/itrat/icycle,mx,icn
      common/totmem/mem
      print=iflags(1).gt.10
      i1=1          !t2
      i2=i1+no2u2   !t4
      i3=i2+nu4     !o4
      i4=i3+nu4     !v
      it=i4+nu4
      if(icycle.eq.2.and.i.eq.2.and.print)write(6,98)mem,it
 98   format('Space usage in t4ppt4: available - ',i10,'   used - ',i10)
      call t4ppt4(i,j,k,l,no,nu,o1(i3),o1(i2),o1(i4))
      return
      end
      SUBROUTINE T4PPT4(I,J,K,L,NO,NU,O4,T4,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D,E,F
      COMMON/NEWOPT/NOPT(6)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION O4(NU,NU,NU,NU),T4(NU,NU,NU,NU),V(NU,NU,NU,NU)
      DATA ZERO/0.0D+0/,HALF/0.5D+0/,FOURTH/0.25D+0/,
     *THIRD/0.3333333333333333333D+0/
      call ienter(71)
      NU4=NU3*NU
      I10=0
      IF(NOPT(1).GT.5)I10=10
      CALL RDOV4(I10,NO,NU,O4,V)
      IAS=IT4(I,J,K,L)
      CALL RDT4(IAS,NU,O4)
      IF (I10.EQ.0)CALL MTRANS(V,NU,7)
      if(i.eq.j.and.k.eq.l)goto 1201
      if(i.eq.j)goto 301
      if(j.eq.k)goto 401
      if(k.eq.l)goto 501
       CALL MATMUL(V,O4,T4,NU2,NU2,NU2,0,0)
      CALL MTRANS(O4,NU,7)
      CALL MTRANS(T4,NU,7)
       CALL MATMUL(V,O4,T4,NU2,NU2,NU2,0,0)
      CALL MTRANS(O4,NU,10)
      CALL MTRANS(T4,NU,10)
      CALL MATMUL(V,O4,T4,NU2,NU2,NU2,0,0)
      CALL MTRANS(V,NU,9)
      CALL MATMUL(O4,V,T4,NU2,NU2,NU2,0,0)
      CALL MTRANS(O4,NU,10)
      CALL MTRANS(T4,NU,10)
      CALL MATMUL(O4,V,T4,NU2,NU2,NU2,0,0)
      CALL MTRANS(O4,NU,7)
      CALL MTRANS(T4,NU,7)
      CALL MATMUL(O4,V,T4,NU2,NU2,NU2,0,0)
      goto 999
 201  continue
      CALL MATMUL(V,O4,T4,NU2,NU2,NU2,0,0)
      CALL MTRANS(O4,NU,7)
      CALL MTRANS(T4,NU,7)
      call vecmul(t4,nu4,half)
      CALL MATMUL(V,O4,T4,NU2,NU2,NU2,0,0)
      CALL MTRANS(O4,NU,10)
      CALL MTRANS(T4,NU,7)
      call symt411(t4,nu,1)
      CALL MTRANS(T4,NU,7)
      CALL MTRANS(V,NU,9)
      call vecmul(t4,nu4,half)
      CALL MATMUL(O4,V,T4,NU2,NU2,NU2,0,0)
      CALL MTRANS(O4,NU,3)
      call symt411(t4,nu,10)
      CALL MTRANS(T4,NU,7)
      CALL MATMUL(O4,V,T4,NU2,NU2,NU2,0,0)
      goto 999
 301  continue
      CALL MATMUL(V,O4,T4,NU2,NU2,NU2,0,0)
      CALL MTRANS(O4,NU,7)
      CALL MTRANS(T4,NU,7)
      call vecmul(t4,nu4,half)
      CALL MATMUL(V,O4,T4,NU2,NU2,NU2,0,0)
      call symt411(t4,nu,8)
      CALL MTRANS(O4,NU,10)
      CALL MTRANS(T4,NU,9)
      call vecmul(t4,nu4,half)
      CALL MATMUL(V,O4,T4,NU2,NU2,NU2,0,0)
      call symt411(t4,nu,8)
      CALL MTRANS(O4,NU,14)
      CALL MTRANS(T4,NU,14)
      CALL MATMUL(V,O4,T4,NU2,NU2,NU2,0,0)
      CALL MTRANS(O4,NU,9)
      CALL MTRANS(T4,NU,9)
      goto 999
 401  continue
      call vecmul(t4,nu4,half)
      CALL MATMUL(V,O4,T4,NU2,NU2,NU2,0,0)
      call symt411(t4,nu,7)
      CALL MTRANS(O4,NU,2)
      CALL MTRANS(T4,NU,10)
      CALL MATMUL(V,O4,T4,NU2,NU2,NU2,0,0)
      CALL MTRANS(O4,NU,9)
      CALL MTRANS(T4,NU,9)
      CALL MATMUL(V,O4,T4,NU2,NU2,NU2,0,0)
      CALL MTRANS(O4,NU,10)
      CALL MTRANS(T4,NU,10)
      call vecmul(t4,nu4,half)
      CALL MATMUL(V,O4,T4,NU2,NU2,NU2,0,0)
      call symt411(t4,nu,11)
      CALL MTRANS(O4,NU,22)
      CALL MTRANS(T4,NU,9)
      goto 999
 501  continue
      CALL MATMUL(V,O4,T4,NU2,NU2,NU2,0,0) 
      CALL MTRANS(O4,NU,7)
      CALL MTRANS(T4,NU,7)
      call vecmul(t4,nu4,half)
      CALL MATMUL(V,O4,T4,NU2,NU2,NU2,0,0)  
      CALL MTRANS(O4,NU,8)
      call symt411(t4,nu,10)
      CALL MTRANS(T4,NU,9)
      call vecmul(t4,nu4,half)
      CALL MATMUL(V,O4,T4,NU2,NU2,NU2,0,0)  
      call symt411(t4,nu,10)
      CALL MTRANS(O4,NU,16)
      CALL MTRANS(T4,NU,11)
      CALL MATMUL(V,O4,T4,NU2,NU2,NU2,0,0)  
      CALL MTRANS(O4,NU,9)
      CALL MTRANS(T4,NU,9)
      goto 999
1201  continue
      CALL MTRANS(O4,NU,7)
      CALL MTRANS(T4,NU,7)
      call vecmul(t4,nu4,half)
      CALL MATMUL(V,O4,T4,NU2,NU2,NU2,0,0)  
      call symt411(t4,nu,8)
      CALL MTRANS(O4,NU,7)
      CALL MTRANS(T4,NU,7)
      CALL MATMUL(V,O4,T4,NU2,NU2,NU2,0,0)  
      CALL MTRANS(O4,NU,7)
      CALL MTRANS(T4,NU,7)
      CALL MTRANS(O4,NU,10)
      CALL MTRANS(T4,NU,9)
      call vecmul(t4,nu4,half)
      CALL MATMUL(V,O4,T4,NU2,NU2,NU2,0,0)  
      call symt411(t4,nu,8)
      CALL MTRANS(O4,NU,14)
      CALL MTRANS(T4,NU,14)
      CALL MATMUL(V,O4,T4,NU2,NU2,NU2,0,0)  
      CALL MTRANS(O4,NU,9)
      CALL MTRANS(T4,NU,9)
      goto 999
 999  continue
      call iexit(71)
      RETURN
      END
      subroutine drt4vot2(i,j,k,l,no,nu,o1)
      implicit double precision (a-h,o-z)
      logical print
      common /flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension o1(1)
      common/itrat/icycle,mx,icn
      common/totmem/mem
      print=iflags(1).gt.10
      i1=1          !t2new
      i2=i1+no2u2   !t4
      i3=i2+nu4     !vt
      i4=i3+nou2    !t2
      it=i4+no2u2
      if(icycle.eq.1.and.i.eq.2.and.print)write(6,98)mem,it
 98   format('Space usage in t4vot2: available - ',i10,'   used - ',i10)
      call t4vot2(i,j,k,l,no,nu,o1(i3),o1(i4),o1(i2))
      return
      end
      SUBROUTINE T4VOT2(I,J,K,L,NO,NU,VT,T2,T4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NVTAP/NTITER,NQ2,NVT
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION VT(1),T2(NO,NU,NU,NO),T4(1)
      call ienter(73)
      CALL RO2HPP(0,NO,NU,VT,T2)
      KKK=IGN(I,J,K,NO)
      CALL RDGEN(NVT,KKK,NOU2,VT)
      CALL MATMUL(VT,T2(1,1,1,L),T4,NU2,NU2,NO,0,1)
      CALL MTRANS(T4,NU,1)
      KKK=IGN(I,J,L,NO)
      CALL RDGEN(NVT,KKK,NOU2,VT)
      CALL MATMUL(VT,T2(1,1,1,K),T4,NU2,NU2,NO,0,1)
      CALL MTRANS(T4,NU,2)
      KKK=IGN(I,K,J,NO)
      CALL RDGEN(NVT,KKK,NOU2,VT)
      CALL MATMUL(VT,T2(1,1,1,L),T4,NU2,NU2,NO,0,1)
      CALL MTRANS(T4,NU,1)
      KKK=IGN(I,K,L,NO)
      CALL RDGEN(NVT,KKK,NOU2,VT)
      CALL MATMUL(VT,T2(1,1,1,J),T4,NU2,NU2,NO,0,1)
      CALL MTRANS(T4,NU,3)
      KKK=IGN(I,L,J,NO)
      CALL RDGEN(NVT,KKK,NOU2,VT)
      CALL MATMUL(VT,T2(1,1,1,K),T4,NU2,NU2,NO,0,1)
      CALL MTRANS(T4,NU,1)
      KKK=IGN(I,L,K,NO)
      CALL RDGEN(NVT,KKK,NOU2,VT)
      CALL MATMUL(VT,T2(1,1,1,J),T4,NU2,NU2,NO,0,1)
      CALL MTRANS(T4,NU,4)
      KKK=IGN(J,K,I,NO)
      CALL RDGEN(NVT,KKK,NOU2,VT)
      CALL MATMUL(VT,T2(1,1,1,L),T4,NU2,NU2,NO,0,1)
      CALL MTRANS(T4,NU,1)
      KKK=IGN(J,K,L,NO)
      CALL RDGEN(NVT,KKK,NOU2,VT)
      CALL MATMUL(VT,T2(1,1,1,I),T4,NU2,NU2,NO,0,1)
      CALL MTRANS(T4,NU,3)
      KKK=IGN(J,L,I,NO)
      CALL RDGEN(NVT,KKK,NOU2,VT)
      CALL MATMUL(VT,T2(1,1,1,K),T4,NU2,NU2,NO,0,1)
      CALL MTRANS(T4,NU,1)
      KKK=IGN(J,L,K,NO)
      CALL RDGEN(NVT,KKK,NOU2,VT)
      CALL MATMUL(VT,T2(1,1,1,I),T4,NU2,NU2,NO,0,1)
      CALL MTRANS(T4,NU,5)
      KKK=IGN(K,L,I,NO)
      CALL RDGEN(NVT,KKK,NOU2,VT)
      CALL MATMUL(VT,T2(1,1,1,J),T4,NU2,NU2,NO,0,1)
      CALL MTRANS(T4,NU,1)
      KKK=IGN(K,L,J,NO)
      CALL RDGEN(NVT,KKK,NOU2,VT)
      CALL MATMUL(VT,T2(1,1,1,I),T4,NU2,NU2,NO,0,1)
      CALL MTRANS(T4,NU,4)
      call iexit(73)
      RETURN
      END
      subroutine drt4wt32(i,j,k,l,no,nu,t)
      implicit double precision (a-h,o-z)
      logical print
      common /flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/itrat/icycle,mx,icn
      common/totmem/mem
      dimension t(*)
      print=iflags(1).gt.10
      no3u3=no3*nu3
      i1=1            !t2
      i2=i1+no2u2     !t4
      i3=i2+nu4       !ti
      i4=i3+nou3      !t3
      i5=i4+no3u3     !t3i
      it=i5+no3u3
      if(icycle.eq.10.and.i.eq.2.and.print)write(6,98)mem,it
 98   format('Space usage in t4wt32me: available - ',i8,'   used - ',i8)
      call t4wt32me(i,j,k,l,no,nu,t(i3),t(i4),t(i5),t(i2))
      i1=1            !t2
      i2=i1+no2u2     !t4
      i3=i2+nu4       !ti
      i4=i3+nu3       !t3
      i5=i4+no3u3     !t3mn
      it=i5+no4*nou
      if(icycle.eq.10.and.i.eq.2.and.print)write(6,99)mem,it
 99   format('Space usage in t4wt32mn: available - ',i8,'   used - ',i8)
      call t4wt32mn(i,j,k,l,no,nu,t(i3),t(i4),t(i5),t(i2))
      return
      end
      SUBROUTINE T4WT32ME(I,J,K,L,NO,NU,TI,T3,T3I,T4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D,E,F
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION TI(NU,NU,NU),T3(NO,NU,NU,NU,NO,NO),T4(NU,NU,NU,NU),
     *          T3I(NU,NU,NO,NU,NO,NO),TSYS(2)
      DATA ZERO/0.0D+0/,TWO/2.0D+0/,HALF/0.5D+0/
      call ienter(74)
      DO 10 I1=1,NO
      DO 10 J1=1,NO
      DO 10 K1=1,NO
      IAS=NO2*(I1-1)+NO*(J1-1)+K1
      CALL RDVT3I(IAS,NU,TI)
      DO 8 C=1,NU
         call veccop(nu2,t3i(1,1,k1,c,i1,j1),ti(1,1,c))
 8    CONTINUE
 10   CONTINUE
      if(i.eq.j)goto 901
      if(j.eq.k)goto 902
      if(k.eq.l)goto 903
      CALL RDT32(K,L,NO,NU,TI,T3)
      CALL MATMUL(T3I(1,1,1,1,I,J),T3,T4,NU2,NU2,NOU,0,0)
      CALL MTRANS(T4,NU,7)
      CALL RDT32(J,L,NO,NU,TI,T3)
      CALL MATMUL(T3I(1,1,1,1,I,K),T3,T4,NU2,NU2,NOU,0,0)
      CALL MTRANS(T4,NU,2)
      CALL RDT32(K,J,NO,NU,TI,T3)
      CALL MATMUL(T3I(1,1,1,1,I,L),T3,T4,NU2,NU2,NOU,0,0)
      CALL RDT32(I,L,NO,NU,TI,T3)
      CALL MTRANS(T4,NU,4)
      CALL MATMUL(T3I(1,1,1,1,J,K),T3,T4,NU2,NU2,NOU,0,0)
      CALL MTRANS(T4,NU,2)
      CALL RDT32(K,I,NO,NU,TI,T3)
      CALL MATMUL(T3I(1,1,1,1,J,L),T3,T4,NU2,NU2,NOU,0,0)
      CALL MTRANS(T4,NU,5)
      CALL RDT32(I,J,NO,NU,TI,T3)
      CALL MATMUL(T3I(1,1,1,1,K,L),T3,T4,NU2,NU2,NOU,0,0)
      CALL MTRANS(T4,NU,9)
      goto 910
 901  continue
      CALL RDT32(K,L,NO,NU,TI,T3)
      CALL MATMUL(T3I(1,1,1,1,I,J),T3,T4,NU2,NU2,NOU,0,0)
      CALL MTRANS(T4,NU,7)
      CALL RDT32(J,L,NO,NU,TI,T3)
      call vecmul(t4,nu4,half)
      CALL MATMUL(T3I(1,1,1,1,I,K),T3,T4,NU2,NU2,NOU,0,0)
      CALL MTRANS(T4,NU,2)
      CALL RDT32(K,J,NO,NU,TI,T3)
      CALL MATMUL(T3I(1,1,1,1,I,L),T3,T4,NU2,NU2,NOU,0,0)
      CALL MTRANS(T4,NU,10)
      call symt411(t4,nu,6)
      CALL MTRANS(T4,NU,9)
      CALL RDT32(I,J,NO,NU,TI,T3)
      CALL MATMUL(T3I(1,1,1,1,K,L),T3,T4,NU2,NU2,NOU,0,0)
      CALL MTRANS(T4,NU,9)
      goto 910
 902  continue
      CALL RDT32(K,L,NO,NU,TI,T3)
      call vecmul(t4,nu4,half)
      CALL MATMUL(T3I(1,1,1,1,I,J),T3,T4,NU2,NU2,NOU,0,0)
      CALL MTRANS(T4,NU,16)
      CALL RDT32(K,I,NO,NU,TI,T3)
      CALL MATMUL(T3I(1,1,1,1,J,L),T3,T4,NU2,NU2,NOU,0,0)
      CALL MTRANS(T4,NU,15)
      call symt411(t4,nu,7)
      CALL MTRANS(T4,NU,10)
      CALL RDT32(K,J,NO,NU,TI,T3)
      CALL MATMUL(T3I(1,1,1,1,I,L),T3,T4,NU2,NU2,NOU,0,0)
      CALL MTRANS(T4,NU,4)
      CALL RDT32(I,L,NO,NU,TI,T3)
      CALL MATMUL(T3I(1,1,1,1,J,K),T3,T4,NU2,NU2,NOU,0,0)
      CALL MTRANS(T4,NU,12)
      goto 910
 903  continue
      CALL RDT32(K,L,NO,NU,TI,T3)
      CALL MATMUL(T3I(1,1,1,1,I,J),T3,T4,NU2,NU2,NOU,0,0)
      CALL MTRANS(T4,NU,7)
      CALL RDT32(J,L,NO,NU,TI,T3)
      call vecmul(t4,nu4,half)
      CALL MATMUL(T3I(1,1,1,1,I,K),T3,T4,NU2,NU2,NOU,0,0)
      CALL MTRANS(T4,NU,8)
      CALL RDT32(I,L,NO,NU,TI,T3)
      CALL MATMUL(T3I(1,1,1,1,J,K),T3,T4,NU2,NU2,NOU,0,0)
      CALL MTRANS(T4,NU,12)
      call symt411(t4,nu,1)
      CALL MTRANS(T4,NU,9)
      CALL RDT32(I,J,NO,NU,TI,T3)
      CALL MATMUL(T3I(1,1,1,1,K,L),T3,T4,NU2,NU2,NOU,0,0)
      CALL MTRANS(T4,NU,9)
      goto 910
 910  continue
      DO 110 I1=1,NO
      DO 110 J1=1,NO
      DO 110 K1=1,NO
      IAS=NO3+NO2*(I1-1)+NO*(J1-1)+K1
      CALL RDVT3I(IAS,NU,TI)
      DO 18 C=1,NU
         call veccop(nu2,t3i(1,1,k1,c,i1,j1),ti(1,1,c))
 18   CONTINUE
 110  CONTINUE
      DO 103 I1=1,NO
      DO 103 J1=1,NO
      DO 103 K1=1,NO
      IF(I1.EQ.J1.AND.J1.EQ.K1)THEN
      CALL ZEROMA(TI,1,NU3)
      GO TO 39
      ENDIF
      CALL RDVT3ONW(I1,J1,K1,NU,TI)
 39   CONTINUE
      DO 113 C=1,NU
      DO 113 B=1,NU
      DO 113 A=1,NU
      T3(I1,A,B,C,J1,K1)=TI(A,B,C)
 113  CONTINUE
 103  CONTINUE
c
      if(i.eq.j)goto 1001
      if(j.eq.k)goto 1002
      if(k.eq.l)goto 1003
      CALL TRT3ALL(NO,NU,T3,23)
      CALL MATMUL(T3I(1,1,1,1,I,K),T3(1,1,1,1,J,L),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,7)
      CALL MATMUL(T3I(1,1,1,1,I,J),T3(1,1,1,1,K,L),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,10)
      CALL MATMUL(T3I(1,1,1,1,I,J),T3(1,1,1,1,L,K),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,1)
      CALL MATMUL(T3I(1,1,1,1,I,K),T3(1,1,1,1,L,J),T4,NU2,NU2,NOU,0,1)
C
      CALL MTRANS(T4,NU,7)
      CALL MATMUL(T3I(1,1,1,1,I,L),T3(1,1,1,1,K,J),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,10)
      CALL MATMUL(T3I(1,1,1,1,I,L),T3(1,1,1,1,J,K),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,16)
      CALL MATMUL(T3I(1,1,1,1,J,L),T3(1,1,1,1,K,I),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,1)
      CALL MATMUL(T3I(1,1,1,1,J,I),T3(1,1,1,1,K,L),T4,NU2,NU2,NOU,0,1)
c
      CALL MTRANS(T4,NU,10)
      CALL MATMUL(T3I(1,1,1,1,J,I),T3(1,1,1,1,L,K),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,1)
      CALL MATMUL(T3I(1,1,1,1,J,K),T3(1,1,1,1,L,I),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,10)
      CALL MATMUL(T3I(1,1,1,1,J,K),T3(1,1,1,1,I,L),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,1)
      CALL MATMUL(T3I(1,1,1,1,J,L),T3(1,1,1,1,I,K),T4,NU2,NU2,NOU,0,1)
C
       CALL MTRANS(T4,NU,11)
      CALL MATMUL(T3I(1,1,1,1,K,L),T3(1,1,1,1,I,J),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,1)
      CALL MATMUL(T3I(1,1,1,1,K,J),T3(1,1,1,1,I,L),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,7)
      CALL MATMUL(T3I(1,1,1,1,K,I),T3(1,1,1,1,J,L),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,10)
      CALL MATMUL(T3I(1,1,1,1,K,I),T3(1,1,1,1,L,J),T4,NU2,NU2,NOU,0,1)
C
      CALL MTRANS(T4,NU,1)
      CALL MATMUL(T3I(1,1,1,1,K,J),T3(1,1,1,1,L,I),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,7)
      CALL MATMUL(T3I(1,1,1,1,K,L),T3(1,1,1,1,J,I),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,8)
      CALL MATMUL(T3I(1,1,1,1,L,K),T3(1,1,1,1,J,I),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,1)
      CALL MATMUL(T3I(1,1,1,1,L,I),T3(1,1,1,1,J,K),T4,NU2,NU2,NOU,0,1)
C
      CALL MTRANS(T4,NU,10)
      CALL MATMUL(T3I(1,1,1,1,L,I),T3(1,1,1,1,K,J),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,1)
      CALL MATMUL(T3I(1,1,1,1,L,J),T3(1,1,1,1,K,I),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,10)
      CALL MATMUL(T3I(1,1,1,1,L,J),T3(1,1,1,1,I,K),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,1)
      CALL MATMUL(T3I(1,1,1,1,L,K),T3(1,1,1,1,I,J),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,16)
      CALL TRT3ALL(NO,NU,T3,23)
      goto 9000
 1001 continue
      CALL TRT3ALL(NO,NU,T3,23)
      call vecmul(t4,nu4,half)
      CALL MATMUL(T3I(1,1,1,1,I,K),T3(1,1,1,1,J,L),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,10)
      CALL MATMUL(T3I(1,1,1,1,I,K),T3(1,1,1,1,L,J),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,8)
      CALL MATMUL(T3I(1,1,1,1,K,I),T3(1,1,1,1,L,J),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,10)
      CALL MATMUL(T3I(1,1,1,1,K,I),T3(1,1,1,1,J,L),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,22)
      CALL MATMUL(T3I(1,1,1,1,I,L),T3(1,1,1,1,K,J),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,10)
      CALL MATMUL(T3I(1,1,1,1,I,L),T3(1,1,1,1,J,K),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,8)
      CALL MATMUL(T3I(1,1,1,1,L,I),T3(1,1,1,1,J,K),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,10)
      CALL MATMUL(T3I(1,1,1,1,L,I),T3(1,1,1,1,K,J),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,5)
      CALL MATMUL(T3I(1,1,1,1,I,J),T3(1,1,1,1,K,L),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,10)
      CALL MATMUL(T3I(1,1,1,1,I,J),T3(1,1,1,1,L,K),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,18)
      CALL MATMUL(T3I(1,1,1,1,K,L),T3(1,1,1,1,J,I),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,8)
      CALL MATMUL(T3I(1,1,1,1,L,K),T3(1,1,1,1,J,I),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,11)
      call symt411(t4,nu,6)
      CALL TRT3ALL(NO,NU,T3,23)
      goto 9000
 1002 continue
      CALL TRT3ALL(NO,NU,T3,23)
      call vecmul(t4,nu4,half)
      CALL MATMUL(T3I(1,1,1,1,I,K),T3(1,1,1,1,J,L),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,10)
      CALL MATMUL(T3I(1,1,1,1,I,K),T3(1,1,1,1,L,J),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,8)
      CALL MATMUL(T3I(1,1,1,1,K,I),T3(1,1,1,1,L,J),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,10)
      CALL MATMUL(T3I(1,1,1,1,K,I),T3(1,1,1,1,J,L),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,22)
      CALL MATMUL(T3I(1,1,1,1,I,L),T3(1,1,1,1,K,J),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,8)
      CALL MATMUL(T3I(1,1,1,1,L,I),T3(1,1,1,1,K,J),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,23)
      CALL MATMUL(T3I(1,1,1,1,K,L),T3(1,1,1,1,J,I),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,8)
      CALL MATMUL(T3I(1,1,1,1,L,K),T3(1,1,1,1,J,I),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,6)
      CALL MATMUL(T3I(1,1,1,1,J,K),T3(1,1,1,1,L,I),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,10)
      CALL MATMUL(T3I(1,1,1,1,J,K),T3(1,1,1,1,I,L),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,1)
      CALL MATMUL(T3I(1,1,1,1,J,L),T3(1,1,1,1,I,K),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,8)
      CALL MATMUL(T3I(1,1,1,1,L,J),T3(1,1,1,1,I,K),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,21)
      call symt411(t4,nu,7)
      CALL TRT3ALL(NO,NU,T3,23)
      goto 9000
 1003 continue
      CALL TRT3ALL(NO,NU,T3,23)
      call vecmul(t4,nu4,half)
      CALL MATMUL(T3I(1,1,1,1,I,K),T3(1,1,1,1,J,L),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,10)
      CALL MATMUL(T3I(1,1,1,1,I,K),T3(1,1,1,1,L,J),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,8)
      CALL MATMUL(T3I(1,1,1,1,K,I),T3(1,1,1,1,L,J),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,10)
      CALL MATMUL(T3I(1,1,1,1,K,I),T3(1,1,1,1,J,L),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,12)
      CALL MATMUL(T3I(1,1,1,1,I,J),T3(1,1,1,1,K,L),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,21)
      CALL MATMUL(T3I(1,1,1,1,K,L),T3(1,1,1,1,J,I),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,6)
      CALL MATMUL(T3I(1,1,1,1,J,L),T3(1,1,1,1,K,I),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,10)
      CALL MATMUL(T3I(1,1,1,1,J,L),T3(1,1,1,1,I,K),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,14)
      CALL MATMUL(T3I(1,1,1,1,K,J),T3(1,1,1,1,I,L),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,10)
      CALL MATMUL(T3I(1,1,1,1,K,J),T3(1,1,1,1,L,I),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,23)
      CALL MATMUL(T3I(1,1,1,1,L,K),T3(1,1,1,1,I,J),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,18)
      CALL MATMUL(T3I(1,1,1,1,J,I),T3(1,1,1,1,K,L),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,12)
      call symt411(t4,nu,1)
      CALL TRT3ALL(NO,NU,T3,23)
 9000 continue
      if(i.eq.j)goto 801
      if(j.eq.k)goto 802
      if(k.eq.l)goto 803
      call symt3in(i,j,no,nu,ti,t3i)
      CALL MATMUL(TI,T3(1,1,1,1,K,L),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,7)
      call symt3in(i,k,no,nu,ti,t3i)
      CALL MATMUL(TI,T3(1,1,1,1,J,L),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,2)
      call symt3in(i,l,no,nu,ti,t3i)
      CALL MATMUL(TI,T3(1,1,1,1,K,J),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,9)
      call symt3in(k,j,no,nu,ti,t3i)
      CALL MATMUL(TI,T3(1,1,1,1,I,L),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,14)
      call symt3in(l,j,no,nu,ti,t3i)
      CALL MATMUL(TI,T3(1,1,1,1,K,I),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,22)
      call symt3in(k,l,no,nu,ti,t3i)
      CALL MATMUL(TI,T3(1,1,1,1,I,J),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,9)
      goto 800
 801  continue
      call symt3in(i,j,no,nu,ti,t3i)
      CALL MATMUL(TI,T3(1,1,1,1,K,L),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,7)
      call symt3in(i,k,no,nu,ti,t3i)
      call vecmul(t4,nu4,half)
      CALL MATMUL(TI,T3(1,1,1,1,J,L),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,2)
      call symt3in(i,l,no,nu,ti,t3i)
      CALL MATMUL(TI,T3(1,1,1,1,K,J),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,10)
      call symt411(t4,nu,6)
      CALL MTRANS(T4,NU,9)
      call symt3in(k,l,no,nu,ti,t3i)
      CALL MATMUL(TI,T3(1,1,1,1,I,J),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,9)
      goto 800
 802  continue
      call symt3in(i,j,no,nu,ti,t3i)
      call vecmul(t4,nu4,half)
      CALL MATMUL(TI,T3(1,1,1,1,K,L),T4,NU2,NU2,NOU,0,1)
c
      CALL MTRANS(T4,NU,14)
      call symt3in(l,j,no,nu,ti,t3i)
      CALL MATMUL(TI,T3(1,1,1,1,I,K),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,5)
      call symt411(t4,nu,7)
      CALL MTRANS(T4,NU,10)
      call symt3in(i,l,no,nu,ti,t3i)
      CALL MATMUL(TI,T3(1,1,1,1,K,J),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,9)
      call symt3in(k,j,no,nu,ti,t3i)
      CALL MATMUL(TI,T3(1,1,1,1,I,L),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,8)
      goto 800
 803  continue
      call symt3in(i,j,no,nu,ti,t3i)
      CALL MATMUL(TI,T3(1,1,1,1,K,L),T4,NU2,NU2,NOU,0,1)
      call vecmul(t4,nu4,half)
      CALL MTRANS(T4,NU,7)
      call symt3in(i,k,no,nu,ti,t3i)
      CALL MATMUL(TI,T3(1,1,1,1,J,L),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,13)
      call symt3in(k,j,no,nu,ti,t3i)
      CALL MATMUL(TI,T3(1,1,1,1,I,L),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,8)
      call symt411(t4,nu,1)
      CALL MTRANS(T4,NU,9)
      call symt3in(k,l,no,nu,ti,t3i)
      CALL MATMUL(TI,T3(1,1,1,1,I,J),T4,NU2,NU2,NOU,0,1)
      CALL MTRANS(T4,NU,9)
 800  continue
      call iexit(74)
      RETURN
      END
      subroutine symt3in(i,j,no,nu,ti,t3i)
      implicit double precision (a-h,o-z)
      integer a,b,c
      dimension ti(nu,nu,no,nu),t3i(nu,nu,no,nu,no,no)
      do 113 k1=1,no
      do 113 a=1,nu
      do 113 b=1,nu
      do 113 c=1,nu
      ti(a,b,k1,c)=t3i(a,b,k1,c,i,j)+t3i(b,a,k1,c,j,i)
 113  continue
      return
      end
      SUBROUTINE T4WT32MN(I,J,K,L,NO,NU,TI,T3,T3MN,T4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D,E,F
      DIMENSION TI(NU,NU,NU),T3(NO,NO,NU,NU,NU,NO),
     *T3MN(NU,NO,NO,NO,NO,NO),T4(NU,NU,NU,NU)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DATA ZERO/0.0D+0/,TWO/2.0D+0/
      call ienter(75)
      DO 10 I1=1,NO
      DO 10 J1=1,NO
      DO 10 K1=1,NO
      IF(I1.EQ.J1.AND.J1.EQ.K1)THEN
      CALL ZEROMA(TI,1,NU3)
      GO TO 9
      ENDIF
      CALL RDVT3ONW(I1,J1,K1,NU,TI)
 9    CONTINUE
      DO 11 A=1,NU
      DO 11 B=1,NU
      DO 11 C=1,NU
      T3(I1,J1,A,B,C,K1)=TI(A,B,C)
 11   CONTINUE
 10   CONTINUE
      CALL RDVT3IMN(NO,NU,TI,T3MN)
      CALL MATMUL(T3MN(1,1,1,I,J,K),T3(1,1,1,1,1,L),T4,NU,NU3,NO2,0,0)
      CALL MTRANS(T4,NU,1)
      CALL MATMUL(T3MN(1,1,1,I,J,L),T3(1,1,1,1,1,K),T4,NU,NU3,NO2,0,0)
      CALL MTRANS(T4,NU,10)
      CALL MATMUL(T3MN(1,1,1,I,K,L),T3(1,1,1,1,1,J),T4,NU,NU3,NO2,0,0)
      CALL MTRANS(T4,NU,12)
      CALL MATMUL(T3MN(1,1,1,L,I,K),T3(1,1,1,1,1,J),T4,NU,NU3,NO2,0,0)
      CALL MTRANS(T4,NU,11)
      CALL MATMUL(T3MN(1,1,1,J,I,K),T3(1,1,1,1,1,L),T4,NU,NU3,NO2,0,0)
      CALL MTRANS(T4,NU,1)
      CALL MATMUL(T3MN(1,1,1,J,I,L),T3(1,1,1,1,1,K),T4,NU,NU3,NO2,0,0)
      CALL MTRANS(T4,NU,11)
      CALL MATMUL(T3MN(1,1,1,K,I,L),T3(1,1,1,1,1,J),T4,NU,NU3,NO2,0,0)
      CALL MTRANS(T4,NU,1)
      CALL MATMUL(T3MN(1,1,1,K,I,J),T3(1,1,1,1,1,L),T4,NU,NU3,NO2,0,0)
      CALL MTRANS(T4,NU,11)
      CALL MATMUL(T3MN(1,1,1,L,I,J),T3(1,1,1,1,1,K),T4,NU,NU3,NO2,0,0)
      CALL MTRANS(T4,NU,3)
      CALL MATMUL(T3MN(1,1,1,L,J,K),T3(1,1,1,1,1,I),T4,NU,NU3,NO2,0,0)
      CALL MTRANS(T4,NU,12)
      CALL MATMUL(T3MN(1,1,1,K,L,J),T3(1,1,1,1,1,I),T4,NU,NU3,NO2,0,0)
      CALL MTRANS(T4,NU,12)
      CALL MATMUL(T3MN(1,1,1,J,K,L),T3(1,1,1,1,1,I),T4,NU,NU3,NO2,0,0)
      CALL MTRANS(T4,NU,20)
      call iexit(75)
      RETURN
      END
      subroutine drtrint3(no,nu,t)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical print
      common /flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension t(*)
      common/itrat/icycle,mx,icn
      common/totmem/mem
      print=iflags(1).gt.10
      i1=1            !vt3m
      i2=i1+nou3      !v
      i3=i2+nu4       !vm
      i4=i3+no3u      !ve
      i5=i4+nou3      !vt
      i6=i5+no4*nu2   !ti
      it=i6+nu3
      if(icycle.eq.1.and.print)write(6,99)mem,it
 99   format('Space usage in trint3: available - ',i8,'   used - ',i8)
      call trint3(no,nu,t(i1),t(i2),t(i3),t(i4),t(i5),t(i6))
      return
      end
      SUBROUTINE TRINT3(NO,NU,VT3M,V,VM,VE,VT,TI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VT3M(NU,NU,NU,NO),V(NU,NU,NU,NU),VM(NO,NU,NO,NO),
     *VE(NU,NU,NU,NO),TI(NU,NU,NU),VT(NU,NU,NO,NO,NO,NO)
      call ienter(46)
      no4u2=no*no*no*no*nu*nu
      CALL TRINT3ME(NO,NU,VT3M,V,VE,TI)
      CALL TRINT3MN(NO,NU,VE,V,VM,TI)
      CALL TRINVTME(NO,NU,VT,VM,TI,VE)
      CALL TRINVTEF(NO,NU,VT,VE,TI)
      call iexit(46)
      RETURN
      END
      SUBROUTINE TRINT3ME(NO,NU,TABCM,V,VE,VTM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SDTQ4
      DIMENSION TABCM(NU,NU,NU,NO),V(NU,NU,NU,NU),VE(NU,NO,NU,NU),
     *VTM(NU,NU,NU)
      COMMON/NEWOPT/NOPT(6)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DATA ZERO/0.0D+0/,HALF/0.5D+00/,TWO/2.0D+00/
      call ienter(47)
      SDTQ4=NOPT(1).GE.6
      I10=0
      IF(SDTQ4)I10=10
      CALL RDVEM2(I10,NO,NU,VTM,VE)
      DO 10 I=1,NO
      DO 10 J=1,i
      CALL ZEROMA(TABCM,1,NOU3)
      DO 9 K=1,NO
      IF(I.EQ.J.AND.J.EQ.K)GOTO 9
      CALL RDVT3ONW(I,J,K,NU,VTM)
      CALL VECCOP(NU3,TABCM(1,1,1,K),VTM)
 9    CONTINUE
      CALL VECMUL(VE,NOU3,TWO)
      CALL MATMUL(TABCM,VE,V,NU2,NU2,NOU,1,0)
      CALL VECMUL(VE,NOU3,HALF)
      CALL TRANMD(TABCM,NU,NU,NU,NO,23)
      CALL MATMUL(TABCM,VE,V,NU2,NU2,NOU,0,1)
      CALL TRANMD(TABCM,NU,NU,NU,NO,312)
      CALL MATMUL(TABCM,VE,V,NU2,NU2,NOU,0,1)
      CALL TRANMD(TABCM,NU,NU,NU,NO,13)
      CALL TRANMD(VE,NU,NO,NU,NU,14)
      CALL MATMUL(TABCM,VE,V,NU2,NU2,NOU,0,1)
      CALL MTRANS(V,NU,7)
      CALL TRANMD(TABCM,NU,NU,NU,NO,23)
      CALL MATMUL(TABCM,VE,V,NU2,NU2,NOU,0,1)
      CALL TRANMD(TABCM,NU,NU,NU,NO,13)
      CALL MTRANS(V,NU,8)
      CALL MATMUL(TABCM,VE,V,NU2,NU2,NOU,0,1)
      CALL MTRANS(V,NU,7)
      CALL MTRANS(V,NU,6)
      CALL TRANMD(VE,NU,NO,NU,NU,14)
      kk=it2(i,j)
      IASV=NO3+(kk-1)*NU+1
      CALL WRVA(IASV,NU,V)
 10   CONTINUE
      call iexit(47)
      RETURN
      END
      SUBROUTINE TRINT3MN(NO,NU,TABCM,V,VM,VTM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SDTQ4
      COMMON/NEWOPT/NOPT(6)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION TABCM(NU,NU,NU,NO),V(NU,NU,NU,NU),VM(NO,NU,NO,NO),
     *VTM(NU,NU,NU)
      call ienter(48)
      DATA ZERO/0.0D+0/,HALF/0.5D+00/,TWO/2.0D+00/
      SDTQ4=NOPT(1).GE.6
      I11=1
      IF (SDTQ4)I11=11
      CALL RDVEM2(I11,NU,NO,VTM,VM)
      IF(I11.EQ.11)CALL TRANMD(VM,NO,NU,NO,NO,34)
      CALL TRANMD(VM,NO,NU,NO,NO,13)
c
      DO 10 I=1,NO
      DO 10 J=1,i
      kk=it2(i,j)
      IASV=NO3+(kk-1)*NU+1
      CALL RDVA(IASV,NU,V)
c
      DO 15 N=1,NO
      CALL ZEROMA(TABCM,1,NOU3)
      DO 9 M=1,NO
      IF(I.EQ.M.AND.N.EQ.M)GOTO 9
      CALL RDVT3ONW(I,M,N,NU,VTM)
      CALL VECCOP(NU3,TABCM(1,1,1,M),VTM)
 9    CONTINUE
      CALL MATMUL(TABCM,VM(1,1,N,J),V,NU3,NU,NO,0,0)
 15   CONTINUE
      CALL MTRANS(V,NU,6)
      DO 151 N=1,NO
      CALL ZEROMA(TABCM,1,NOU3)
      DO 91 M=1,NO
      IF(J.EQ.M.AND.N.EQ.M)GOTO 91
      CALL RDVT3ONW(J,M,N,NU,VTM)
      CALL VECCOP(NU3,TABCM(1,1,1,M),VTM)
 91   CONTINUE
      CALL MATMUL(TABCM,VM(1,1,N,I),V,NU3,NU,NO,0,0)
 151  CONTINUE
c
      CALL MTRANS(V,NU,6)
      kk=it2(i,j)
      IASV=NO3+(kk-1)*NU+1
      CALL WRVA(IASV,NU,V)
 10   CONTINUE
      call iexit(48)
      RETURN
      END
      SUBROUTINE TRINVTEF(NO,NU,VT,VE,VTM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SDTQ4
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON/NEWOPT/NOPT(6)
      COMMON/NVTAP/NTITER,NQ2,NVT
      DIMENSION VE(NU,NU,NU,NO),VTM(NU,NU,NU),VT(NU,NU,NO,NO)
      call ienter(49)
      SDTQ4=NOPT(1).GE.6
      I10=0
      IF (SDTQ4)I10=10
      CALL RDVEM4(I10,NO,NU,VTM,VE)
      CALL TRANMD(VE,NU,NU,NU,NO,312)
      DO 500 I=1,NO
      DO 500 J=1,NO
      DO 500 K=1,NO
      kk=no2*(i-1)+no*(j-1)+k
      call rdgen(nvt,kk,nou2,vt(1,1,1,k))
      IF(I.EQ.J.AND.J.EQ.K) GOTO 500
      CALL RDVT3ONW(I,J,K,NU,VTM)
      CALL MATMUL(VTM,VE,VT(1,1,1,K),NU,NOU,NU2,0,0)
      CALL TRANMD(VTM,NU,NU,NU,1,12)
      kk=no2*(i-1)+no*(j-1)+k
      call wrgen(nvt,kk,nou2,vt(1,1,1,k))
 500  CONTINUE
      call iexit(49)
      RETURN
      END
      SUBROUTINE TRINVTME(NO,NU,VT,VM,VTM,TABCM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SDTQ4
      COMMON/NVTAP/NTITER,NQ2,NVT
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON/NEWOPT/NOPT(6)
      DIMENSION VM(NU,NO,NO,NO),
     *VTM(NU,NU,NU),VT(NU,NU,NO,NO),TABCM(NU,NU,NU,NO)
      DATA ZERO/0.0D+0/,HALF/0.5D+00/,TWO/2.0D+00/
      call ienter(50)
      SDTQ4=NOPT(1).GE.6
      NO4U2=NO3U*NOU
      CALL ZEROMA(VT,1,NO2U2)
      I11=1

      IF (SDTQ4)I11=11
      CALL RDVEM1(I11,NU,NO,VTM,VM)
      IF(I11.EQ.11) CALL TRANMD(VM,NU,NO,NO,NO,34)
      DO 5001 I=1,NO
      DO 5001 J=1,NO
      DO 101 K=1,NO
      IF(I.EQ.J.AND.J.EQ.K)THEN
      CALL ZEROMA(VTM,1,NU3)
      GO TO 91
      ENDIF
      CALL RDVT3ONW(I,J,K,NU,VTM)
 91   CONTINUE
      CALL VECCOP(NU3,TABCM(1,1,1,K),VTM)
 101  CONTINUE
      DO 5002 K=1,NO
      CALL MATMUL(TABCM,VM(1,1,1,K),VT(1,1,1,k),NU2,NO,NOU,1,0)
      CALL TRANMD(TABCM,NU,NU,NU,NO,23)
      CALL MATMUL(TABCM,VM(1,1,1,K),VT(1,1,1,k),NU2,NO,NOU,0,1)
      CALL TRANMD(TABCM,NU,NU,NU,NO,23)
      CALL TRANMD(VM(1,1,1,K),NU,NO,NO,1,23)
      CALL VECMUL(VM(1,1,1,K),NO2U,HALF)
      CALL MATMUL(TABCM,VM(1,1,1,K),VT(1,1,1,k),NU2,NO,NOU,0,1)
      CALL VECMUL(VM(1,1,1,k),NO2U,two)
      CALL TRANMD(VM(1,1,1,K),NU,NO,NO,1,23)
 5002 CONTINUE
      CALL TRANMD(VM(1,1,1,j),NU,NO,NO,1,23)
      DO 503 K=1,NO
      DO 1012 M=1,NO
      IF(I.EQ.K.AND.K.EQ.M)THEN
      CALL ZEROMA(VTM,1,NU3)
      GO TO 9110
      ENDIF
      CALL RDVT3ONW(I,K,M,NU,VTM)
 9110 CONTINUE
      CALL VECCOP(NU3,TABCM(1,1,1,M),VTM)
 1012 CONTINUE
      CALL TRANMD(TABCM,NU,NU,NU,NO,23)
      CALL MATMUL(TABCM,VM(1,1,1,J),VT(1,1,1,k),NU2,NO,NOU,0,1)
      CALL TRANMD(TABCM,NU,NU,NU,NO,23)
      kk=no2*(i-1)+no*(j-1)+k
      call wrgen(nvt,kk,nou2,vt(1,1,1,k))
 503  CONTINUE
      CALL TRANMD(VM(1,1,1,j),NU,NO,NO,1,23)
 5001 CONTINUE
      call iexit(50)
      RETURN
      END
      SUBROUTINE VT3EF(NO,NU,TI,VOE,T3,VT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,E,F
      DIMENSION T3(NU,NU,NU,NO,NO,NO),VOE(NU,NU,NO,NO),TI(NU,NU,NU),
     *     VT(NU,NO,NO,NO,NO,NO)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON/NEWTAPE/NT2T4,NT4INT
      COMMON/NEWT3/NT3INT,LT3INT
      DATA ZERO/0.0D+0/,TWO/2.0D+00/,HALF/0.5D+0/
      call ienter(83)
      NO3U3=NO3*NU3
      CALL RO2PPH(1,NO,NU,TI,VOE)
      DO 10 I=1,NO
      DO 10 J=1,NO
      DO 10 K=1,NO
      IF(I.EQ.J.AND.J.EQ.K)THEN
      CALL ZEROMA(TI,1,NU3)
      GO TO 9
      ENDIF
      CALL RDVT3ONW(I,J,K,NU,TI)
 9    CONTINUE
      DO 11 A=1,NU
      DO 11 B=1,NU
      DO 11 C=1,NU
      T3(A,B,C,I,J,K)=TI(A,B,C)
 11   CONTINUE
 10   CONTINUE
      DO 51 I=1,NO
      DO 51 J=1,NO
      DO 51 K=1,NO
      CALL MATMUL(T3(1,1,1,I,J,K),VOE,VT(1,1,1,I,J,K),NU,NO2,NU2,1,0)
 51   CONTINUE
      CALL WRVT3IMN(NO,NU,TI,VT)
      call iexit(83)
      RETURN
      END
      subroutine drvt3int(no,nu,t)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical print
      common /flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/itrat/icycle,mx,icn
      common/totmem/mem
      dimension t(*)
      print=iflags(1).gt.10
      no3u3=no3*nu3
      i1=1          !ti
      i2=i1+nu3     !voe
      i3=i2+no2u2   !t3
      i4=i3+no3u3   !vt
      it=i4+no4*nou
      if(icycle.eq.1.and.print)write(6,99)mem,it
 99   format('Space usage in vt3int: available - ',i8,'   used - ',i8)
      call vt3int(no,nu,t(i1),t(i2),t(i3),t(i4))
      return
      end
      SUBROUTINE VT3INT(NO,NU,TI,VOE,T3,VT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T3(NO,NO,NO,NU,NU,NU),VOE(NO,NU,NU,NO),TI(NU,NU,NU),
     *VT(NU,NO,NO,NO,NO,NO)
      call ienter(84)
      CALL VT3ME(NO,NU,TI,VOE,T3,VT)
      CALL VT3EF(NO,NU,TI,VOE,T3,VT)
      call iexit(84)
      RETURN
      END
      SUBROUTINE VT3ME(NO,NU,TI,VOE,T3,VT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,E,F
c      REAL CPTIME,CPTIM
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION T3(NU,NU,NU,NO,NO,NO),VOE(NU,NO,NU,NO),TI(NU,NU,NU),
     *VT(NU,NO,NO,NO,NO,NO)
      COMMON/RESLTS/CMP(30)
      COMMON/NEWTAPE/NT2T4,NT4INT
      COMMON/NEWT3/NT3INT,LT3INT
      DATA ZERO/0.0D+0/,TWO/2.0D+00/,HALF/0.5D+0/
      call ienter(85)
      NO3U3=NO3*NU3
      CALL RO2PHP(1,NO,NU,TI,VOE)
      DO 10 I=1,NO
      DO 10 J=1,NO
      DO 10 K=1,NO
      IF(I.EQ.J.AND.J.EQ.K)THEN
      CALL ZEROMA(TI,1,NU3)
      GO TO 9
      ENDIF
      CALL RDVT3ONW(I,J,K,NU,TI)
 9    CONTINUE
      DO 11 A=1,NU
      DO 11 B=1,NU
      DO 11 C=1,NU
      T3(A,B,C,K,I,J)=TI(A,B,C)
 11   CONTINUE
 10   CONTINUE
      CALL SYMT21(VOE,NU,NO,NU,NO,13)
      DO 51 I=1,NO
      DO 51 J=1,NO
      DO 51 K=1,NO
      CALL MATMUL(T3(1,1,1,1,I,J),VOE(1,1,1,K),TI,NU2,NU,NOU,1,0)
      CALL VECMUL(TI,NU3,HALF)
      IAS=NO2*(I-1)+NO*(J-1)+K
      CALL WRVT3I(IAS,NU,TI)
 51   CONTINUE
      CALL DESM21(VOE,NU,NO,NU,NO,13)
c
      CALL TRANMD(VOE,NU,NO,NU,NO,13)
      DO 151 I=1,NO
      DO 151 J=1,NO
      CALL TRANMD(T3(1,1,1,1,I,J),NU,NU,NU,NO,23)
      DO 152 K=1,NO
      CALL MATMUL(T3(1,1,1,1,I,J),VOE(1,1,1,K),TI,NU2,NU,NOU,1,1)
      CALL VECMUL(TI,NU3,HALF)
      IAS=NO3+NO2*(I-1)+NO*(J-1)+K
      CALL WRVT3I(IAS,NU,TI)
 152  CONTINUE
      CALL TRANMD(T3(1,1,1,1,I,J),NU,NU,NU,NO,23)
 151  CONTINUE
      CALL TRANMD(VOE,NU,NO,NU,NO,13)
c
      DO 510 I=1,NO
      DO 510 J=1,NO
      CALL TRANMD(T3(1,1,1,1,I,J),NU,NU,NU,NO,23)
      DO 511 K=1,NO
      CALL MATMUL(T3(1,1,1,1,I,J),VOE(1,1,1,K),TI,NU2,NU,NOU,1,1)
      CALL VECMUL(TI,NU3,HALF)
      IAS=NO2*(I-1)+NO*(J-1)+K
      CALL WRVT3ID(IAS,NU,TI,VT)
 511  CONTINUE
      CALL TRANMD(T3(1,1,1,1,I,J),NU,NU,NU,NO,312)
      DO 512 K=1,NO
      CALL MATMUL(T3(1,1,1,1,I,J),VOE(1,1,1,K),TI,NU2,NU,NOU,1,1)
      CALL VECMUL(TI,NU3,HALF)
      IAS=NO2*(I-1)+NO*(J-1)+K
      CALL WRVT3ID(IAS,NU,TI,VT)
 512  CONTINUE
 510  CONTINUE
      call iexit(85)
      RETURN
      END
      subroutine drvt4cd(no,nu,t)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical print
      common /flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/itrat/icycle,mx,icn
      common/totmem/mem
      dimension t(*)
      print=iflags(1).gt.10
      nou4=no*nu4
      i1=1          !ti
      i2=i1+nu3     !o4
      i3=i2+nou4    !v
      i4=i3+nu4     !t4
      i5=i4+nu4     !voe
      it=i5+no2u2
      if(icycle.eq.1.and.print)write(6,99)mem,it
 99   format('Space usage in vt4cd: available - ',i8,'   used - ',i8)
      call vt4cd(no,nu,t,t(i2),t(i3),t(i4),t(i5))
      return
      end
      SUBROUTINE VT4CD(NO,NU,TI,O4,V,T4,VOE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL IEJ,NE1,NE2
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION TI(NU,NU,NU),O4(NU,NU,NU,NU,NO),
     *V(NU,NU,NU,NU),VOE(NU,NO,NU,NO),T4(NU,NU,NU,NU)
      call ienter(86)
      CALL RO2PHP(1,NO,NU,TI,VOE)
      CALL TRANMD(VOE,NU,NO,NU,NO,13)
      nou4=no*nu4
      DO 1 I=1,NO
         DO 1 J=1,i
      kk=it2(i,j)
      IASV=NO3+(kk-1)*NU+1
      IEJ=I.EQ.J
      CALL RDVA(IASV,NU,V)
      DO 11 N=1,NO
      NE1=I.EQ.N
      NE2=J.EQ.N
      DO 10 M=1,NO
      IF (IEJ.AND.(NE2.OR.J.EQ.M).OR.M.EQ.N.AND.(NE1.OR.NE2))THEN
      CALL ZEROMA(O4(1,1,1,1,M),1,NU4)
      GO TO 10
      ENDIF
      CALL RDIJKLT3(I,J,M,N,NO,NU,T4,O4(1,1,1,1,M))
 10   CONTINUE
      CALL MATMUL(O4,VOE(1,1,1,N),V,NU3,NU,NOU,0,1)
 11   continue
      kk=it2(i,j)
      IASV=NO3+(kk-1)*NU+1
      CALL WRVA(IASV,NU,V)
 1    CONTINUE
      call iexit(86)
      RETURN
      END
      subroutine drvt4kl(no,nu,t)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical print
      common /flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/itrat/icycle,mx,icn
      common/totmem/mem
      dimension t(*)
      print=iflags(1).gt.10
      nou4=no*nu4
      i1=1          !ti
      i2=i1+nu3     !o4
      i3=i2+nou4    !t4
      i4=i3+nu4     !vt
      i5=i4+no2u2   !voe
      it=i5+no2u2
      if(icycle.eq.1.and.print)write(6,99)mem,it
 99   format('Space usage in vt4kl: available - ',i8,'   used - ',i8)
      call vt4kl(no,nu,t,t(i2),t(i3),t(i4),t(i5))
      return
      end
      SUBROUTINE VT4KL(NO,NU,TI,O4,T4,VT,VOE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL  IEJ,JEK,KEM
      COMMON/NVTAP/NTITER,NQ2,NVT
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION TI(NU,NU,NU),O4(NU,NU,NU,NU,NO),T4(NU,NU,NU,NU),
     *VT(NU,NU,NO,NO),VOE(NU,NU,NO,NO)
      data half/0.5d+0/
      call ienter(87)
      nou4=no*nu4
      no4u2=no4*nu2
      CALL RO2PPH(1,NO,NU,TI,VOE)
      call vecmul(voe,no2u2,half)
      CALL TRANMD(VOE,NU,NU,NO,NO,34)
      DO 11 I=1,NO
      DO 11 J=1,NO
      IEJ=I.EQ.J
      DO 11 K=1,NO
      JEK=J.EQ.K
      DO 110 M=1,NO
      KEM=K.EQ.M
      IF (IEJ.AND.(JEK.OR.J.EQ.M).OR.KEM.AND.(I.EQ.K.OR.JEK))THEN
      CALL ZEROMA(O4(1,1,1,1,M),1,NU4)
      GO TO 110
      ENDIF
      CALL RDIJKLT3(I,J,K,M,NO,NU,T4,O4(1,1,1,1,M))
 110  CONTINUE
      kk=no2*(i-1)+no*(j-1)+k
      call rdgen(nvt,kk,nou2,vt(1,1,1,k))
      CALL MATMUL(O4,VOE,VT(1,1,1,k),NU2,NO,NOU2,0,0)
      kk=no2*(i-1)+no*(j-1)+k
      call wrgen(nvt,kk,nou2,vt(1,1,1,k))
 11   CONTINUE
      call iexit(87)
      RETURN
      END
      subroutine drwdexc(no,nu,t)
      implicit double precision(a-h,o-z)
      logical print
      common /flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/itrat/icycle,mx,icn
      common/totmem/nmem
      dimension t(*)
      print=iflags(1).gt.10
      no2u2h=no2u2/2+1
      nou3h=nou3/2+1
      i1=1              !o1
      i2=i1+nou         !ti
      ni3=nu3
      if(nu3.lt.no2u2h)ni3=no2u2h
      i3=i2+ni3         !voe
      ni4=no2u2
      if(no2u2.lt.nou3h)ni4=nou3h
      i4=i3+ni4       !wdppph
      itot=i4+nou3
      if(icycle.eq.1.and.print)write(6,99)nmem,itot
 99   format('Space usage in wdexc: available - ',i8,'   used - ',i8)
      call wdexc(no,nu,t,t(i2),t(i3),t(i4))
      return
      end
      SUBROUTINE WDEXC(NO,NU,O1,TI,VOE,WDPPPH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/pak/intg,intr
      common/wpak/nfr(15),nsz(15)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION O1(1),TI(1),VOE(1),WDPPPH(1)
      call ienter(88)
      call ro1(nou,1,o1)
      CALL RO2HPP(1,NO,NU,TI,VOE)
      CALL RDVEM4(0,NO,NU,TI,WDPPPH)
      call trt1(no,nu,ti,o1)
      CALL MATMUL(O1,VOE,WDPPPH,NU,NOU2,NO,0,1)
      CALL WRVEM(10,NO,NU,TI,WDPPPH)
      call store(nou3,nnn,wdppph,voe)
      nsz(8)=nnn
      call vripak(intg,8,nnn,wdppph,voe)
C
      call trt1(nu,no,ti,o1)
      CALL RO2HPP(1,NO,NU,TI,VOE)
      CALL RDVEM4(1,NU,NO,TI,WDPPPH)
      CALL TRANSQ(VOE,NOU)
      CALL MATMUL(O1,VOE,WDPPPH,NO,NO2U,NU,0,0)
      CALL WRVEM(11,NU,NO,TI,WDPPPH)
      call store(no3u,nnn,wdppph,voe)
      nsz(9)=nnn
      call vripak(intg,9,nnn,wdppph,voe)
      call iexit(88)
      RETURN
      END
      subroutine drgetfock(no,nu,t,eh,ep)
      implicit double precision(a-h,o-z)
      logical print
      common /flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/itrat/icycle,mx,icn
      common/totmem/nmem
      dimension t(*),eh(no),ep(nu)
      print=iflags(1).gt.10
      n=no+nu
      i1=1
      i2=i1+nu2
      i3=i2+n*n
      i4=i3+n*n
      itot=i4+n*n
      if(icycle.eq.1.and.print)write(6,99)nmem,itot
 99   format('Space usage in getfock: available - ',i8,'   used - ',i8)
      call getfock(no,nu,n,t,t(i2),t(i3),t(i4),eh,ep)
      return
      end
      SUBROUTINE GETFOCK(NO,NU,N,FPP,FOCK,C,SCR,OEH,OEP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      COMMON/ORBINF/NBASIS,NCOMP
      DIMENSION FPP(NU,NU),FOCK(N,N),C(N,N),SCR(N,N),OEH(NO),OEP(NU)
      NLAST=5*NO+2*NU
      CALL ZEROMA(FPP,1,NU2)
      call wr(nall4,nlast+1,no2,fpp)
      call wr(nall4,nlast+2,nou,fpp)
      call wr(nall4,nlast+3,nu2,fpp)
      RETURN
      N2=N*N
      CALL GETREC(20,'JOBARC','FOCKA',iintfp*N2,FOCK)
      CALL GETREC(20,'JOBARC','SCFEVCA0',iintfp*N2,C)
      CALL CFCVEC(N,FOCK,C,SCR)
      call cpfock(no,no,0,0,fpp,fock,oeh)
      call wr(nall4,nlast+1,no2,fpp)
      call cpfock(no,nu,0,no,fpp,fock,oeh)
      call wr(nall4,nlast+2,nou,fpp)
      call cpfock(nu,nu,no,no,fpp,fock,oep)
      call wr(nall4,nlast+3,nu2,fpp)
      RETURN
      END
      subroutine cpfock(n1,n2,iad,jad,n,fock,fp,oe)
      implicit double precision(a-h,o-z)
      dimension fp(n1,n2),fock(n,n),oe(1)
      do 10 i=1,n1
      do 10 j=1,n2
         fp(i,j)=fock(i+iad,j+jad)
 10   continue
      if (iad.eq.jad)then
         do 1 i=1,n1
            oe(i)=fp(i,i)
 1       continue
      endif
      return
      end
      subroutine drtestd2(no,nu,t,eh,ep)
      implicit double precision(a-h,o-z)
      logical print
      common /flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/itrat/icycle,mx,icn
      common/totmem/nmem
      dimension t(*),eh(no),ep(nu)
      print=iflags(1).gt.10
      i1=1
      i2=i1+nou2
      i3=i2+no2u2
      itot=i3+no2u2
      if(icycle.eq.1.and.print)write(6,99)nmem,itot
 99   format('Space usage in testd2: available - ',i8,'   used - ',i8)
      call testd2(no,nu,t,t(i2),t(i3),eh,ep)
      return
      end
      SUBROUTINE TESTD2(NO,NU,ti,T2VO,VOE,OEH,OEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      COMMON/ITRAT/ITER,MAXIT,ICNV
      DIMENSION T2VO(NO,NU,NU,NO),VOE(NO,NU,NU,NO),ti(no,nu,nu),
     *oeh(no),oep(nu)
      DATA ZERO/0.0D+0/,TWO/2.0D+0/
      call ro2hpp(0,no,nu,ti,t2vo)
      call ro2hpp(1,no,nu,ti,voe)
      E2=ZERO
      E2T=ZERO
      EWR=ZERO
      ET2=ZERO
      E2TS=ZERO
      DO 21 I=1,NO
      DO 21 J=1,NO
      DO 21 A=1,NU
         DO 21 B=1,NU
      DENOM=OEH(I)+OEH(J)-OEP(A)-OEP(B)
      E2 =E2 + VOE(I,A,B,J)*(TWO*VOE(I,A,B,J)-VOE(I,B,A,J))/DENOM
      ET2=ET2+VOE(I,A,B,J)*(TWO*T2VO(I,A,B,J)-T2VO(I,B,A,J))/DENOM
      E2T=E2T+VOE(I,A,B,J)*(TWO*T2VO(I,A,B,J)-T2VO(I,B,A,J))
      EWR=EWR+T2VO(I,A,B,J)*(TWO*T2VO(I,A,B,J)-T2VO(I,B,A,J))*DENOM
   21 CONTINUE
      IF(ITER.EQ.1)WRITE(6,10)E2
 10   FORMAT('SECOND ORDER ENERGY:',F15.10)
 11   FORMAT('WIGNER RULE DOUBLES,ERG:',2F15.10)
      RETURN
      END
      SUBROUTINE RDVEM4Anew(IN,NO,NU,TI,VEM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C
      DIMENSION TI(NU,NU,NU),VEM(NU,NU,NU,NO)
      DATA ZERO/0.0D+0/
      call ienter(110)
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
      IF(IN.GE.10)I11=NO+NU
      INO=1
      IF(IN.EQ.1.OR.IN.EQ.11)INO=0
      NLAST=3*NNO+NNU+I11
      DO 1 I=1,NO
      IASV=NLAST+NU*INO+I
      READ(NALL4,REC=IASV)TI
      DO 1 A=1,NU
      DO 1 B=1,NU
      DO 1 C=1,NU
      VEM(A,B,C,I)=TI(A,B,C)
 1    CONTINUE
      call iexit(110)
      return
      END
      SUBROUTINE RO2HPHnew(INO,NO,NU,TI,T2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      COMMON /NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION TI(NU,NU,NO),T2(NO,NU,NO,NU)
      call ienter(111)
      NLAST=0
      NOU2=NO*NU*NU
      INN=INO
      IF(INO.GE.10) THEN
      INN=INO-10
      NLAST=6*NO+4*NU+7
      ENDIF
      DO 1 I=1,NO
      IASV=NLAST+I+INN*NO
      READ(NALL4,REC=IASV)TI
      DO 1 J=1,NO
      DO 1 A=1,NU
      DO 1 B=1,NU
      T2(I,A,J,B)=TI(A,B,J)
 1    CONTINUE
      NO2U2=NO*NO*NU*NU
      call iexit(111)
      RETURN
      END
      SUBROUTINE SYMT411(V,NU,ID)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      DIMENSION V(NU,NU,NU,NU)
      DATA TWO/2.0D+0/
      if (id.eq.1)goto 1
      if (id.eq.6)goto 6
      if (id.eq.7)goto 7
      if (id.eq.8)goto 8
      if (id.eq.8)goto 9
      if (id.eq.10)goto 10
      if (id.eq.11)goto 11
 1    CONTINUE
      DO 100 A=1,NU
      DO 100 B=1,NU
      DO 100 C=1,NU
      DO 100 D=1,C
      x=v(a,b,c,d)+v(a,b,d,c)
      V(A,B,C,D)=x
      V(A,B,D,C)=x
  100 CONTINUE
      GO TO 1000
 6    CONTINUE
      DO 101 A=1,NU
      DO 101 B=1,A
      DO 101 C=1,NU
      DO 101 D=1,NU
      x=v(a,b,c,d)+v(b,a,c,d)
      V(A,B,C,D)=x
      V(b,a,c,d)=x
 101  CONTINUE
      GO TO 1000
 7    CONTINUE
      DO 102 A=1,NU
      DO 102 B=1,NU
      DO 102 C=1,B
      DO 102 D=1,NU
      x=v(a,b,c,d)+v(a,c,b,d)
      V(A,B,C,D)=x
      V(a,c,b,d)=x
 102  CONTINUE
      GO TO 1000
 8    CONTINUE
      DO 103 A=1,NU
      DO 103 B=1,NU
      DO 103 C=1,A
      DO 103 D=1,NU
      x=v(a,b,c,d)+v(c,b,a,d)
      V(A,B,C,D)=x
      V(c,b,a,d)=x
 103  CONTINUE
      GO TO 1000
 9    CONTINUE
      DO 203 A=1,NU
      DO 203 B=1,NU
      DO 203 C=1,A
      limd=nu
      if(a.eq.c)limd=b
      DO 203 D=1,LIMD
      x=v(a,b,c,d)+v(c,d,a,b)
      V(A,B,C,D)=x
      V(c,d,a,b)=x
 203  CONTINUE
      GO TO 1000
 10   CONTINUE
      DO 104 A=1,NU
      DO 104 B=1,NU
      DO 104 C=1,NU
      DO 104 D=1,B
      x=v(a,b,c,d)+v(a,d,c,b)
      V(A,B,C,D)=x
      V(a,d,c,b)=x
 104  CONTINUE
      GO TO 1000
 11   CONTINUE
      DO 105 A=1,NU
      DO 105 B=1,NU
      DO 105 C=1,NU
      DO 105 D=1,A
      x=v(a,b,c,d)+v(d,b,c,a)
      V(A,B,C,D)=x
      V(d,b,c,a)=x
 105  CONTINUE
      GO TO 1000
 1000 CONTINUE
      END

