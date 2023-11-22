      subroutine drcct1(no,nu,t,eh,ep)
      implicit double precision (a-h,o-z)
      logical print
      common/flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/itrat/icycle,mx,icn
      common/totmem/nmem
      dimension t(*),eh(no),ep(nu)
      print=iflags(1).gt.10
      i1=1
      i2=i1+nou
      i3=i2+nou
      i4=i3+nu3
      i5=i4+no2u2
      itot=i5+nou3
      if(icycle.eq.1.and.print)write(6,99)nmem,itot
      call cct1(no,nu,t,t(i2),t(i3),t(i4),t(i5),eh,ep)
 99   format('Space usage in cct1: available - ',i8,'   used - ',i8)
      return
      end
      SUBROUTINE CCT1(NO,NU,O1,T1,TI,T2,VE,OEH,OEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SDT1,SDT
      COMMON/NEWOPT/NOPT(6)
      common/unpak/nt3,lnt3
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON /ITRAT/ITER,MAXIT,ICCNV
      DIMENSION O1(1),T1(1),TI(1),T2(1),VE(1),OEH(NO),OEP(NU)
      DATA TWO/2.0D+00/,HALF/0.5D+00/
      call ienter(2)
      SDT1=NOPT(1).GT.0
      SDT=NOPT(1).GT.1
      NLAST=5*NO+2*NU
c      call getl2(no,nu,ti,t2)
c      call chksum('getl2   ',t2,no2u2)
c      stop'done'
      call ro1(nou,1,o1)
      call rf(nou,nlast+4,t2)
      call rf(no2,nlast+5,ve)
      CALL MATMULT(t2,O1,ve,NO,NO,NU,0,1)
c 7
      CALL TMATMUL(ve,O1,    T1,No,Nu,NO,1,1)
c 8
      call rf(nu2,nlast+6,ve)
      CALL MATMUL(O1,  ve,   T1,NO,NU,NU,  0,0)
c 9
      CALL RO2HPH(0,NO,NU,TI,ve)
      CALL SYMT21(ve,NO,NU,NO,NU,13)
      CALL MATMULT(ve,    t2,T1,NOU,1,NOU,0,0)
c 10
      CALL RO2HPH(1,NO,NU,TI,VE)
      CALL VECMUL(O1,NOU,TWO)
      CALL MATMUL(VE,   O1,  T1,NOU,1,NOU,0,0)
      CALL VECMUL(O1,NOU,HALF)
c 11
      CALL RO2HPH(2,NO,NU,TI,VE)
      CALL MATMUL(VE,O1,     T1,NOU,1,NOU,0,1)
c 12
      CALL RO2HPP(0,NO,NU,TI,T2)
      CALL INSITU(NO,NU,NU,NO,TI,T2,13)
      CALL RDVEM4(0,NO,NU,TI,VE)
      CALL SYMT21(VE,NU,NU,NU,NO,23)
      CALL TMATMULT(T2,VE,T1,No,Nu,NOU2,0,0)
c 13
      CALL RDVEM1(1,NU,NO,TI,VE)
      CALL SYMT21(VE,NU,NO,NO,NO,23)
      CALL TMATMULT(ve,T2,T1,No,Nu,NO2U,0,1)
      IF (ITER.GT.1.AND.SDT1)then
        if(sdt)then
           CALL T1WT3(NO,NU,TI,T1,t2,VE)
        else
           call rd(nt3,1,nou,ti)
           call vecadd(t1,ti,nou)
        endif
      endif
      call rf(nou,nlast+2,ve)
      CALL VECADD(T1,ve,NOU)
      CALL DENMT1(T1,O1,OEH,OEP,NO,NU)
c      call zeroma(t1,1,nou)
      CALL WO1(NOU,2,T1)
      call iexit(2)
      RETURN
      END
      subroutine RF(n,irec,a)
      implicit double precision(a-h,o-z)
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      dimension a(n)
      READ(NALL4,REC=IREC)A
      RETURN
      END
      subroutine RO1(n,irec,a)
      implicit double precision(a-h,o-z)
      COMMON/NEWCCSD/NTT2
      dimension a(n)
      READ(NTT2,REC=IREC)A
      RETURN
      END
      subroutine WO1(n,irec,a)
      implicit double precision(a-h,o-z)
      COMMON/NEWCCSD/NTT2
      dimension a(n)
      WRITE(NTT2,REC=IREC)A
      RETURN
      END
      subroutine drcct2(no,nu,t,eh,ep)
      implicit double precision (a-h,o-z)
      logical print
      common/flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/wpak/nfr(30),nsz(30)
      common/totmem/nmem
      common/itrat/icycle,mx,icn
      dimension t(1),eh(no),ep(nu)
      print=iflags(1).gt.10
      i1=1                !t2
      i2=i1+no2u2         !ti
      i3=i2+nu3           !o2
      i4=i3+no2u2         !io2
      i5=i4+no2u2/2+1     !v
c      ni5=nsz(1)
      ni5=nu4
      if(nou3.gt.ni5)ni5=nou3
      i6=i5+ni5           !iv
      ni6=nsz(1)/2+1
      nni6=nou3/2+1
      if(nni6.gt.ni6)ni6=nni6
      i7=i6+ni6           !nda
      i8=i7+nu2/2+1       !ndb
      itot=i8+ nu2/2+1     
      if(icycle.eq.1.and.print)write(6,99)nmem,itot
 99   format('Space usage in cct2: available - ',i8,'   used - ',i8)
      call cct2(no,nu,t,t(i2),t(i3),t(i4),t(i5),t(i6),t(i7),t(i8),
     &eh,ep)
      return
      end
      SUBROUTINE CCT2(NO,NU,T2,TI,O2,io2,V,iv,nda,ndb,oeh,oep)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SDT1,SDT
      COMMON/ENERGY/ ECORR(500,2),IXTRLE(500)
      common/pak/intg,intr
      common/lnt/lrecb,lrecw
      common/wpak/nfr(30),nsz(30)
      COMMON/ENRSDT/ECCSDT
      COMMON/ITRAT/ITER,MAXIT,ICCNV
      COMMON/NEWOPT/NOPT(6)
      common/newccsd/ntt2
      common/unpak/nt3,lnt3
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION TI(1),O2(no,nu,nu,no),T2(1),V(1),OEH(NO),OEP(NU),
     &iv(1),io2(1),nda(1),ndb(1)
      DATA TWO/2.0D+00/,HALF/0.5D+00/
      call ienter(3)
      SDT1=NOPT(1).GT.0
      SDT =NOPT(1).GT.1
      NLAST=5*NO+2*NU
      call zeroma(t2,1,no2u2)
      CALL RO2HHP(0,NO,NU,TI,O2)
c 31
      call ienter(111)
      CALL RF(NU2,NLAST+6,V)
      CALL MATMUL(O2,V,T2,NO2U,NU,NU,1,0)
      CALL IEXIT(111)
c 32
      call ienter(112)
      CALL RO1(NOU,1,v)
      call adt12(1,no,nu,v,o2)
c      call reapak(intg,1,v,iv)
      call rdov4(0,no,nu,ti,v)
      call mtrans(v,nu,7)
c      call scattr(nu4,nsz(1),v,iv)
      call sym2v(nu,v)
      call sym2o2(no,nu,o2)
      nx=nu*(nu+1)/2
      call matmul(o2,v,t2,no2,nu2,nx,0,0)
      CALL IEXIT(112)
      CALL RO2HHP(0,NO,NU,TI,O2)
      CALL RO1(NOU,1,v)
      call adt12(1,no,nu,v,o2)
      CALL VECMUL(O2,NO2U2,HALF)
c 33
      call ienter(113)
      CALL RDOV4(11,NU,NO,TI,V)
      call mtrans(V,no,9)
      CALL MATMUL(V,O2 ,T2,NO2,NU2,NO2,0,0)
      CALL VECMUL(O2,NO2U2,TWO)
      CALL RO1(NOU,1,v)
      call adt12(2,no,nu,v,o2)
      CALL IEXIT(113)
c 34
C get fhh
      CALL TRANMD(T2,NO,NO,NU,NU,34)
      CALL INSITU(NO,NO,NU,NU,TI,T2,13)
      CALL TRANSQ(T2,NOU)
      CALL RO2HPP(0,NO,NU,TI,O2)
      call ienter(114)
      CALL RF(NO2,NLAST+5,V)
      CALL MATMUL(O2,V,T2,NOU2,NO,NO,0,1)
      CALL IEXIT(114)
      CALL RO2HHP(0,NO,NU,TI,O2)
      CALL INSITU(NO,NO,NU,NU,TI,O2,23)
c 35
      call ienter(115)
      CALL SYMT21(O2,NO,NU,NO,NU,13)
      call reapak(intg,10,v,iv)
      call scattr(no2u2,nsz(10),v,iv)
      call matmul(o2,v,t2,nou,nou,nou,0,0)
      CALL DESM21(O2,NO,NU,NO,NU,13)
      CALL IEXIT(115)
      call ro2hph(0,no,nu,ti,o2)
c 36
      call ienter(116)
      call reapak(intg,11,v,iv)
      call scattr(no2u2,nsz(11),v,iv)
      call matmul(o2,v,t2,nou,nou,nou,0,1)
      CALL IEXIT(116)
c 37
      call ienter(117)
      CALL TRANMD(O2,NO,NU,No,Nu,13)
      CALL TRANMD(T2,NO,NU,NU,NO,23)
      call matmul(o2,v,t2,nou,nou,nou,0,1)
      CALL TRANMD(T2,NO,NU,NU,NO,23)
      CALL IEXIT(117)
c 38
      call ienter(118)
      CALL RO1(NOU,1,O2)
      call reapak(intg,12,v,iv)
      call scattr(nou3,nsz(12),v,iv)
      call matmul(o2,v,t2,no,nou2,nu,0,0)
      CALL IEXIT(118)
c 39
      call ienter(119)
      CALL TRANSQ(T2,NOU)
      call trt1(no,nu,ti,o2)
      call reapak(intg,13,v,iv)
      call scattr(no3u,nsz(13),v,iv)
      CALL MATMUL(O2,V,T2,NU,NO2U,NO,0,1)
      CALL TRANSQ(T2,NOU)
      CALL IEXIT(119)
      IF (ITER.GT.1.AND.SDT1) THEN
         if(sdt)then
            CALL INSITU(NO,NU,NU,NO,TI,T2,13)
            call drt2wt3(no,nu,t2,ti)
         else
            do 109 i=1,no
               call rd(nt3,i+1,nou2,o2(1,1,1,i))
 109        continue
            call vecadd(t2,o2,no2u2)
         endif
      ENDIF
      CALL SYMETR(T2,NO,NU)

      IF (ITER.GT.1.AND.SDT) THEN
         CALL INSITU(NO,NU,NU,NO,TI,T2,13)
         CALL TRANMD(T2,NU,NU,NO,NO,34)
         CALL T2FT3(NO,NU,T2,TI,V)
         CALL TRANMD(T2,NU,NU,NO,NO,34)
         CALL INSITU(Nu,NU,NO,NO,TI,T2,13)
      ENDIF
      CALL RO2HPP(1,NO,NU,TI,O2)
      CALL VECADD(T2,O2,NO2U2)
      CALL T2DEN(T2,TI,OEH,OEP,NO,NU)
      CALL WO2HPP(NO,NU,TI,T2)
c 
      CALL SYMT21(O2,NO,NU,NU,NO,23)
      et2=ddot(no2u2,t2,1,o2,1)
      CALL RO1(NOU,2,T2)
      CALL RF(NOU,NLAST+2,V)
      CALL VECMUL(V,NOU,TWO)
      CALL MATMUL(T2,O2,V,1,NOU,NOU,0,0)
      CALL TRT1(NO,NU,TI,T2)
      et1=ddot(nou,t2,1,v,1)
      ECCSDT=ET1+ET2
      ecorr(iter+1,1)=et2+et1
      ecorr(iter+1,2)=et2
      call iexit(3)
      RETURN
      END
      SUBROUTINE INITT1(NO,NU,T1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/NEWCCSD/NTT2
      INTEGER A
      DIMENSION T1(NO,NU)
      DATA ZERO/0.0D+0/
      DO 10 I=1,NO
      DO 10 A=1,NU
      T1(I,A)=ZERO
 10   CONTINUE
      WRITE(NTT2,REC=1)T1
      RETURN
      END
      SUBROUTINE INSI12(N1,N2,N3,TI,O2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TI(N1,N2,N3),O2(N2,N1,N3)
      N=N1*N2*N3
      CALL VECCOP(N,TI,O2)
      DO 10 I=1,N1
      DO 10 J=1,N2
      DO 10 K=1,N3
      O2(J,I,K)=TI(I,J,K)
 10   CONTINUE
      RETURN
      END
      SUBROUTINE INSI13(N1,N2,N3,TI,O2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TI(N1,N2,N3),O2(N3,N2,N1)
      N=N1*N2*N3
      CALL VECCOP(N,TI,O2)
      DO 10 I=1,N1
      DO 10 J=1,N2
      DO 10 K=1,N3
      O2(K,J,I)=TI(I,J,K)
 10   CONTINUE
      RETURN
      END
      SUBROUTINE INSI23(N1,N2,N3,TI,O2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TI(N1,N2,N3),O2(N1,N3,N2)
      N=N1*N2*N3
      CALL VECCOP(N,TI,O2)
      DO 10 I=1,N1
      DO 10 J=1,N2
      DO 10 K=1,N3
      O2(I,K,J)=TI(I,J,K)
 10   CONTINUE
      RETURN
      END
      SUBROUTINE INSITU(N1,N2,N3,N4,TI,O2,IC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TI(N1,N2,N3),O2(N1,N2,N3,N4)
      call ienter(8)
      DO 10 I=1,N4
      IF (IC.EQ.12)GOTO 12
      IF (IC.EQ.23)GOTO 23
      GOTO 13
 12   CALL INSI12(N1,N2,N3,TI,O2(1,1,1,I))
      GOTO 10
 23   CALL INSI23(N1,N2,N3,TI,O2(1,1,1,I))
      GOTO 10
 13   CALL INSI13(N1,N2,N3,TI,O2(1,1,1,I))
 10   CONTINUE
      call iexit(8)
      RETURN
      END
      subroutine drintrit2(no,nu,t)
      implicit double precision (a-h,o-z)
      logical print
      common/flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/totmem/nmem
      common/itrat/icycle,mx,icn
      dimension t(*)
      print=iflags(1).gt.10
      nou3h=nou3/2+1
      i1=1                !o1
      i2=i1+nou           !vm
      i3=i2+no3u          !voe
      i4=i3+no2u2         !ve
      i5=i4+nou3          !ti
      ni5=nu3
      if(nu3.lt.nou3h)ni5=nou3h
      itot=i5+ni5
      if(icycle.eq.1.and.print)write(6,99)nmem,itot
      call intrit2(no,nu,t,t(i2),t(i3),t(i4),t(i5))
 99   format('Space usage in intrit2: available - ',i8,'   used - ',i8)
      return
      end
      SUBROUTINE INTRIT2(NO,NU,O1,vm,VOE,VE,ti)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/pak/intg,intr
      common/wpak/nfr(30),nsz(30)
      DIMENSION O1(1),TI(1),VOE(1),VM(1),VE(1)
      DATA ZERO/0.0D+0/,TWO/2.0D+00/,HALF/0.5D+00/
      call ienter(15)
      CALL RO1(NOU,1,O1)
      CALL RDVEM4(0,NO,NU,TI,VE)
      CALL RO2HPP(1,NO,NU,TI,VOE)
C 20
      CALL TRANMD(VE,NU,NU,NU,NO,12)
      call matmul(o1,ve,voe,no,nou2,nu,0,0)
      CALL TRANMD(VE,NU,NU,NU,NO,12)
c 21 5 diagramow
      CALL TRT1(NO,NU,TI,O1)
      CALL TRANMD(VOE,NO,NU,NU,NO,1234)
      call matmul(o1,voe,ve,nu,nou2,no,0,1)
      CALL TRANMD(VE,NU,NU,NU,NO,231)
C 22
      CALL RO2HPP(2,NO,NU,TI,VOE)
      call matmul(o1,voe,ve,nu,nou2,no,0,1)      
      CALL TRANMD(VE,NU,NU,NU,NO,13)
      call store(nou3,nnn,ve,ti)
      nsz(12)=nnn
      call vripak(intg,12,nnn,ve,ti)
      call scattr(nou3,nnn,ve,ti)
C 27
      call rdvem3(0,no,nu,ti,ve)
      call tranmd(ve,nu,nu,no,nu,14)
      call ro2hhp(0,no,nu,ti,voe)
      CALL RDVEM4(1,NU,NO,TI,VM)
      call tranmd(vm,no,no,no,nu,312)
      call matmul(voe,ve,vm,no2,nou,nu2,0,0)
      call tranmd(vm,no,no,no,nu,231)
      call store(no3u,nnn,vm,ti)
      nsz(13)=nnn
      call vripak(intg,13,nnn,vm,ti)
      call iexit(15)
      RETURN
      END

      SUBROUTINE RDOVT2(INO,NO,NU,TI,V)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER A,B,C
      COMMON/NEWCCSD/NTT2
      DIMENSION TI(NU,NU,NU),V(NU,NU,NU,NU)
      DATA ZERO/0.0D+0/
      NNO=NO
      NNU=NU
      IF (INO.EQ.1) THEN
      NNO=NU
      NNU=NO
      ENDIF
      NLAST=4*NNO+NNU+2
      DO 1 I=1,NU
      IASV=NLAST+NNU*INO+I
      READ(NTT2,REC=IASV)TI
      DO 1 A=1,NU
      DO 1 B=1,NU
      DO 1 C=1,NU
      V(A,B,C,I)=TI(A,B,C)
 1    CONTINUE
      RETURN
      END
      SUBROUTINE RDT1HP(NO,NU,T1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/NEWCCSD/NTT2
      DIMENSION T1(NO,NU)
      READ(NTT2,REC=1)T1
      RETURN
      END
      SUBROUTINE RDT1PH(NO,NU,T1,T1PH)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER A
      COMMON/NEWCCSD/NTT2
      DIMENSION T1(NO,NU),T1PH(NU,NO)
      READ(NTT2,REC=1)T1
      DO 1 I=1,NO
      DO 1 A=1,NU
      T1PH(A,I)=T1(I,A)
 1    CONTINUE
      RETURN
      END



      SUBROUTINE RDVEMT2(INO,NO,NU,TI,VEM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C
      COMMON/NEWCCSD/NTT2
      DIMENSION TI(NU,NU,NU),VEM(NU,NU,NU,NO)
      DATA ZERO/0.0D+0/
      NNO=NO
      NNU=NU
      IF (INO.EQ.1) THEN
      NNO=NU
      NNU=NO
      ENDIF
      NLAST=3*NNO+2
      DO 1 I=1,NO
      IASV=NLAST+NU*INO+I
      READ(NTT2,REC=IASV)TI
      DO 1 A=1,NU
      DO 1 B=1,NU
      DO 1 C=1,NU
      VEM(A,B,C,I)=TI(A,B,C)
 1    CONTINUE
      NOU3=NO*NU*NU*NU
      RETURN
      END
      SUBROUTINE RDVOET2(NO,NU,TI,VOE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER A,B
      COMMON/NEWCCSD/NTT2
      DIMENSION TI(NU,NU,NO),VOE(NO,NU,NO,NU)
      NLAST=5*NO+2*NU+2
      DO 2 I=1,NO
      IASV=NLAST+I
      READ(NTT2,REC=IASV)TI
      DO 1 A=1,NU
      DO 1 B=1,NU
      DO 1 J=1,NO
      VOE(I,A,J,B)=TI(A,B,J)
 1    CONTINUE
 2    CONTINUE
      RETURN
      END
      SUBROUTINE RDINT2(NTT2,NLAST,NO,NOU2,VOE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION VOE(NOU2,NO)
      DO 2 I=1,NO
         call rd(ntt2,nlast+i,nou2,voe(1,i))
 2    CONTINUE
      RETURN
      END
      SUBROUTINE RDVOT2(NO,NU,TI,VOE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER A,B
      COMMON/NEWCCSD/NTT2
      DIMENSION TI(NU,NU,NO),VOE(NO,NU,NO,NU)
      NLAST=6*NO+2*NU+2
      DO 2 I=1,NO
      IASV=NLAST+I
      READ(NTT2,REC=IASV)TI
      DO 1 A=1,NU
      DO 1 B=1,NU
      DO 1 J=1,NO
      VOE(I,A,J,B)=TI(A,B,J)
 1    CONTINUE
 2    CONTINUE
      RETURN
      END
      SUBROUTINE T1WT3ijk(i,j,k,NO,NU,T1,VOE,ti,t3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION t1(nu,no),VOE(Nu,NU,No,NO),t3(1),ti(1)
      call ienter(129)
      call smt3four(NU,T3,ti)
      call matmul(t3,voe(1,1,j,k),t1(1,i),nu,1,nu2,0,0)
      call iexit(129)
      return
      end
      SUBROUTINE drt1wt3ij(i,j,NO,NU,T1,VOE,ti,t3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical sdt
      COMMON/NEWOPT/NOPT(6)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION t1(1),VOE(1),t3(1),ti(1)
      sdt=nopt(1).gt.1
      if(sdt)CALL RDVT3ONW(i,j,j,NU,TI)
      call T1WT3ijk(i,j,j,NO,NU,T1,VOE,ti,t3)
      call trant3(ti,nu,2)
      call T1WT3ijk(j,i,j,NO,NU,T1,VOE,ti,t3)
      call trant3(ti,nu,1)
      call T1WT3ijk(j,j,i,NO,NU,T1,VOE,ti,t3)
      return
      end
      SUBROUTINE drt1wt3jk(j,k,NO,NU,T1,VOE,ti,t3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical sdt
      common /newopt/nopt(6)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION t1(1),VOE(1),t3(1),ti(1)
      sdt=nopt(1).gt.1
      if(sdt)CALL RDVT3ONW(j,j,k,NU,TI)
      call T1WT3ijk(j,j,k,NO,NU,T1,VOE,ti,t3)
      call trant3(ti,nu,1)
      call T1WT3ijk(j,k,j,NO,NU,T1,VOE,ti,t3)
      call trant3(ti,nu,2)
      call T1WT3ijk(k,j,j,NO,NU,T1,VOE,ti,t3)
      return
      end
      SUBROUTINE drt1wt3ijk(i,j,k,NO,NU,T1,VOE,ti,t3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical sdt
      common /newopt/nopt(6)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION t1(1),VOE(1),t3(1),ti(1)
      sdt=nopt(1).gt.1
      if(sdt)CALL RDVT3ONW(i,j,k,NU,TI)
      call T1WT3ijk(i,j,k,NO,NU,T1,VOE,ti,t3)
      call trant3(ti,nu,1)
      call T1WT3ijk(i,k,j,NO,NU,T1,VOE,ti,t3)
      call trant3(ti,nu,4)
      call T1WT3ijk(j,i,k,NO,NU,T1,VOE,ti,t3)
      call trant3(ti,nu,1)
      call T1WT3ijk(j,k,i,NO,NU,T1,VOE,ti,t3)
      call trant3(ti,nu,5)
      call T1WT3ijk(k,i,j,NO,NU,T1,VOE,ti,t3)
      call trant3(ti,nu,1)
      call T1WT3ijk(k,j,i,NO,NU,T1,VOE,ti,t3)
      return
      end
      SUBROUTINE T1WT3(NO,NU,TI,T1,VOE,t3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D,E,F
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION TI(NU,NU,NU),VOE(Nu,NU,No,NO),T1(Nu,No),t3(nu,nu,nu)
      DATA ZERO/0.0D+0/,TWO/2.0D+00/,HALF/0.5D+00/,FOUR/4.0D+0/
      call ienter(52)
      call trt1(no,nu,ti,t1)
      CALL RO2pph(1,NO,NU,TI,VOE)
      DO 10 I=1,NO
         i1=i-1
      DO 10 j=1,i1
         j1=j-1
      DO 10 k=1,j1
      call drt1wt3ijk(i,j,k,NO,NU,T1,VOE,ti,t3)
 10   CONTINUE
      DO 11 I=1,NO
         i1=i-1
      DO 11 j=1,i1
      call drt1wt3ij(i,j,NO,NU,T1,VOE,ti,t3)
 11   continue
      DO 12 j=1,NO
         j1=j-1
      DO 12 k=1,j1
      call drt1wt3jk(j,k,NO,NU,T1,VOE,ti,t3)
 12      continue
      call trt1(nu,no,ti,t1)
      call iexit(52)
      RETURN
      END
      SUBROUTINE smt3four(NU,T3,V3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D,E,F
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension t3(nu,nu,nu),v3(nu,nu,nu)
      data two/2.0d+0/
      do 1 a=1,nu
      do 1 b=1,nu
      do 1 c=1,nu
         t3(a,b,c)=(v3(a,b,c)-v3(b,a,c))*TWO-v3(a,c,b)+v3(b,c,a)
 1    continue
      return
      end
      SUBROUTINE T2FT3(NO,NU,T2,T3,FPH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION T2(NU,NU,NO,NO),T3(NU,NU,NU),FPH(NU,NO)
      call ienter(53)
      CALL RDFPH(NO,NU,T3,FPH)
      DO 500 I=1,NO
      DO 500 J=1,NO
      DO 400 M=1,NO
      IF (I.EQ.J.AND.J.EQ.M)GOTO 400
      CALL RDVT3ONW(I,J,M,NU,T3)
      CALL SYMT3(T3,NU,6)
      CALL MATMUL(T3,FPH(1,M),T2(1,1,I,J),NU2,1,NU,0,0)
 400  CONTINUE
 500  CONTINUE
      call iexit(53)
      RETURN
      END
      subroutine drt2wt3(no,nu,t2,t)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical print
      common/flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/totmem/nmem
      common/itrat/icycle,mx,icn
      dimension t2(1),t(1)
      print=iflags(1).gt.10
      i1=1            !t3
      i2=i1+nu3       !vm
      i3=i2+no3u      !tmp
      i4=i3+nu3       !vepak
      i5=i4+nou3      !ive
      itot=i5+nou3/2+1
      if(icycle.eq.1.and.print)write(6,99)nmem,itot
 99   format('Space usage in drt2wt3: available - ',i8,'   used - ',i8)
      call t2wt3(no,nu,t2,t(i1),t(i2),t(i3),t(i4),t(i5))
      return
      end
      SUBROUTINE T2WT3(NO,NU,t2,t3,vm,tmp,vepak,ive)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D,E,F
      LOGICAL SDT
      common/wpak/nfr(30),nsz(30)
      common/pak/intg,intr
      COMMON/NEWOPT/NOPT(6)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION T2(NU,NU,NO,NO),T3(NU,NU,NU),VM(NO,NU,NO,NO),tmp(1),
     &vepak(1),ive(1)
      DATA TWO/2.0D+0/,half/0.5d+0/
      call ienter(56)
      SDT=NOPT(1).GT.1
      I10=0
      I11=1
      IF(SDT)THEN
         I10=10
         CALL RDVEM4(I10,NO,NU,T3,vepak)
         nsz(8)=nnve
         I11=11
      else
         call rdvem4(0,no,nu,t3,vepak)
      ENDIF
      CALL RDVEM1(I11,NU,NO,T3,VM)
      CALL TRANMD(VM,Nu,No,NO,NO,23)
      call ienter(129)
      DO 11 I=1,NO
      i1=i-1
      DO 11 J=1,i1
      j1=j-1
      DO 11 M=1,j1
      call rdvt3onw(i,j,m,nu,t3)
      call drt2wt3ijk(j,m,i,no,nu,t2,t3,vm,tmp,vepak,ive)
   11 CONTINUE
      DO 12 I=1,NO
      i1=i-1
      DO 12 J=1,i1
      call rdvt3onw(i,j,j,nu,t3)
      call drt2wt3ij(i,j,no,nu,t2,t3,vm,tmp,vepak,ive)
 12   CONTINUE
      DO 13 J=1,no
      j1=j-1
      DO 13 M=1,j1
      call rdvt3onw(j,j,m,nu,t3)
      call drt2wt3jm(j,m,no,nu,t2,t3,vm,tmp,vepak,ive)
 13   CONTINUE
      call iexit(129)
      CALL INSITU(NU,NU,NO,NO,T3,T2,13)
      call iexit(56)
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
      SUBROUTINE WO2HPP(NO,NU,TI,T2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER A,B
      COMMON/NEWCCSD/NTT2
      DIMENSION TI(NU,NU,NO),T2(NO,NU,NU,NO)
      NOU2=NO*NU*NU
      NLAST=NO+2
      DO 1 I=1,NO
      DO 2 J=1,NO
      DO 2 A=1,NU
      DO 2 B=1,NU
      TI(A,B,J)=T2(I,A,B,J)
 2    CONTINUE
      IASV=NLAST+I
      WRITE(NTT2,REC=IASV)TI
 1    CONTINUE
      RETURN
      END

      SUBROUTINE WRVEMT2(INO,NO,NU,TI,VEM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER A,B,C
      COMMON/NEWCCSD/NTT2
      DIMENSION TI(NU,NU,NU),VEM(NU,NU,NU,NO)
      DATA ZERO/0.0D+0/
      nou3=no*nu*nu*nu
      NNO=NO
      NNU=NU
      IF (INO.EQ.1) THEN
      NNO=NU
      NNU=NO
      ENDIF
      NLAST=3*NNO+2
      DO 2 I=1,NO
      IASV=NLAST+NU*INO+I
      DO 1 A=1,NU
      DO 1 B=1,NU
      DO 1 C=1,NU
      TI(A,B,C)=VEM(A,B,C,I)
 1    CONTINUE
      WRITE(NTT2,REC=IASV)TI
 2    CONTINUE
      RETURN
      END
      SUBROUTINE WRVOET2(NO,NU,TI,VOE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER A,B
      COMMON/NEWCCSD/NTT2
      DIMENSION TI(NU,NU,NO),VOE(NO,NU,NU,NO)
      NLAST=5*NO+2*NU+2
      DO 2 I=1,NO
      DO 1 A=1,NU
      DO 1 B=1,NU
      DO 1 J=1,NO
      TI(A,B,J)=VOE(I,A,B,J)
 1    CONTINUE
      IASV=NLAST+I
      WRITE(NTT2,REC=IASV)TI
 2    CONTINUE
      RETURN
      END
      SUBROUTINE WRINT2(NTT2,NLAST,NO,NOU2,VOE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION VOE(NOU2,NO)
      DO 2 I=1,NO
         call wr(ntt2,nlast+i,nou2,voe(1,i))
 2    CONTINUE
      RETURN
      END
      SUBROUTINE WRVOOT2(NO,NU,TI,VOE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER A,B
      COMMON/NEWCCSD/NTT2
      DIMENSION TI(NU,NU,NO),VOE(NO,NU,NU,NO)
      NLAST=6*NO+2*NU+2
      DO 2 I=1,NO
      DO 1 A=1,NU
      DO 1 B=1,NU
      DO 1 J=1,NO
      TI(A,B,J)=VOE(I,A,B,J)
 1    CONTINUE
      IASV=NLAST+I
      WRITE(NTT2,REC=IASV)TI
 2    CONTINUE
      RETURN
      END
      SUBROUTINE WRVOT2(INO,NO,NU,TI,V)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER A,B,C
      COMMON/NEWCCSD/NTT2
      DIMENSION TI(NU,NU,NU),V(NU,NU,NU,NU)
      DATA ZERO/0.0D+0/
      NNO=NO
      NNU=NU
      IF (INO.EQ.1) THEN
      NNO=NU
      NNU=NO
      ENDIF
      NLAST=4*NNO+NNU+2
      DO 2 I=1,NU
      IASV=NLAST+NNU*INO+I
      DO 1 A=1,NU
      DO 1 B=1,NU
      DO 1 C=1,NU
      TI(A,B,C)=V(A,B,C,I)
 1    CONTINUE
      WRITE(NTT2,REC=IASV)TI
 2    CONTINUE
      RETURN
      END
      subroutine mmf(ni,nj,nk,nn,a,b,ib,c,iadd,iplus)
      implicit double precision(a-h,o-z)
      dimension a(ni,1),b(1),c(ni,1),ib(1)
      data onem/-1.0d+0/
      call ienter(127)
      nij=ni*nj
      if(iplus.eq.1)call vecmul(c,nij,onem)
      do 20 jk=1,nn
         k=mod(ib(jk),nk)
         j=ib(jk)/nk+1
         if(k.ne.0)goto 11
         k=nk
         j=j-1
 11      continue
         x=b(jk)
         do 10 i=1,ni
            c(i,j)=c(i,j)+a(i,k)*x
 10      continue
 20   continue
      if(iplus.eq.1)call vecmul(c,nij,onem)
      call iexit(127)
      return
      end
      subroutine tr(ni,nj,a,b)
      implicit double precision(a-h,o-z)
      dimension a(1),b(1)
      nij=ni*nj
      call dcopy(nij,a,1,b,1)
      do 10 i1=1,nij
         i=mod(i1,ni)
         j=i1/ni+1
         if (i.eq.0)then
            i=ni
            j=j-1
         endif
         ji=(i-1)*nj+j
         a(ji)=b(i1)
 10   continue
      return
      end
      subroutine scattr(n,nn,b,ib)
      implicit double precision(a-h,o-z)
      dimension b(*),ib(*)
      data zero/0.0d+0/
      call ienter(126)
      id1=n-ib(nn)
      do 10 i=nn,1,-1
         ibi=ib(i)
         b(ibi)=b(i)
         do 11 j=1,id1
            b(ibi+j)=zero
 11      continue
         id1=ib(i)-ib(i-1)-1
 10   continue
      b(ib(1))=b(1)
      ib1=ib(1)-1
      do 12 i=1,ib1
         b(i)=zero
 12   continue
      call iexit(126)
      return
      end
      subroutine store(n,nn,b,ib)
      implicit double precision(a-h,o-z)
      dimension b(*),ib(*)
      data tr/0.1d-12/
      call ienter(125)
      in=0
      do 10 i =1,n
         x=b(i)
         if (dabs(x).gt.tr)then
            in=in+1
            b(in)=x
            ib(in)=i
         endif
 10   continue
      nn=in
      call iexit(125)
      return
      end
      subroutine storei(ioff,nj,ni,nn,b,ib)
      implicit double precision(a-h,o-z)
      common/isizes/ive(50)
      dimension b(*),ib(*)
      data tr/0.1d-12/
      call ienter(125)
      in=0
      ive(ioff+1)=1
      i=0
      do 10 ii =1,ni
         iveii=0
         do 9 ij =1,nj
            i=i+1
            x=b(i)
            if (dabs(x).gt.tr)then
               iveii=iveii+1
               in=in+1
               b(in)=x
               ib(in)=i
            endif
 9       continue
         ive(ioff+ii+1)=ive(ioff+ii)+iveii
 10   continue
      nn=in
 99   format(6f12.8)
      call iexit(125)
      return
      end
            
         
      subroutine rrv4(irec,n,v)
      implicit double precision(a-h,o-z)
      common/newio/nvv,nhh,npp,nvo,nt3,no1,nu1
      dimension v(n)
      read(nvv,rec=irec)v
      return
      end
      subroutine tranmdpak(iv,n1,n2,n3,n4,n,itransposition)
      implicit double precision (a-h,o-z)
      dimension iv(n)
      nu3=n1*n2*n3
      nu2=n1*n2
      nu=n1
      do 10 ii=1,n
         ijkl=iv(ii)
         ijk=mod(ijkl,nu3)
         l=ijkl/nu3+1
         if(ijk.eq.0)then
            ijk=nu3
            l=l-1
         endif
         ij=mod(ijk,nu2)
         k=ijk/nu2+1
         if(ij.eq.0)then
            ij=nu2
            k=k-1
         endif
         i=mod(ij,nu)
         j=ij/nu+1
         if(i.eq.0)then
            i=nu
            j=j-1
         endif
         if (itransposition.eq.12)goto 12
         if (itransposition.eq.13)goto 13
         if (itransposition.eq.14)goto 14
         if (itransposition.eq.23)goto 23
         if (itransposition.eq.24)goto 24
         if (itransposition.eq.34)goto 34
 12      continue
         iv(ii)=(l-1)*nu3+(k-1)*nu2+(i-1)*nu+j
         goto 10
 13      continue
         iv(ii)=(l-1)*nu3+(i-1)*nu2+(j-1)*nu+k
         goto 10
 14      continue
         iv(ii)=(i-1)*nu3+(k-1)*nu2+(j-1)*nu+l
         goto 10
 23      continue
         itest =(l-1)*nu3+(k-1)*nu2+(j-1)*nu+i
         iv(ii)=(l-1)*nu3+(j-1)*nu2+(k-1)*nu+i
c         write(6,*)'old,new',ijkl,itest
         goto 10
 24      continue
         iv(ii)=(j-1)*nu3+(k-1)*nu2+(l-1)*nu+i
         goto 10
 34      continue
         iv(ii)=(k-1)*nu3+(l-1)*nu2+(j-1)*nu+i
         goto 10
 10   continue
      return
      end
      subroutine vripak(nv4,nint,n,vec,ivec)
      implicit double precision (a-h,o-z)
      common/wpak/nfr(30),nsz(30)
      dimension vec(*),ivec(*)
      common/isizes/ive(50)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/lnt/lrecb,lrecw
      lrec=lrecw
      nrec=n/lrec+1
      noff=0
      if (nrec.eq.1)then
         irec=nfr(nint)
         call ww(irec,n,vec,ivec)
      else
      do 10 ir=1,nrec-1
         irec=nfr(nint)+ir-1
         call ww(irec,lrec,vec(noff+1),ivec(noff+1))
      noff=noff+lrec
 10   continue
      irest=n-(nrec-1)*lrec
      irec=nfr(nint)+nrec-1
         call ww(irec,irest,vec(noff+1),ivec(noff+1))
      endif
      isafe=nrec/100+1
      if(nint.ne.15)then
         nfr(nint+1)=nfr(nint)+nrec+isafe
      endif
      return
      end
      subroutine ww(ir,n,v,iv)
      implicit double precision (a-h,o-z)
      dimension v(n),iv(n)
      common/pak/nv4,intr
      common/lnt/lrecb,lrecw
c      common/pak/nv4,lrecb,lrecw
      write(nv4,rec=ir)v,iv
      return
      end
      subroutine rr(ir,n,v,iv)
      implicit double precision (a-h,o-z)
      dimension v(n),iv(n)
      common/pak/nv4,intr
      common/lnt/lrecb,lrecw
c      common/pak/nv4,lrecb,lrecw
      read(nv4,rec=ir)v,iv
      return
      end
      subroutine reapak(nv4,nint,vec,ivec)
      implicit double precision (a-h,o-z)
      dimension vec(1),ivec(1)
      common/wpak/nfr(30),nsz(30)
      common/lnt/lrecb,lrecw
      n=nsz(nint)
      lrec=lrecw
      nrec=n/lrec+1
      noff=0
      if (nrec.eq.1)then
         irec=nfr(nint)
         call rr(irec,n,vec,ivec)
      else
      do 10 ir=1,nrec-1
         irec=nfr(nint)+ir-1
         call rr(irec,lrec,vec(noff+1),ivec(noff+1))
      noff=noff+lrec
 10   continue
      irest=n-(nrec-1)*lrec
      irec=nfr(nint)+nrec-1
         call rr(irec,irest,vec(noff+1),ivec(noff+1))
      endif
      return
      end
      SUBROUTINE ADT12(igoto,NO,NU,T1,T2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer a,b
      dimension t1(no,nu),t2(no,no,nu,nu)
      do 10 i=1,no
      do 10 j=1,no
      do 10 a=1,nu
      do 10 b=1,nu
         goto(1,2)igoto
 1       continue
         t2(i,j,a,b)=t2(i,j,a,b)+t1(i,a)*t1(j,b)
         goto 10
 2       continue
         t2(i,j,a,b)=t2(i,j,a,b)-t1(i,a)*t1(j,b)
 10   continue
      return
      end
      subroutine drintquat2(no,nu,t)
      implicit double precision (a-h,o-z)
      logical print
      common/flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/totmem/nmem
      common/itrat/icycle,mx,icn
      dimension t(*)
      print=iflags(1).gt.10
      i1=1
      i2=i1+nou
      i3=i2+nu3
      i4=i3+no2u2
      i5=i4+no2u2
      i6=i5+nou3
      i7=i6+nou3/2+1
      itot=i7+nu4
c      write(6,*)'nmem',nmem
      if(icycle.eq.1.and.print)write(6,99)nmem,itot
 99   format('Space usage in intquat2: available - ',i8,'   used - ',i8)
      call intquat2(no,nu,t,t(i2),t(i3),t(i4),t(i5),t(i6),t(i7))
      return
      end
      SUBROUTINE INTQUAT2(NO,NU,O1,TI,O2,T2,ve,ive,v)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical sdt
      common/newopt/nopt(6)
      common/pak/intg,intr
      common/wpak/nfr(30),nsz(30)
      common/newccsd/ntt2
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION O1(1),TI(1),O2(no,no,nu,nu),T2(1),VE(nu,nu,no,nu),
     *ive(1),v(nu,nu,nu,nu)
      DATA TWO/2.0D+00/,HALF/0.5D+00/,ONE/1.0D+0/
c      write(6,*)'nsz1-6',(nsz(iik),iik=1,8)
c      write(6,*)'intquat2 entered'
      sdt=nopt(1).ge.2
      call ienter(11)
c pppp intermediate
      if(.not.sdt)goto 300
      CALL RDOV4(0,NO,NU,TI,V)
      CALL RDT1HP(NO,NU,O1)
      CALL RDVEM3(0,NO,NU,TI,VE)
      CALL RDVEM3(-10,NO,NU,TI,VE)
      CALL VECMUL(VE,NOU3,HALF)
      CALL RO2pph(0,NO,NU,TI,T2)
      CALL RO2hhp(1,NO,NU,TI,o2)
      do 102 i=1,nu
         call trant3(v(1,1,1,i),nu,1)
         if(sdt)CALL MATMUL(t2,o2(1,1,1,i),V(1,1,1,i),NU2,NU,NO2,0,0)
         call trant3(v(1,1,1,i),nu,1)
         call matmul(ve(1,1,1,i),o1,v(1,1,1,i),nu2,nu,no,0,1)
         call trt1(no,nu,ti,o1)
         call insitu(nu,nu,no,nu,ti,ve,13)
         call tranmd(ve,no,nu,nu,nu,24)
         call matmul(o1,ve(1,1,1,i),v(1,1,1,i),nu,nu2,no,0,1)
         call trt1(nu,no,ti,o1)
         call tranmd(ve,no,nu,nu,nu,24)
         call insitu(no,nu,nu,nu,ti,ve,13)
         call trant3(v(1,1,1,i),nu,1)
 102  continue
      if(.not.sdt)call mtrans(v,nu,3)
      CALL WROV4(0,NO,NU,TI,V)
 300  continue
c      write(6,*)'pppp done'
c pppp intermediate done
      CALL RO2PPH(0,NO,NU,TI,O2)
      CALL TRANMD(O2,NU,NU,NO,NO,34)
      CALL RO1(NOU,1,O1)
c hhhh intermediate
      call ienter(122)
      CALL RDOV4(1,NU,NO,TI,T2)
      CALL MTRANS(T2,NO,8)
c 1
      CALL RO2HHP(1,NO,NU,TI,VE)
      CALL TRANMD(VE,NO,NO,NU,NU,34)
      CALL MATMUL(VE,O2,T2,NO2,NO2,NU2,0,0)
      CALL VECMUL(T2,NO4,HALF)
      call tranmd(T2,no,no,no,no,231)
c 2
      CALL RDVEM4(1,NU,NO,TI,VE)
      CALL RDVEM4(-11,NU,NO,TI,VE)
      CALL VECMUL(VE,NO3U,HALF)
      CALL MATMULT(VE,O1,T2,NO3,NO,NU,0,0)
      call symv(T2,no)
      call tranmd(T2,no,no,no,no,312)
c     attention HHHH INT  written on disk as T2(m,n,i,j) (m-i,n-j)
       CALL WROV4(1,NU,NO,TI,T2)
      call iexit(122)
c hhhh intermediate done; compute voe intermediates
c      write(6,*)'hhhh done'
      call ienter(123)
      CALL RO2HPP(1,NO,NU,TI,T2)
      CALL RDVEM4(10,NO,NU,TI,VE)
      CALL TRANMD(VE,NU,NU,NU,NO,12)
c 1
      CALL MATMUL(O1,VE,T2,NO,NOU2,NU,0,0)
      CALL TRANMD(T2,NO,NU,NU,NO,1234)
      CALL TRANSQ(T2,NOU)
C now we have phhp with last two active
      CALL RDVEM4(1,NU,NO,TI,VE)
c 2
      call trt1(no,nu,ti,o1)
      CALL MATMUL(O1,VE,T2,NU,NO2U,NO,0,1)
      CALL TRANSQ(T2,NOU)
      CALL TRANMD(T2,NO,NU,NU,NO,1234)
      call wrvoep(1,no,nu,ti,t2)
C now we have hpph with last two active
      CALL RO2HPH(0,NO,NU,TI,O2)
c      write(6,*)'entering reapak 3'
      call reapak(intg,3,ve,ive)
c      write(6,*)'exiting reapak 3'
      call scattr(no2u2,nsz(3),ve,ive)
      CALL VeCMUL(O2,NO2U2,HALF)
      CALL SYMT21(o2,no,nu,no,nu,13)
c 3
      call matmul(o2,ve,t2,nou,nou,nou,0,0)
      CALL desm21(o2,no,nu,no,nu,13)
c 4
      CALL TRANMD(ve,NO,NU,NU,NO,23)
      call matmul(o2,ve,t2,nou,nou,nou,0,1)
      call tranmd(t2,no,nu,nu,no,1234)
      call store(no2u2,nnint,t2,ive)
      nsz(10)=nnint
      call vripak(intg,10,nnint,t2,ive)
      call scattr(no2u2,nnint,t2,ive)
      call rovoep(1,no,nu,ti,ve)
      call tranmd(t2,no,nu,nu,no,1234)
      call vecsub(t2,ve,no2u2)
      call vecadd(ve,t2,no2u2)
c now t2 zawiera wynik dwoch ostatnich matmuli
      call vecadd(ve,t2,no2u2)
      call wrvoe(no,nu,ti,ve)
      call iexit(123)
c      write(6,*)'voe done'
c now vo
      call ienter(124)
c 1
      CALL RO2HPP(2,NO,NU,TI,T2)
      CALL TRANSQ(T2,NOU)
      CALL RDVEM4(1,NU,NO,TI,VE)
      CALL TRANMD(VE,NO,NO,NO,NU,13)
      call ro1(nou,1,o1)
c 1
      CALL TMATMUL(O1,VE,T2,NU,NO2U,NO,0,1)
      CALL TRANSQ(T2,NOU)
      CALL TRANMD(T2,NO,NU,NU,NO,1234)
c tu zapamietujemy t2
      call wrvoep(2,no,nu,ti,t2)
      CALL TRANMD(O2,NO,NU,No,Nu,13)
      call reapak(intg,3,ve,ive)
      call scattr(no2u2,nsz(3),ve,ive)
      CALL TRANMD(ve,NO,NU,Nu,No,23)
      call matmul(o2,ve,t2,nou,nou,nou,1,1)
c 2
      call rovoep(2,no,nu,ti,o2)
      call vecadd(o2,t2,no2u2)
      CALL TRANMD(o2,NO,NU,NU,NO,23)
c 3
      CALL RDVEM4(10,NO,NU,TI,VE)
      CALL TRANMD(VE,NU,NU,NU,NO,13)
      call ro1(nou,1,o1)
      call matmul(o1,ve,o2,no,nou2,nu,0,0)
      call tranmd(o2,no,nu,nu,no,14)
      call store(no2u2,nnint,o2,ive)
      nsz(11)=nnint
      call vripak(intg,11,nnint,o2,ive)
      call scattr(no2u2,nnint,o2,ive)
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call tranmd(o2,no,nu,nu,no,1234)
      call vecadd(o2,t2,no2u2)
      call wrvo(no,nu,ti,o2)
      call ro2hpp(2,no,nu,ti,t2)
      call rdvem4(0,no,nu,ti,ve)
      call tranmd(ve,nu,nu,nu,no,13)
      call matmul(o1,ve,t2,no,nou2,nu,0,0)
      call wrvoep(3,no,nu,ti,t2)
      call rovoep(3,no,nu,ti,t2)
c      write(6,*)'vo done'
      call iexit(124)
      call iexit(11)
      RETURN
      END
      SUBROUTINE CHKSUM(A,T,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T(N)
      DATA ZERO/0.0D+0/
c      return
      X=ZERO
      Y=ZERO
      ZY=ZERO
      XT=ZERO
      DO 10 I=1,N
      X=X+DABS(T(I))
      XT=XT+T(I)
 10   CONTINUE
      DO 20 I=2,N,2
      Z=T(I-1)+T(I)
      Y=Y+Z*Z
 20   CONTINUE
      DO 30 I=5,N,5
      Z=T(I-4)+T(I-3)+T(I-2)+T(I-1)+T(I)    
      ZY=ZY+Z*Z
 30   CONTINUE
 100  FORMAT(1X,A8,'ABS:',D15.10)
C      WRITE(6,100)A,X
 200  FORMAT(A8,' ',D16.11,' ',D16.11,' ',D16.11,' ',D16.10)
      WRITE(6,200)A,X,Y,ZY,XT
      RETURN
      END
