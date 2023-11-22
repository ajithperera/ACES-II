      subroutine drt4(cr)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON /newio/nto2,nthh,ntpp,ntvo,ntt3,no,nu
      common/info/nocco(2),nvrto(2)
      dimension cr(1)
      common/actadres/iacta
      no=nocco(1)
      nu=nvrto(1)
      iacta=1
      call setsiz(no,nu)
      call opfit4
      nnev=2*(no+nu+1)
      call zeroma(cr,1,nu3)
      ii0=1               !o3
      ii1=ii0+nu3+100   !w1
      ii2=ii1+nu3+100   !w2
      ii3=ii2+nu3+100   !t3
      ii4=ii3+nu3+100   !isym
      ii5=ii4+400+10000    !iadw2
      ii6=ii5+8+10000     !iadblk
      ii7=ii6+8+10000     !len
      call gett3s(no,nu,cr,cr(ii1),cr(ii2),cr(ii3),cr(ii4),cr(ii5),
     &cr(ii6),cr(ii7))
      call t4(no,nu,cr,cr(nnev+1))
      return
      end 
      SUBROUTINE t4(no,nu,ev,cr)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
#include "machsp.com"
      logical smlv,prt
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/iopos/icrsiz,icsc1(4)
      common/flags/iflags(100)
      common/en51/ec1
      common/avlspace/maxcor
      common/syminf/nstart,nirrep,irrepa(255),irrepb(255),dirprd(8,8)
      common/icdacc/idaccr,idaccm
      DIMENSION ev(1),cr(1)
      equivalence(icllvl,iflags(2))
      call zclock('start   ',0)
      prt=iflags(1).ge.10
      smlv=icllvl.eq.11.or.icllvl.eq.13.or.icllvl.eq.14.or.
     *icllvl.eq.15.or.icllvl.eq.19.or.
     *icllvl.eq.21.or.icllvl.eq.22
 99   format(/,77('-'),/,22X,'NONITERATIVE T4 CONTRIBUTION',/,77('-'),/)
      write(6,99)
      ec1=0.0d+0
      call getrec(20,'JOBARC','SCFEVALA',iintfp*(no+nu),ev)
      imno3u=no3u
      if(nu3.gt.no3u)imno3u=nu3
      imnu3=no2u2
      if(nu3.gt.no2u2)imnu3=nu3
      ito1=6*nu2
      ito2=2*no2u2+imnu3+nu3
      ito3=8*no2u2+imnu3
      ito4=5*no2u2+3*imnu3
      ito5=4*no2u2+imnu3+no4
      ito6=2*no2u2+imnu3+nou3+nou+no3u
      if(nirrep.gt.1)write(6,996)ito1
      if(prt)then
      write(6,991)ito2
      write(6,992)ito3
      write(6,993)ito4
      write(6,994)ito5
      write(6,995)ito6
      endif
      maxcor=(icrsiz-mod(icrsiz,2))/2
      maxcor=iflags(36)/2
      imaxs=ito1
      if(ito2.gt.imaxs)imaxs=ito2
      if(ito3.gt.imaxs)imaxs=ito3
      if(ito4.gt.imaxs)imaxs=ito4
      if(ito5.gt.imaxs)imaxs=ito5
      if(ito6.gt.imaxs)imaxs=ito6
      if(imaxs.gt.maxcor)then
      write(6,997)imaxs,maxcor
      endif
 997  format('Insuficient memory: needed -',i10,',  available -',i10) 
      call zeroma(cr,1,no2u2)
      io0=1
      io1=io0+no2u2 +100000 
      io2=io1+nou   +100000 
      io3=io2+imnu3 +100000 
      io4=io3+nou3  +100000 
      io5=io4+no2u2 +100000 
      io6=io5+no2u2 +100000
      io7=io6+no2u2 +100000
      io8=io7+no2u2 +100000
      ito=io5+no3u
      call t2wwt3(no,nu,cr(io1),cr(io2),cr(io3),cr(io4),cr(io5),cr,ev,
     &ev(no+1),cr(io6),cr(io7),cr(io8))
      if(icllvl.eq.18)call auxioo(.TRUE.)
      close(unit=45,status='delete')
      io0=1
      io1=io0+nu2
      io2=io1+nu2
      io3=io2+nu2
      io4=io3+nu2
      io5=io4+nu2
      ito=io5+nu2  
      call rdsymv(nu,cr,cr(io1),cr(io2),cr(io3),cr(io4),cr(io5),233)
      if(smlv)then
         call symvt(nu,cr,cr(io1),cr(io2),cr(io3))
      endif
      istat=ishell('test -f MOABCD')
      io0=1              !o2
      io1=io0+no2u2      !a
      call dgr22sym(no,nu,cr,cr(io1),ev,ev(no+1))
      io0=no2u2+1           !ti
      io1=io0+nu3     !o2
      io2=io1+no2u2   !vo2
      io3=io2+no2u2   ! v
      ito=io3+imnu3   ! ito
      call ro2hpp(0,no,nu,cr(io0),cr(io3))
      call rank1(no,nu,cr(io0),cr(io1),cr(io2),cr(io3),ev,ev(no+1))
      io0=1           ! t2n
      io1=io0+no2u2   ! fpp
      io2=io1+no2u2   ! ti
      io3=io2+no2u2   ! o2
      io4=io3+no2u2   ! t2
      io5=io4+no2u2   ! vo2
      io6=io5+no2u2   ! v
      io7=io6+imnu3   ! voe
      io8=io7+no2u2   ! vo
      ito=io8+no2u2
      call zeroma(cr,1,ito)
      call ro2hpp(0,no,nu,cr(io1),cr(io3))
      call t2hpt2sq(no,nu,cr(io1),cr(io2),cr(io3),cr(io4),cr(io5),
     &                       cr(io6),cr(io7),cr(io8), cr)
      io1=no2u2+1
      io2=no2u2+io1
      call symetr(cr,no,nu)
      call engy1(no,nu,cr(io1),cr,cr(io2),'t2hpt2sq')
      call zeroma(cr,1,no2u2)
      io0=1             !t2n
      io1=io0+no2u2     !ti
      io2=io1+imnu3     !o2
      io3=io2+no2u2     !t2
      io4=io3+imnu3     !voe
      io5=io4+no2u2     !v
      io6=io5+imnu3     !o4
      io7=io6+imnu3     !vo2
      ito=io7+no2u2
      call t2ppt2sq(no,nu,cr(io1),cr(io2),cr(io3),cr(io4),
     &cr(io5),cr(io6),cr(io7),cr,ev,ev(no+1))
      io1=no2u2+1
      io2=no2u2+io1
      call symetr(cr,no,nu)
      call engy1(no,nu,cr(io1),cr,cr(io2),'t2ppt2sq')
      call zeroma(cr,1,no2u2)
      io0=1             !t2n
      io1=io0+no2u2     !ti
      io2=io1+imnu3     !o2
      io3=io2+no2u2     !t2
      io4=io3+no2u2     !vo2
      io5=io4+no2u2     !v
      ito=io5+no4
      call t2hht2sq(no,nu,cr(io1),cr(io2),cr(io3),cr(io4),cr(io5),cr)
      io0=1             ! t2n
      io1=io0+no2u2     ! ti
      io2=io1+nou2      ! vo
      ito=io2+no2u2
      call symetr(cr,no,nu)
      call engy(no,nu,cr(io1),cr,cr(io2),'t2hht2sq')
 991  format('Fact. 5th order quads: rank1  requires',i8,' words')
 992  format('Fact. 5th order quads: t2hpt2 requires',i8,' words')
 993  format('Fact. 5th order quads: t2ppt2 requires',i8,' words')
 994  format('Fact. 5th order quads: t2hht3 requires',i8,' words')
 995  format('Fact. 5th order quads: t2wwt3 requires',i8,' words')
 996  format('Fact. 5th order quads: rdsym  requires',i8,' words')
      close(unit=35,status='delete')
      close(unit=42,status='delete')
      close(unit=43,status='delete')
      close(unit=44,status='delete')
      close(unit=46,status='delete')
      RETURN
      END

      SUBROUTINE T1WT3INT(NO,NU,TI,T1,VOE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D,E,F
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION TI(NU,NU,NU),VOE(NU,NU,NO,NO),T1(Nu,No)
      DATA ZERO/0.0D+0/,TWO/2.0D+00/,HALF/0.5D+00/,FOUR/4.0D+0/
      call insitu(no,nu,nu,no,ti,voe,13)
      call tranmd(voe,nu,nu,no,no,12)
      call symt21(voe,nu,nu,no,no,12)
      call zeroma(ti,1,nu3)
      DO 10 I=1,NO
      DO 10 M=1,NO
      DO 10 N=1,NO
      IF(I.EQ.M.AND.M.EQ.N)goto 10
      CALL RDVT3ONW(I,M,N,NU,TI)
      call matmul(ti,voe(1,1,m,n),t1(1,i),nu,1,nu2,0,0)
      call trant3(ti,nu,2)
      call matmul(ti,voe(1,1,m,n),t1(1,i),nu,1,nu2,0,1)
 10   continue
      call trt1(nu,no,ti,t1)
      call desm21(voe,nu,nu,no,no,12)
      call tranmd(voe,nu,nu,no,no,12)
      call insitu(nu,nu,no,no,ti,voe,13)
cc      call iexit(52)
      RETURN
      END
      SUBROUTINE T2FT3INT(NO,NU,T2,T3,FPH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION T2(NU,NU,NO,NO),T3(NU,NU,NU),FPH(NU,NO)
      DO 500 I=1,NO
      DO 500 J=1,NO
      DO 400 M=1,NO
      IF (I.EQ.J.AND.J.EQ.M)GOTO 400
      CALL RDVT3ONW(I,J,M,NU,T3)
      CALL SYMT3(T3,NU,6)
      CALL MATMUL(T3,FPH(1,M),T2(1,1,I,J),NU2,1,NU,0,0)
 400  CONTINUE
 500  CONTINUE
      RETURN
      END
      SUBROUTINE T2WT3P(NO,NU,T2,T3,VE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,nou,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION T2(NU,NU,NO,NO),T3(NU,NU,NU),VE(NU,NU,NU,NO)
      DATA TWO/2.0D+0/
      DO 11 I=1,NO
      DO 11 J=1,NO
      DO 12 M=1,NO
      IF (I.EQ.J.AND.I.EQ.M) GO TO 12
      CALL RDVT3ONW(I,M,J,NU,T3)
      CALL SYMT3(T3,NU,1)
      CALL MATMUL(VE(1,1,1,M),T3,T2(1,1,J,I),NU,NU,NU2,0,0)
   12 CONTINUE
   11 CONTINUE
      RETURN
      END
      SUBROUTINE T2WT3H(NO,NU,T2,T3,VM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,nou,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION T2(NU,NU,NO,NO),T3(NU,NU,NU),VM(NO,NU,NO,NO)
      DATA TWO/2.0D+0/
      DO 13 J=1,NO
      DO 14 M=1,NO
      DO 14 N=1,NO
      IF (M.EQ.N.AND.M.EQ.J) GO TO 14
      CALL RDVT3ONW(M,N,J,NU,T3)
      CALL VECMUL(T3,NU3,TWO)
      CALL MATMUL(VM(1,1,M,N),T3,T2(1,1,1,J),NO,NU2,NU,0,1)
      CALL RDVT3ONW(N,M,J,NU,T3)
      CALL MATMUL(VM(1,1,M,N),T3,T2(1,1,1,J),NO,NU2,NU,0,0)
      CALL TRANMD(T2,NO,NU,NU,NO,23)
      CALL RDVT3ONW(J,M,N,NU,T3)
      CALL MATMUL(VM(1,1,M,N),T3,T2(1,1,1,J),NO,NU2,NU,0,0)
      CALL TRANMD(T2,NO,NU,NU,NO,23)
   14 CONTINUE
   13 CONTINUE
      RETURN
      END

      SUBROUTINE VMWT3(NO,NU,VM,T3,VO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,nou,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION VM(NO,NU,NO,NO),T3(NU,NU,NU),VO(NO,NU,NU,NO)
      DATA TWO/2.0D+0/
      DO 11 I=1,NO
      DO 11 J=1,NO
      DO 12 M=1,NO
      IF (I.EQ.J.AND.I.EQ.M) GO TO 12
      CALL RDVT3ONW(I,M,J,NU,T3)
      CALL SYMT3(T3,NU,1)
      CALL MATMUL(VO(1,1,1,M),T3,VM(1,1,J,I),NO,NU,NU2,0,0)
   12 CONTINUE
   11 CONTINUE
      RETURN
      END
      SUBROUTINE VEWT3(NO,NU,VE,T3,VO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,nou,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION VE(NU,NU,NU,NO),T3(NU,NU,NU),VO(NU,NU,NO,NO)
      DATA TWO/2.0D+0/
      DO 13 J=1,NO
      DO 14 M=1,NO
      DO 14 N=1,NO
      IF (M.EQ.N.AND.M.EQ.J) GO TO 14
      CALL RDVT3ONW(M,N,J,NU,T3)
      CALL VECMUL(T3,NU3,TWO)
      CALL MATMUL(VO(1,1,M,N),T3,VE(1,1,1,J),NU,NU2,NU,0,1)
      CALL RDVT3ONW(N,M,J,NU,T3)
      CALL MATMUL(VO(1,1,M,N),T3,VE(1,1,1,J),NU,NU2,NU,0,0)
      CALL TRANMD(VE,NU,NU,NU,NO,23)
      CALL RDVT3ONW(J,M,N,NU,T3)
      CALL MATMUL(VO(1,1,M,N),T3,VE(1,1,1,J),NU,NU2,NU,0,0)
      CALL TRANMD(VE,NU,NU,NU,NO,23)
   14 CONTINUE
   13 CONTINUE
      RETURN
      END
      subroutine opfit4
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON /newio/nto2,nthh,ntpp,ntvo,ntt3,no,nu
      nword=8
      nthh=42
      ntpp=43
      ntvo=44
      ntt3=45
      llt1=nword*8192
      llt2=nword*13724
      lo2=nword*nou2
      lvo=nword*nou2
      lhh=nword*no2
      lpp=nword*nu2
      lt3=68536
      close(39,status='delete')
      open(nthh,recl=lhh,access='direct',form='unformatted')
      open(ntpp,recl=lpp,access='direct',form='unformatted')
      open(ntvo,recl=lvo,access='direct',form='unformatted')
      open(ntt3,recl=lt3,access='direct',form='unformatted')
      open(46,recl=llt1 ,access='direct',form='unformatted')
      open(35,recl=llt1 ,access='direct',form='unformatted')
      open(36,recl=llt1 ,access='direct',form='unformatted')
      open(37,recl=llt2 ,access='direct',form='unformatted')
      open(38,recl=llt2 ,access='direct',form='unformatted')
      open(39,recl=llt2 ,access='direct',status='new',
     *form='unformatted')
      return
      end
      subroutine setsiz(nh,np)
      implicit double precision(a-h,o-z)
      COMMON/NN/NH2,NH3,NH4,NP2,NP3,NP4,NHP,NH2P,NH3P,NHP2,NHP3,NH2P2
      NH2=NH*NH
      NH3=NH*NH2
      NH4=NH*NH3
      NP2=NP*NP
      NP3=NP*NP2
      NP4=NP*NP3
      NHP=NH*NP
      NH2P=NH2*NP
      NH3P=NH3*NP
      NHP2=NH*NP2
      NHP3=NH*NP3
      NH2P2=NH2*NP2
      return
      end

