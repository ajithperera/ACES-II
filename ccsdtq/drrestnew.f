      subroutine drrest(iopt,no,nu,t)
      implicit double precision (a-h,o-z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/newopt/nopt(6)
      dimension t(1)
      i1=no2u2+1
      if(iopt.eq.0)then
      call reares(no,nu,t,t(i1))
      else
      call wrires(no,nu,t,t(i1)) 
      endif
      return
      end
      subroutine reapakres(iadr,nint,vec,ivec)
      implicit double precision (a-h,o-z)
      dimension vec(1),ivec(1)
      common/restart/nresf,nresl
      common/ressiz/isz(1000)
      n=nint
      lrecw=nresl/12
      lrec=lrecw
      nrec=n/lrec+1
      noff=0
      if (nrec.eq.1)then
         irec=iadr
         call rrres(irec,n,vec,ivec)
      else
      do 10 ir=1,nrec-1
         irec=iadr+ir-1
         call rrres(irec,lrec,vec(noff+1),ivec(noff+1))
         noff=noff+lrec
 10   continue
      irest=n-(nrec-1)*lrec
      irec=iadr+nrec-1
      call rrres(irec,irest,vec(noff+1),ivec(noff+1))
      endif
      iadr=irec+1
      return
      end

      subroutine reares(no,nu,t,it)
      implicit double precision (a-h,o-z)
      logical irest,t32
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/restart/nresf,nresl
      common/rstrt/irest,ccopt(10),t32
      common/ressiz/isz(1000)
      common/newopt/nopt(6)
      common/treshold/tsh
      dimension t(1),it(1)
      data zero/0.0d+0/
      read(10,rec=1)isz
      iad=2
      call flush(6)
      nnn=isz(1)
      call reapakres(iad,nnn,t,it)
      if (tsh.gt.zero)call select(nnn,t)
      call scattr(nou,nnn,t,it)
      call wo1(nou,1,t)
      nnn=isz(2)
      call reapakres(iad,nnn,t,it)
      if (tsh.gt.zero)call select(nnn,t)
      call scattr(no2u2,nnn,t,it)
      call wo2new(no,nu,it,t)
 99   continue
      return
      end
      subroutine rrres(ir,n,v,iv)
      implicit double precision (a-h,o-z)
      dimension v(n),iv(n)
      common/restart/nresf,nresl
      read(nresf,rec=ir)v,iv
      return
      end
      subroutine select(n,a)
      implicit double precision (a-h,o-z)
      common/treshold/tsh
      dimension a(n)
      data zero/0.0d+0/
      do 10 i=1,n
      if(dabs(a(i)).lt.tsh)a(i)=zero
  10  continue
      return
      end
      subroutine vripakres(iadr,nint,n,vec,ivec)
      implicit double precision (a-h,o-z)
      common/wpak/nfr(30),nsz(30)
      common/restart/nresf,nresl
      dimension vec(*),ivec(*)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      lrecw=nresl/12
      lrec=lrecw
      nrec=n/lrec+1
      noff=0
      if (nrec.eq.1)then
         irec=iadr
         call wwres(irec,n,vec,ivec)
      else
         do 10 ir=1,nrec-1
            irec=iadr+ir-1
            call wwres(irec,lrec,vec(noff+1),ivec(noff+1))
            noff=noff+lrec
 10      continue
         irest=n-(nrec-1)*lrec
         irec=iadr+nrec-1
         call wwres(irec,irest,vec(noff+1),ivec(noff+1))
      endif
      iadr=irec+1
      return
      end
      subroutine wrires(no,nu,t,it)
      implicit double precision (a-h,o-z)
      logical irest,t32
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/restart/nresf,nresl
      common/rstrt/irest,ccopt(10),t32
      common/ressiz/isz(1000)
      common/newopt/nopt(6)
      dimension t(1),it(1)
      iad=2
      call rdt1hp(no,nu,t)
      call store(nou,nnn,t,it)
      isz(1)=nnn
CSS      write(6,*)'nnn',nnn
      call flush(6)
      call vripakres(iad,nou,nnn,t,it)
      call ro2hpp(0,no,nu,it,t)
      call store(no2u2,nnn,t,it)
      isz(2)=nnn
      call vripakres(iad,no2u2,nnn,t,it)
      isz(1000)=0
      write(10,rec=1)isz
      return
      end
      subroutine wwres(ir,n,v,iv)
      implicit double precision (a-h,o-z)
      dimension v(n),iv(n)
      common/restart/nresf,nresl
      write(nresf,rec=ir)v,iv
      return
      end
      subroutine drmovo2(no,nu,t)
      implicit double precision (a-h,o-z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension t(*)
      i1=1
      i2=i1+nou2
      i3=i2+no2u2
      call movo2(no,nu,t(i1),t(i2),t(i3))
      return
      end
      SUBROUTINE MOVO2(NO,NU,TI,O2,T2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      COMMON /NEWT4/NT4,NTO4,LT4,NALL4,LL4
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/dift2/difmax
      COMMON /NEWCCSD/NTT2
      DIMENSION TI(NU,NU,NO),T2(NO,NU,NU,NO),O2(1)
      data onem/-1.0d+0/
      call ienter(18)
      nlast=no+2
      call ro2hpp(0,no,nu,ti,o2)
      IOFF=1
      DO 1 I=1,NO
      ias=nlast+i
      read(ntt2,rec=ias)ti
      IOFF=IOFF+NO
      IASV=I
       do 2 j=1,no
       do 2 a=1,nu
       do 2 b=1,nu
         t2(i,a,b,j)=ti(a,b,j)
 2     continue
      WRITE(NALL4,REC=IASV)TI
 1    CONTINUE
      CALL VADD(O2,T2,O2,NO2U2,ONEM)
      MAXDIF=ISAMAX(NO2U2,O2,1)
      DIFMAX=O2(MAXDIF)
      call iexit(18)
      RETURN
      END
      SUBROUTINE MOVT1(NO,NU,T1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NEWCCSD/NTT2
      DIMENSION T1(NO,NU)
      call ienter(17) 
      READ(NTT2,REC=2)T1
      WRITE(NTT2,REC=1)T1
      call iexit(17)
      RETURN
      END
      SUBROUTINE wrva0(NO,NU,v)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION v(1)
      common/nvtap/ntiter,nq2,nvt
      call zeroma(v,1,nu4)
      do 10 i=1,no
         do 11 j=1,i
            kk=it2(i,j)
            iasv=no3+(kk-1)*nu+1
            call wrva(iasv,nu,v)
 11      continue
 10   continue
      do 20 i=1,no
         do 21 j=1,no
            do 22 k=1,no
               kk=no2*(i-1)+no*(j-1)+k
               call wrgen(nvt,kk,nou2,v)
 22         continue
 21      continue
 20   continue
      RETURN
      END
