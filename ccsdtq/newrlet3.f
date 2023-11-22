      subroutine drnewrlet3(no,nu,t)
      implicit double precision (a-h,o-z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/newopt/nopt(6)
      dimension t(1)
      nl=it3(no,no,no-1)
      nlt4=it4(no,no,no-1,no-1)
      l1=nou2
      l2=no2u2+nou
      if(nopt(1).ge.2)then
      l1=nu3
      l2=no2u2+nou+nl*nu3
      endif
      if(nopt(1).eq.8)then
      endif
      io1=1
      io2=io1+l1
      io3=io2+l2
      io4=io3+l2
      call newrlet3(no,nu,t,t(io2),t(io3),t(io4))
      return
      end 
      subroutine newrlet3(no,nu,ti,o2,t2,crle)
      implicit double precision (a-h,o-z)
      logical irest
      common/rlenew/nrle0,nrle,irle
      common/newopt/nopt(6)
      common/rstrt/irest,ccopt(10),t32
      common/nt3t3/nto3,ntt3,lt3
      common/rlemat/xr(1000)
      common/itrat/iter,ix1,ix2
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension ti(1),t2(1),o2(1),crle(1)
      if(iter.eq.1.or.irest)then
         if(irest)irest=.false.
         call zeroma(xr,1,1000)
         nrle0=5
         nrle=1
         irle=1
      endif
      nl=it3(no,no,no-1)
      nlt4=it4(no,no,no-1,no-1)
      lent4=nlt4*nu4
      lent3=nl*nu3
      len=nou+no2u2
      lentot=len
      if(nopt(1).ge.2)lentot=len+lent3
      io=2*irle-1
      id=2*irle
      call ro1(nou,1,o2)
      call ro2hpp(0,no,nu,ti,o2(nou+1))
      call wrrle(11,io,len,o2)
      call ro1(nou,2,t2)
      call rt2hpp(no,nu,ti,t2(nou+1))
      call vecsub(t2,o2,len)
      call wrrle(11,id,len,t2)
      ioff=1
      iot3=(irle-1)*nl
      if(nopt(1).ge.2)then
         do 50 i=1,no
            do 49 j=1,i
               j1=j
               if(i.eq.j)j1=j-1
               do 48 k=1,j1
                  kkk=it3(i,j,k)
                  if(iter.eq.1)then
                     call zeroma(o2,1,nu3)
                  else
                     call rdvt3o(kkk,nu,o2)
                  endif
                  call rdvt3n(kkk,nu,t2)
                  ir=iot3+kkk
                  call wrrle(12,ir,nu3,o2)
                  call vecsub(t2,o2,nu3)
                  call wrrle(13,ir,nu3,t2)
 48            continue
 49         continue
 50      continue
      endif
      do 10 i=1,nrle
         id=2*i
         call rdrle(11,id,len,o2)
         if(nopt(1).ge.2)call rtot(13,no,nu,i,o2(len+1))
         do 11 j=1,i
            jd=2*j
            ij=i+(j-1)*nrle0
            ji=j+(i-1)*nrle0
            call rdrle(11,jd,len,t2)
            if(nopt(1).ge.2)call rtot(13,no,nu,j,t2(len+1))
            x=sdot(lentot,o2,1,t2,1)      
            xr(ji)=x
            xr(ij)=x
 11      continue
 10   continue
      call dodiis(xr,o2,crle,nrle0,nrle)
      call zeroma (o2,1,lentot)
      do 20 i=1,nrle
      id=2*i
      fx=crle(i)
      call rdrle(11,id,len,t2)
      if(nopt(1).ge.2)call rtot(13,no,nu,i,t2(len+1))
      call saxpy(lentot,fx,t2,1,o2,1)
 20   continue
      if (nrle.eq.1)goto 90
      do 30 i=1,nrle
      fx=crle(i)
      io=2*i-1
      call rdrle(11,io,len,t2)
      if(nopt(1).ge.2)call rtot(12,no,nu,i,t2(len+1))
      call saxpy(lentot,fx,t2,1,o2,1)
 30   continue
      call wo1(nou,2,o2) 
      call wo2hpp(no,nu,ti,o2(nou+1))
      if(nopt(1).ge.2)call wtot(ntt3,no,nu,1,o2(len+1))
 90   continue
      nrle=nrle+1
      if(nrle.gt.nrle0)nrle=nrle0
      irle =irle+1
      if(irle.gt.nrle)irle=1
 1009 format(5f15.10)
      return
      end
      subroutine stot4(irle,no,nu,o4,t4)
      implicit double precision (a-h,o-z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/itrat/iter,ix1,ix2
      dimension o4(1),t4(1)
      nlt4=it4(no,no,no-1,no-1)
      iot4=(irle-1)*nlt4
      kkkk=0
      do 50 i=1,no
         do 49 j=1,i
            j1=j
            if(i.eq.j)j1=j-1
            do 48 k=1,j1
               k1=k
               if(k.eq.j)k1=k-1
               do 47 l=1,k1
                  kkkk=kkkk+1
                  kkk=it4(i,j,k,l)
                  write(6,*)'irle,i,j,k,l,kk,it4:',irle,i,j,k,l,
     *kkkk,kkk
                  if(iter.eq.1)then
                     call zeroma(o4,1,nu4)
                  else
                     call reado4(kkk,nu,o4)
                  endif
                  call readt4(kkk,nu,t4)
                  ir=iot4+kkk
                  write(6,*)'iot4,kkk,ir:',iot4,kkk,ir
                  call wrrle(14,ir,nu4,o4)
                  call vecsub(t4,o4,nu4)
                  call wrrle(15,ir,nu4,t4)
 47            continue
 48         continue
 49      continue
 50   continue
      return
      end
      subroutine reado4(ir,nu,o4)
      implicit double precision (a-h,o-z)
      dimension o4(nu,nu,nu,nu)
      common/newt4/nt4,no4,lt4,nall4,ll4
      read(no4,rec=ir)o4
      return
      end
      subroutine readt4(ir,nu,t4)
      implicit double precision (a-h,o-z)
      dimension t4(nu,nu,nu,nu)
      common/newt4/nt4,no4,lt4,nall4,ll4
      read(nt4,rec=ir)t4
      return
      end
      subroutine wtott4(no,nu,irle,t4)
      implicit double precision (a-h,o-z)
      COMMON/NN/NO2,NO3,NnO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/newt4/nt4,no4,lt4,nall4,ll4
      dimension t4(1)
      ioff=1
      do 10 i=1,no
      do 10 j=1,i
      do 10 k=1,j
      if(i.eq.k)goto 10
      do 11 l=1,k
      if (j.eq.l)goto11
      kkkk=it4(i,j,k,l)
      call wrrle(nt4,kkkk,nu4,t4(ioff))
      ioff=ioff+nu4
 11   continue
 10   continue
      return
      end
      subroutine rdrle(nf,ir,n,a)
      implicit double precision (a-h,o-z)
      dimension a(n)
      read(nf,rec=ir)a
      return
      end
      SUBROUTINE RT2HPP(NO,NU,TI,T2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER A,B
      COMMON/NEWCCSD/NTT2
      DIMENSION TI(NU,NU,NO),T2(NO,NU,NU,NO)
      NOU2=NO*NU*NU
      NLAST=NO+2
      DO 1 I=1,NO
      IASV=NLAST+I
      READ(NTT2,REC=IASV)TI
      DO 2 J=1,NO
      DO 2 A=1,NU
      DO 2 B=1,NU
      T2(I,A,B,J)=TI(A,B,J)
 2    CONTINUE
 1    CONTINUE
      RETURN
      END
      subroutine rtot(nfil,no,nu,irle,t)
      implicit double precision (a-h,o-z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension t(1)
      nl=it3(no,no,no-1)
      ioff=1
      io=(irle-1)*nl
      do 10 i=1,no
      do 10 j=1,i
      j1=j
      if(i.eq.j)j1=j-1
      do 10 k=1,j1
      kkk=it3(i,j,k)
      ir=io+kkk
      call rdrle(nfil,ir,nu3,t(ioff))
      ioff=ioff+nu3
 10   continue
      return
      end
      subroutine rtott4(nfil,no,nu,irle,t)
      implicit double precision (a-h,o-z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension t(1)
      nl=it4(no,no,no-1,no-1)
      ioff=1
      io=(irle-1)*nl
      do 10 i=1,no
      do 10 j=1,i
      do 10 k=1,j
      do 10 l=1,k
      if(i.eq.k.or.j.eq.l)goto 10
      kkk=it4(i,j,k,l)
      ir=io+kkk
      call rdrle(nfil,ir,nu4,t(ioff))
      ioff=ioff+nu4
 10   continue
      return
      end

      SUBROUTINE WO2new(NO,NU,TI,T2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER A,B
      COMMON /NEWT4/NT4,NTO4,LT4,NALL4,LL4
      DIMENSION TI(NU,NU,NO),T2(NO,NU,NU,NO)
      DO 1 I=1,NO
      DO 2 J=1,NO
      DO 2 A=1,NU
      DO 2 B=1,NU
      TI(A,B,J)=T2(I,A,B,J)
 2    CONTINUE
      IASV=I
      WRITE(nall4,REC=IASV)TI
 1    CONTINUE
      RETURN
      END
      subroutine wrrle(nf,ir,n,a)
      implicit double precision (a-h,o-z)
      dimension a(n)
      write(nf,rec=ir)a
      return
      end

      subroutine wtot(nfil,no,nu,irle,t)
      implicit double precision (a-h,o-z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension t(1)
      nl=it3(no,no,no-1)
      ioff=1
      io=(irle-1)*nl
      do 10 i=1,no
      do 10 j=1,i
      j1=j
      if(i.eq.j)j1=j-1
      do 10 k=1,j1
      kkk=it3(i,j,k)
      ir=io+kkk
      call wrrle(nfil,ir,nu3,t(ioff))
      ioff=ioff+nu3
 10   continue
      return
      end
