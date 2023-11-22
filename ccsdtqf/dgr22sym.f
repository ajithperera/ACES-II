      subroutine dgr22sym(no,nu,o2,a,oeh,oep)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical smlv
      integer dirprd,dissyv,dissyt,dis233
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/file49/ibg(8),lrc
      common/timeinfo/timein,timenow,timetot,timenew
      common/sympop/irpdpd(8,22),isytyp(2,500),njunk(18)
      DIMENSION O2(nu,no,no,nu),a(1),oeh(no),oep(nu)
      common/syminf/nstart,nirrep,irrepa(255),irrepb(255),dirprd(8,8)
      common/flags/iflags(100)
      common/avlspace/maxcor
      equivalence(icllvl,iflags(2))
      DATA ZERO/0.0D+0/,ONE/1.0D+0/TWO/2.0D+0/,half/0.5d+0/
c      call timer(tcpu,tsys,0)
      smlv=icllvl.eq.11.or.icllvl.eq.13.or.icllvl.eq.14.or.
c     *icllvl.eq.15.or.icllvl.eq.16.or.icllvl.eq.17.or.icllvl.eq.19.or.
     *icllvl.eq.15.or.icllvl.eq.19.or.
     *icllvl.eq.21.or.icllvl.eq.22
c      call timer(1)
      nr1=1+nu2
      nr2=nr1+nu2
      nr3=nr2+nu2
      nr4=nr3+nu2
      kkk=1
      do 1 irp=1,nirrep
      dissyv=irpdpd(irp,isytyp(1,233))
      numsyv=irpdpd(irp,isytyp(2,233))
      if(smlv)dissyv=numsyv
      ibg(irp)=kkk
      nrc=dissyv/256+1
      kkk=kkk+nrc*numsyv
 1    continue
c
      n1=no2u2+1
      call ro2hpp(1,no,nu,a(n1),a)
      call vecmul(a,no2u2,half)
      call t2den(a,a(n1),oeh,oep,no,nu)
      call insitu(no,nu,nu,no,a(n1),a,13)
      call tranmd(a,nu,nu,no,no,12)
      call ro2hpp(0,no,nu,a(n1),o2)
      call transq(o2,nou)
      call insitu(nu,no,no,nu,a(n1),o2,13)
c
      n2=n1+nu3
      n3=n2+nu2
      do 101 i=1,nu
         call matmulsk(a,o2(1,1,1,i),a(n1),nu2,nu,no2,1,0)
         call trant3(a(n1),nu,1)
         call symfa(37,nu,i,a(n3),a(n2),a(n1))
         call rpakv(35,nu,i,a(n2),a(n3),a(n1))
         call symfa(38,nu,i,a(n3),a(n2),a(n1))
 101  continue
      call dr22a(nirrep,a)
      call rdsymo4pak(nu,a,a(nr1),a(nr2),a(nr3),a(nr4),39,36)
      n1=nu3+1
      n2=n1+no2u2
      n3=n2+nu2
      call insitu(no,no,nu,nu,a(n1),o2,13)
      call zeroma(a(n1),1,no2u2)
      do 11 ia=1,nu
         call rpakv(36,nu,ia,a(n2),a(n3),a)
         call trant3(a,nu,5)
         call matmulsk(a,o2(1,1,1,ia),a(n1),nu2,no2,nu,0,0)
 11   continue
      call insitu(nu,nu,no,no,a(n2),a(n1),13)
      call symetr(a(n1),no,nu)
      call engy1(no,nu,a,a(n1),a(n2))
      call rdsymo4pak(nu,a,a(nr1),a(nr2),a(nr3),a(nr4),37,36)
      n2=1+nu3
      n3=n2+nu2
      do 12 ia=1,nu
         call rpakv(36,nu,ia,a(n2),a(n3),a)
         call trant3(a,nu,3)
         call symfa(37,nu,ia,a(n3),a(n2),a)
 12   continue
      call dr22b(nirrep,a)
      do 13 ia=1,nu
         call rpakv(36,nu,ia,a(n2),a(n3),a)
         call trant3(a,nu,3)
         call symt21(a,nu,nu,nu,1,13)
         call symfa(37,nu,ia,a(n3),a(n2),a)
 13   continue
      n1=nu3+1
      n2=n1+nu2
      do 123 ia=1,nu
         kk=(ia-1)*nu3
         call rpakv(35,nu,ia,a(n1),a(n2),a)
         call trant3(a,nu,4)
         call symfa(38,nu,ia,a(n2),a(n1),a)
 123  continue
      call dr22c(nirrep,a)
      close(unit=37,status='delete')
      close(unit=38,status='delete')
      call tranmd(o2,nu,no,no,nu,14)
      call rdsymo4pak(nu,a,a(nr1),a(nr2),a(nr3),a(nr4),39,36)
      n1=nu3+1
      n2=n1+no2u2
      n3=n2+nu2
      call zeroma(a(n1),1,no2u2)
      do 10 ia=1,nu
         call rpakv(36,nu,ia,a(n2),a(n3),a)
         call trant3(a,nu,5)
         call matmulsk(a,o2(1,1,1,ia),a(n1),nu2,no2,nu,0,0)
 10   continue
      call insitu(nu,nu,no,no,a(n2),a(n1),13)
      call symetr(a(n1),no,nu)
      call engy1(no,nu,a,a(n1),a(n2))
      close(unit=39,status='delete')
      close(unit=36,status='delete')
c      call timer(tcpu,tsys,0)
c      write(6,2348)tcpu
c      call timer(1)
c      write(6,2348)timenew
c 2348 format('Diagram nr.22 required', f15.2,'  seconds')
      return
      end
      subroutine fulmat(irp,ndis,npart,nsiz,ntot,nt1,nt2,a)
      implicit double precision (a-h,o-z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension a(1)
      data two/2.0d+0/
      call getlst   (a(nt1),1,npart,1,irp,233)
      call trsym12b (a(nt1),irp,1,npart,nsiz)
      call rdlistnew(37,irp,1,ndis,a,nsiz)
      call tmatmul(a,a(nt1),a(nt2),ndis,npart,nsiz,1,0)
      call wrlistnew(39,irp,1,ndis,a(nt2),nsiz)
c
      call rdlistnew(36,irp,1,ndis,a,nsiz)
      call tmatmul(a,a(nt1),a(nt2),ndis,ndis,nsiz,1,0)
      call rdlistnew(37,irp,1,ndis,a(nt1),nsiz)
      call vecsuba(a,a(nt1),ntot,two)
      call rdlistnew(38,irp,1,ndis,a(nt1),nsiz)
      call tmatmul(a,a(nt1),a(nt2),ndis,ndis,nsiz,0,1)
      call wrlistnew(40,irp,1,ndis,a(nt2),nsiz)
      return
      end
      subroutine fulmat1(irp,ndis,npart,nsiz,ntot,nt1,nt2,a)
      implicit double precision (a-h,o-z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension a(1)
      data two/2.0d+0/
      call rdlistnew(38,irp,1,npart,a(nt1),nsiz)
      call rdlistnew(37,irp,1,ndis,a,nsiz)
      call tmatmul(a,a(nt1),a(nt2),ndis,npart,nsiz,1,0)
      call wrlistnew(39,irp,1,ndis,a(nt2),nsiz)
      return
      end
      subroutine fulmat2(irp,ndis,npart,nsiz,ntot,nt1,nt2,a)
      implicit double precision (a-h,o-z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension a(1)
      data two/2.0d+0/
      call rdlistnew(38,irp,1,npart,a(nt1),nsiz)
      call rdlistnew(37,irp,1,ndis,a,nsiz)
      call tmatmul(a,a(nt1),a(nt2),ndis,ndis,nsiz,1,0)
      call wrlistnew(39,irp,1,ndis,a(nt2),nsiz)
      return
      end
      subroutine fulmat3(irp,ndis,npart,nsiz,ntot,nt1,nt2,a)
      implicit double precision (a-h,o-z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension a(1)
      data two/2.0d+0/
      call rdlistnew(37,irp,1,ndis,a,nsiz)
      call rdlistnew(38,irp,1,ndis,a(nt1),nsiz)
      call rdlistnew(39,irp,1,ndis,a(nt2),nsiz)
      call tmatmul(a,a(nt1),a(nt2),ndis,ndis,nsiz,0,1)
      call wrlistnew(39,irp,1,ndis,a(nt2),nsiz)
      return
      end
      subroutine partmat1(irp,ndis,nf1,npart,nsiz,ntot,nt1,nt2,a)
      implicit double precision (a-h,o-z)
      dimension a(1)
      data two/2.0d+0/
      call rdlistnew(38,irp,nf1,npart,a(nt1),nsiz)
      call tmatmul(a,a(nt1),a(nt2),ndis,npart,nsiz,1,0)
      call wrlistnew(39,irp,nf1,npart,a(nt2),nsiz)
      return
      end
      subroutine partmat2(irp,ndis,nf1,npart,nsiz,ntot,nt1,nt2,a)
      implicit double precision (a-h,o-z)
      dimension a(1)
      data two/2.0d+0/
      call rdlistnew(38,irp,nf1,npart,a(nt1),nsiz)
      call tmatmul(a,a(nt1),a(nt2),ndis,npart,nsiz,1,0)
      call wrlistnew(39,irp,nf1,npart,a(nt2),nsiz) 
      return
      end
      subroutine partmat3(irp,ndis,nf1,npart,nsiz,ntot,nt1,nt2,a)
      implicit double precision (a-h,o-z)
      dimension a(1)
      data two/2.0d+0/
      call rdlistnew(39,irp,nf1,npart,a(nt2),nsiz)
      call rdlistnew(38,irp,nf1,npart,a(nt1),nsiz)
      call tmatmul(a,a(nt1),a(nt2),ndis,npart,nsiz,0,1)
      call wrlistnew(39,irp,nf1,npart,a(nt2),nsiz)
      return
      end

      subroutine partmat(irp,ndis,npart,nsiz,ntot,nt1,nt2,a)
      implicit double precision (a-h,o-z)
      dimension a(1)
      data two/2.0d+0/
      call rdlistnew(37,irp,1,ndis,a,nsiz)
      do 100 i=1,ndis,npart
         call partmat1(irp,ndis,i,npart,nsiz,ntot,nt1,nt2,a)
 100  continue
      call rdlistnew(36,irp,1,ndis,a,nsiz)
      do 200 i=1,ndis
         call partmat2(irp,ndis,i,npart,nsiz,ntot,nt1,nt2,a)
 200  continue
      call rdlistnew(37,irp,1,ndis,a(nt1),nsiz)
      call vecsuba(a,a(nt1),ntot,two)
      do 300 i=1,ndis
         call partmat3(irp,ndis,i,npart,nsiz,ntot,nt1,nt2,a)
 300  continue
      return
      end
      subroutine partmata(irp,ndis,npart,nsiz,ntot,nt1,nt2,a)
      implicit double precision (a-h,o-z)
      dimension a(1)
      data two/2.0d+0/
      call rdlistnew(37,irp,1,ndis,a,nsiz)
      do 100 i=1,ndis,npart
         call partmat1(irp,ndis,i,npart,nsiz,ntot,nt1,nt2,a)
 100  continue
      return
      end
      subroutine partmatb(irp,ndis,npart,nsiz,ntot,nt1,nt2,a)
      implicit double precision (a-h,o-z)
      dimension a(1)
      data two/2.0d+0/
      call rdlistnew(37,irp,1,ndis,a,nsiz)
      do 200 i=1,ndis
         call partmat2(irp,ndis,i,npart,nsiz,ntot,nt1,nt2,a)
 200  continue
      return
      end
      subroutine partmatc(irp,ndis,npart,nsiz,ntot,nt1,nt2,a)
      implicit double precision (a-h,o-z)
      dimension a(1)
      data two/2.0d+0/
      call rdlistnew(37,irp,1,ndis,a,nsiz)
      do 300 i=1,ndis
         call partmat3(irp,ndis,i,npart,nsiz,ntot,nt1,nt2,a)
 300  continue
      return
      end
      subroutine dr22a(nirrep,a)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer dis233
      logical smlv
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/sympop/irpdpd(8,22),isytyp(2,500),njunk(18)
      common/avlspace/maxcor
      common/flags/iflags(100)
      dimension a(1)
      equivalence(icllvl,iflags(2))
      smlv=icllvl.eq.11.or.icllvl.eq.13.or.icllvl.eq.14.or.
c     *icllvl.eq.15.or.icllvl.eq.16.or.icllvl.eq.17.or.icllvl.eq.19.or.
     *icllvl.eq.15.or.icllvl.eq.19.or.
     *icllvl.eq.21.or.icllvl.eq.22
      do 1000 irp=1,nirrep
         num233=irpdpd(irp,isytyp(2,233))
         dis233=irpdpd(irp,isytyp(1,233))
         if(smlv)dis233=num233
         ntotsiz=num233*dis233
         if(maxcor.gt.3*ntotsiz)then
cc            write(6,*)'fulmat1 selected'
            nt1=ntotsiz+1
            nt2=ntotsiz+nt1
            call fulmat1(irp,num233,num233,dis233,ntotsiz,nt1,nt2,a)
         else
            if (maxcor.gt.ntotsiz+2*no2u2)then
cc               write(6,*)'partmata selected'
               nt1=ntotsiz+1
               nt2=dis233+nt1
               call partmata(irp,num233,1,dis233,ntotsiz,nt1,nt2,a)
            endif
         endif
 1000 continue
      return
      end
      subroutine dr22b(nirrep,a)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer dis233
      logical smlv
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/sympop/irpdpd(8,22),isytyp(2,500),njunk(18)
      common/avlspace/maxcor
      common/flags/iflags(100)
      dimension a(1)
      equivalence(icllvl,iflags(2))
      smlv=icllvl.eq.11.or.icllvl.eq.13.or.icllvl.eq.14.or.
c     *icllvl.eq.15.or.icllvl.eq.16.or.icllvl.eq.17.or.icllvl.eq.19.or.
     *icllvl.eq.15.or.icllvl.eq.19.or.
     *icllvl.eq.21.or.icllvl.eq.22
      do 1000 irp=1,nirrep
         num233=irpdpd(irp,isytyp(2,233))
         dis233=irpdpd(irp,isytyp(1,233))
         if(smlv)dis233=num233
         ntotsiz=num233*dis233
         if(maxcor.gt.3*ntotsiz)then
ccc            write(6,*)'fulmat2 selected'
            nt1=ntotsiz+1
            nt2=ntotsiz+nt1
            call fulmat2(irp,num233,num233,dis233,ntotsiz,nt1,nt2,a)
         else
            if (maxcor.gt.ntotsiz+2*no2u2)then
ccc               write(6,*)'partmatb selected'
               nt1=ntotsiz+1
               nt2=dis233+nt1
               call partmatb(irp,num233,1,dis233,ntotsiz,nt1,nt2,a)
            endif
         endif
 1000 continue
      return
      end
      subroutine dr22c(nirrep,a)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer dis233
      logical smlv
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/sympop/irpdpd(8,22),isytyp(2,500),njunk(18)
      common/avlspace/maxcor
      common/flags/iflags(100)
      dimension a(1)
      equivalence(icllvl,iflags(2))
      smlv=icllvl.eq.11.or.icllvl.eq.13.or.icllvl.eq.14.or.
c     *icllvl.eq.15.or.icllvl.eq.16.or.icllvl.eq.17.or.icllvl.eq.19.or.
     *icllvl.eq.15.or.icllvl.eq.19.or.
     *icllvl.eq.21.or.icllvl.eq.22
      do 1000 irp=1,nirrep
         num233=irpdpd(irp,isytyp(2,233))
         dis233=irpdpd(irp,isytyp(1,233))
         if(smlv)dis233=num233
ccc         write(6,*)'num233,dis233:',num233,dis233
         ntotsiz=num233*dis233
         if(maxcor.gt.3*ntotsiz)then
ccc            write(6,*)'fulmat3 selected'
            nt1=ntotsiz+1
            nt2=ntotsiz+nt1
            call fulmat3(irp,num233,num233,dis233,ntotsiz,nt1,nt2,a)
         else
            if (maxcor.gt.ntotsiz+2*no2u2)then
cc               write(6,*)'partmatc selected'
               nt1=ntotsiz+1
               nt2=dis233+nt1
               call partmatc(irp,num233,1,dis233,ntotsiz,nt1,nt2,a)
            endif
         endif
 1000 continue
      return
      end





