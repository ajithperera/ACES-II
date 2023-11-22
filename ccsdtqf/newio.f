      subroutine rdlist(nfil,irec,a,na,nrec,nres)
      implicit double precision(a-h,o-z)
      common/blokre/blk(256)
      dimension a(1)
      ioff=1
      iirec=irec
      if(nfil.eq.40)then
      endif
      do 10 i=1,nrec-1
         call rrdlist(nfil,iirec,blk,na)
         call veccop(na,a(ioff),blk)
         iirec=irec+i
         ioff=ioff+na
 10   continue
         call rrdlist(nfil,iirec,blk,na)
      call veccop(nres,a(ioff),blk)
      return
      end
      subroutine rdlistnew(nfile,irp,nfst,ndis,a,sdis)
      implicit double precision (a-h,o-z)
      integer sdis
      common/file49/ibg(8),lrec
      dimension a(1)
      ircc=ibg(irp)
      nrec=sdis/lrec+1
      ntot=nfst+ndis-1
      nres=mod(sdis,lrec)
      iof=1
      do 10 idis=nfst,ntot
         irec=ircc+(idis-1)*nrec
         call rdlist(nfile,irec,a(iof),lrec,nrec,nres)
         iof=iof+sdis
 10   continue
      return
      end
      subroutine rrdlist(nfil,irec,a,na)
      implicit double precision(a-h,o-z)
      dimension a(na)
      if(nfil.eq.40)then
c      write(6,*)'nfil,irec,na:',nfil,irec,na
      endif
      read(nfil,rec=irec)a
      return
      end
      subroutine setsym
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer dissyt,dissyv,dirprd,pop,vrt,disnum,a,b,c,d,dis233
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/sym/pop(8,2),vrt(8,2),ntaaa(6)
      common/sympop/irpdpd(8,22),isytyp(2,500),njunk(18)
      common/prevnu/idn(8,8)
      common/syminf/nstart,nirrep,irrepa(255),irrepb(255),dirprd(8,8)
      kk=0
      do 1 i=1,nirrep
         do 2 j=1,vrt(i,1)
            kk=kk+1
 2       continue
 1    continue
      ir1=0
      ir2=0
      ir3=0
      ir4=0
      ir5=0
      ir6=0
      ir7=0
      ir8=0
      do 3 i=1,nirrep
         do 4 j=1,nirrep
         irr=dirprd(i,j)
         goto (11,12,13,14,15,16,17,18)irr
 11      continue
         idn(i,j)=ir1
         ir1=ir1+vrt(i,1)*vrt(j,1)
         goto 4
 12      continue
         idn(i,j)=ir2
         ir2=ir2+vrt(i,1)*vrt(j,1)
         goto 4
 13      continue
         idn(i,j)=ir3
         ir3=ir3+vrt(i,1)*vrt(j,1)
         goto 4
 14      continue
         idn(i,j)=ir4
         ir4=ir4+vrt(i,1)*vrt(j,1)
         goto 4
 15      continue
         idn(i,j)=ir5
         ir5=ir5+vrt(i,1)*vrt(j,1)
         goto 4
 16      continue
         idn(i,j)=ir6
         ir6=ir6+vrt(i,1)*vrt(j,1)
         goto 4
 17      continue
         idn(i,j)=ir7
         ir7=ir7+vrt(i,1)*vrt(j,1)
         goto 4
 18      continue
         idn(i,j)=ir8
         ir8=ir8+vrt(i,1)*vrt(j,1)
         goto 4
 4    continue
 3    continue
      do 8 i=1,nirrep
         do 8 j=1,nirrep
c            write(6,234)i,j,dirprd(i,j),idn(i,j)
 8       continue
 234     format('i,j,prod,idn(i,j):',4i4)
c         stop
         return
         end
      subroutine symfa(nfile,nu,a,isym,ti,v)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer dissyt,dissyv,dirprd,pop,vrt,disnum,a,b,c,d,dis233
      logical smlv
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION TI(1),v(nu,nu,nu),isym(nu,2)
      common/sym/pop(8,2),vrt(8,2),ntaaa(6)
      common/prevnu/idn(8,8)
      common/sympop/irpdpd(8,22),isytyp(2,500),njunk(18)
      common/avlspace/maxcor
      common/file49/ibg(8),lrec
      common/syminf/nstart,nirrep,irrepa(255),irrepb(255),dirprd(8,8)
      common/newio/npp1,nhh,npp,ntvoe,nto3,no1,nu1
      common/flags/iflags(100)
      equivalence(icllvl,iflags(2))
      call setsym
      smlv=icllvl.eq.11.or.icllvl.eq.13.or.icllvl.eq.14.or.
     *icllvl.eq.15.or.icllvl.eq.19.or.
     *icllvl.eq.21.or.icllvl.eq.22
      kk=0
      do 1 i=1,nirrep
         do 2 j=1,vrt(i,1)
            kk=kk+1
            isym(kk,1)=i      
            isym(kk,2)=j      
 2       continue
 1    continue
      lrec=256
      ira=isym(a,1)
      ika=isym(a,2)
      do 3 b=1,nu
         irb=isym(b,1)                  
         ikb=isym(b,2)                  
         irrepab=dirprd(ira,irb)          
         dis233=irpdpd(irrepab,isytyp(1,233))
         num233=irpdpd(irrepab,isytyp(2,233))
         if(smlv)dis233=num233
         disnum=idn(ira,irb)+(ika-1)*vrt(irb,1)+ikb  
         isymcd=0
         do 4 c=1,nu
            irc=isym(c,1)
            ikc=isym(c,2)
            nvrtc=vrt(irc,1)
            do 5 d=1,nu
               ird=isym(d,1)
               ikd=isym(d,2)
               if(dirprd(irc,ird).eq.irrepab)then
                  isymcd=isymcd+1
                  ti(isymcd)=v(c,d,b)
               endif
 5          continue
 4       continue
         if(nfile.eq.38)then
         endif
         call wrlistnew(nfile,irrepab,disnum,1,ti,dis233)
 3    continue
      return
      end
      subroutine trsym12b(v,irp,nfirst,ndis,nsiz)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer dissyt,dissyv,dirprd,pop,vrt,disnum,a,b,c,d,dis233
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION v(1)
      common/sym/pop(8,2),vrt(8,2),ntaaa(6)
      common/prevnu/idn(8,8)
      common/sympop/irpdpd(8,22),isytyp(2,500),njunk(18)
      common/avlspace/maxcor
      common/file49/ibg(8),lrec
      common/syminf/nstart,nirrep,irrepa(255),irrepb(255),dirprd(8,8)
      common/newio/npp1,nhh,npp,ntvoe,nto3,no1,nu1
      do 10 i=nfirst,ndis
         iofdis=(i-1)*nsiz
         iof=0
         do 1 jrp=1,nirrep
            krp=dirprd(irp,jrp)
            if(jrp.eq.krp)then
               nvrt=vrt(jrp,1)
               iof=idn(jrp,krp)
               do 11 ii=1,nvrt
                  do 12 jj=1,ii-1
                     iijj=(ii-1)*nvrt+jj+iof+iofdis
                     jjii=(jj-1)*nvrt+ii+iof+iofdis
                     x=v(iijj)
                     v(iijj)=v(jjii)
                     v(jjii)=x
 12               continue
 11            continue
            endif
            if(jrp.lt.krp)then
               nvrtj=vrt(jrp,1)
               nvrtk=vrt(krp,1)
               iofjk=idn(jrp,krp)
               iofkj=idn(krp,jrp)
               do 13 jj=1,nvrtj
                  do 14 kk=1,nvrtk
                     jjkk=(jj-1)*nvrtk+kk+iofjk+iofdis
                     kkjj=(kk-1)*nvrtj+jj+iofkj+iofdis
                     x=v(jjkk)
                     v(jjkk)=v(kkjj)
                     v(kkjj)=x
 14               continue
 13            continue
            endif
 1       continue
 10   continue
 235     format(2f15.10)
      return
      end
      subroutine wrlist(nfil,irec,a,na,nrec)
      implicit double precision(a-h,o-z)
      dimension a(1)
      ioff=1
      do 10 i=1,nrec
         iirec=irec+i-1
         call wwrlist(nfil,iirec,a(ioff),na)
         ioff=ioff+na
 10   continue
      return
      end
      subroutine wrlistnew(nfile,irp,nfst,ndis,a,sdis)
      implicit double precision (a-h,o-z)
      integer sdis
      common/file49/ibg(8),lrec
      dimension a(1)
      ircc=ibg(irp)
      nrec=sdis/lrec+1
      ntot=nfst+ndis-1
      iof=1
      do 10 idis=nfst,ntot
         irec=ircc+(idis-1)*nrec
         call wrlist(nfile,irec,a(iof),lrec,nrec)
         iof=iof+sdis
 10   continue
      return
      end
      subroutine wwrlist(nfil,irec,a,na)
      implicit double precision(a-h,o-z)
      dimension a(na)
      if(nfil.eq.40)then
      endif
      write(nfil,rec=irec)a
      return
      end
