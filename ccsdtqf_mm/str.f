      SUBROUTINE RDVT3ONWx(I,J,K,NU,V,iv)
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      COMMON/newio/ni1,ni2,ni3,ni4,no3,no1,nu1
      DIMENSION V(NU,NU,NU)
      call ienter(41)
      nnrec=0
 489  format('i,j,k:',i3,5x,3i6,3x,3i3)
      IF (I.GT.J)THEN
      IF(J.GT.K) THEN
         IAS=IT3(I,J,K)+nnrec
      READ(NO3,REC=IAS)V
      ELSE
      IF(I.GT.K)THEN
      IAS=IT3(I,K,J)+nnrec
      READ(NO3,REC=IAS)V
      CALL TRANMD(V,NU,NU,NU,1,23)
      ELSE
      IAS=IT3(K,I,J)+nnrec
      READ(NO3,REC=IAS)V
      CALL TRANMD(V,NU,NU,NU,1,312)
      ENDIF
      ENDIF
      ELSE
      IF (I.GT.K)THEN
      IAS=IT3(J,I,K)+nnrec
      READ(NO3,REC=IAS)V
      CALL TRANMD(V,NU,NU,NU,1,12)
      ELSE
      IF(J.GT.K) THEN
      IAS=IT3(J,K,I)+nnrec
      READ(NO3,REC=IAS)V
      CALL TRANMD(V,NU,NU,NU,1,231)
      ELSE
      IAS=IT3(K,J,I)+nnrec
      READ(NO3,REC=IAS)V
      CALL TRANMD(V,NU,NU,NU,1,13)
      ENDIF
      ENDIF
      ENDIF
      call iexit(41)
      RETURN
      END
      subroutine rpakt3old(ipn,n,t3)
      implicit double precision (a-h,o-z)
      dimension t3(n)
      parameter(lrc=2728)
      common/bufor/buf(lrc),ibuf(lrc)
      common/wpakt3/nfr(2000),nent(2000)
      call zeroma(t3,1,n)
      irec=nfr(ipn)
 10   read(45,rec=irec)buf,ibuf,lrec
      if(lrec.eq.0)goto 12
      do 11 ir=1,lrec
         int3=ibuf(ir)
         t3(int3)=buf(ir)
 11   continue
      irec=irec+1
      goto 10
 12   continue
      return
      end
      subroutine rpakt3(ipn,n,t3)
      implicit double precision (a-h,o-z)
      dimension t3(n)
      parameter(lrc=2728)
      common/bufor/buf(lrc),ibuf(lrc)
      common/wpakt3/nfr(2000)
      call zeroma(t3,1,n)
      irec=nfr(ipn)
 10   continue
      read(45,rec=irec)buf,ibuf,lrec
c 10   read(45,rec=irec)buf,ibuf,lrec
      if(lrec.eq.0)goto 12
      do 11 ir=1,lrec
         int3=ibuf(ir)
         t3(int3)=buf(ir)
 11   continue
      irec=irec+1
      goto 10
 12   continue
      return
      end
      subroutine izo(ia,n)
      implicit double precision(a-h,o-z)
      dimension ia(n)
      do 10 i=1,n
         ia(i)=0
 10   continue
      return
      end
      subroutine vpakt3(ipn,ii,j,k,n,t3)
      implicit double precision(a-h,o-z)
      common/bufor/buf(2728),ibuf(2728)
      common/wpakt3/nfr(2000),nent(2000)
      common/actadres/iacta
      common/count/kkkk
      dimension t3(n)
      data tr/0.1d-12/
      call ienter(125)
      call zeroma(buf,1,2728)
      call izo(ibuf,2728)
      icount=0
      in=0
      icurrec=0
      if(nent(ipn).gt.1)then
         irec=nfr(ipn)
      else
         irec=iacta
         nfr(ipn)=iacta
      endif
      do 10 i =1,n
         x=t3(i)
         if (dabs(x).gt.tr)then
            icount=icount+1
            in=in+1
            buf(in)=x
            ibuf(in)=i
         endif
         if (in.eq.2728)then
            write(45,rec=irec)buf,ibuf,in
            in=0
            icurrec=icurrec+1
            irec=irec+1
         endif
 10   continue
      icurrec=icurrec+1
      write(45,rec=irec)buf,ibuf,in
      if(in.ne.0)then
         irec=irec+1
         icurrec=icurrec+1
         in=0
         write(45,rec=irec)buf,ibuf,in
      endif
      if(icount.gt.kkkk)then
         ikk=kkkk/2728
         icc=icount/2728
         write(6,*)'no of reserved records:',ikk,' actual:',icc
      endif
      if(nent(ipn).eq.1)then
         ikk=kkkk/2728+2
         if(mod(kkkk,2728).eq.0)ikk=ikk-1
         ico=ikk-icurrec
         iacta=irec+ico+3
      endif
      call iexit(125)
      return
      end
      subroutine vpakt3new(ipn,ii,j,k,n,t3)
      implicit double precision(a-h,o-z)
      common/bufor/buf(2728),ibuf(2728)
      common/wpakt3/nfr(2000)
      common/actadres/iacta
      dimension t3(n)
      data tr/0.1d-12/
      call ienter(125)
      call zeroma(buf,1,2728)
      call izo(ibuf,2728)
      icount=0
      in=0
      icurrec=0
         irec=iacta
         nfr(ipn)=iacta
      do 10 i =1,n
         x=t3(i)
         if (dabs(x).gt.tr)then
            icount=icount+1
            in=in+1
            buf(in)=x
            ibuf(in)=i
         endif
         if (in.eq.2728)then
            write(45,rec=irec)buf,ibuf,in
            in=0
            icurrec=icurrec+1
            irec=irec+1
         endif
 10   continue
      icurrec=icurrec+1
      write(45,rec=irec)buf,ibuf,in
      if(in.ne.0)then
         irec=irec+1
         icurrec=icurrec+1
         in=0
         write(45,rec=irec)buf,ibuf,in
      endif
      iacta=irec+1
      call iexit(125)
      return
      end
      subroutine vripak(nint,n,vec,ivec)
      implicit double precision (a-h,o-z)
      common/t3adr/nfr(2000)
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
      subroutine ww(ir,n,v,iv)
      implicit double precision (a-h,o-z)
      dimension v(n),iv(n)
      common/pak/nv4,intr
      common/lnt/lrecb,lrecw
      write(nv4,rec=ir)v,iv
      return
      end
      subroutine rr(ir,n,v,iv)
      implicit double precision (a-h,o-z)
      dimension v(n),iv(n)
      common/pak/nv4,intr
      common/lnt/lrecb,lrecw
      read(nv4,rec=ir)v,iv
      return
      end
