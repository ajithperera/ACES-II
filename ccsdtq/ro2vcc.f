      SUBROUTINE ro2hppnew(NO,NU,TI,T2,EH,EP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON /NEWT4/NT4,NTO4,LT4,NALL4,LL4
      common/pak/intg,intr
      common/wpak/nfr(30),nsz(30)
      DIMENSION TI(nu,nu,no),T2(no,nu,nu,no),EH(NO),EP(NU)
      call ro2hppsm(1,no,nu,ti,t2)
      do 10 i=1,no
         do 11 j=1,no
         do 11 a=1,nu
         do 11 b=1,nu
            den=eh(i)+eh(j)-ep(a)-ep(b)
            ti(a,b,j)=t2(j,b,a,i)/den
 11      continue
         irec=i
         write(nall4,rec=irec)ti
 10   continue
      call store(no2u2,non,t2,ti)
      nsz(7)=non
      call vripak(intg,7,non,t2,ti)
      return
      end
      SUBROUTINE rrghpp(NO,NU,TI,T2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON /NEWT4/NT4,NTO4,LT4,NALL4,LL4
      common/pak/intg,intr
      common/wpak/nfr(30),nsz(30)
      DIMENSION TI(nu,nu,no),T2(no,nu,nu,no)
      call ro2hppsm(1,no,nu,ti,t2)
      do 10 i=1,no
         do 11 j=1,no
         do 11 a=1,nu
         do 11 b=1,nu
            ti(a,b,j)=t2(j,b,a,i)
 11      continue
         irec=i+no
         write(nall4,rec=irec)ti
 10   continue
      call store(no2u2,non,t2,ti)
      nsz(3)=non
      call vripak(intg,3,non,t2,ti)
      return
      end
      SUBROUTINE rldhpp(NO,NU,TI,T2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON /NEWT4/NT4,NTO4,LT4,NALL4,LL4
      common/wpak/nfr(30),nsz(30)
      common/pak/intg,intr
      DIMENSION TI(nu,nu,no),T2(no,nu,nu,no)
      call ro2hppsm(2,no,nu,ti,t2)
      do 10 i=1,no
         do 11 j=1,no
         do 11 a=1,nu
         do 11 b=1,nu
            ti(a,b,j)=t2(j,b,a,i)
 11      continue
         irec=2*no+i
         write(nall4,rec=irec)ti
 10   continue
      call store(no2u2,non,t2,ti)
      nsz(4)=non
      call vripak(intg,4,non,t2,ti)
      return
      end
      SUBROUTINE rhhhh(no,nu,ti,v)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON /NEWT4/NT4,NTO4,LT4,NALL4,LL4
      common/wpak/nfr(30),nsz(30)
      common/pak/intg,intr
      DIMENSION TI(no,no,no),v(no,no,no,no)
      call rdhhhhsm(no,ti,v)
      nlast=4*no+nu
      do 10 i=1,no
         irec=nlast+i
         call wr(nall4,irec,no3,v(1,1,1,i))
 10   continue
      call store(no4,non,v,ti)
      nsz(2)=non
      call vripak(intg,2,non,v,ti)
      return
      end
      subroutine wr(nf,ir,n,v)
      implicit double precision (a-h,o-z)
      dimension v(n)
      write(nf,rec=ir)v
      return
      end
      subroutine rd(nf,ir,n,v)
      implicit double precision (a-h,o-z)
      dimension v(n)
      read(nf,rec=ir)v
      return
      end
      SUBROUTINE rhhhp(no,nu,ti,v)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON /NEWT4/NT4,NTO4,LT4,NALL4,LL4
      common/wpak/nfr(30),nsz(30)
      common/pak/intg,intr
      DIMENSION ti(1),v(no,no,no,nu)
      call rdvemsm(1,nu,no,ti,v)
      nlast=4*no
      do 10 a=1,nu
         irec=nlast+a
         call wr(nall4,irec,no3,v(1,1,1,a))
 10   continue
      call store(no3u,non,v,ti)
      nsz(5)=non
      call vripak(intg,5,non,v,ti)
      return
      end
      SUBROUTINE rhppp(no,nu,ti,v)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON /NEWT4/NT4,NTO4,LT4,NALL4,LL4
      common/wpak/nfr(30),nsz(30)
      common/pak/intg,intr
      DIMENSION ti(1),v(nu,nu,nu,no)
      call rdvemsm(0,no,nu,ti,v)
      nlast=3*no
      do 10 i=1,no
         irec=nlast+i
         call wr(nall4,irec,nu3,v(1,1,1,i))
 10   continue
      call tranmd(v,nu,nu,nu,no,23)
 99   format(f12.8)
      call storei(0,nu3,no,non,v,ti)
      nsz(6)=non
      call vripak(intg,6,non,v,ti)
      return
      end
      subroutine rpppp(no,nu,ti,v,ifr,inr)
      implicit double precision(a-h,o-z)
      integer a,b
      dimension ti(1),v(nu,nu),ifr(1),inr(1)
      common/newt4/nt4,not4,lt4,nall4,ll4
      common/newio/nvv,nhh,npp,nvo,nt3,no1,nu1
      common/pak/intg,intr
      common/wpak/nfr(30),nsz(30)
      data tr/0.1d-12/
      nu2=nu*nu
      nu3=nu2*nu
      call rifr(34,ifr,inr)
      do 10 i=1,nu
      do 10 j=1,nu
c         ij=(i-1)*nu+j
         ij=(j-1)*nu+i
         call readnu2(i,j,nu,v,ifr,inr)
         write(nvv,rec=ij)v
 10   continue
      kk=0
      ijkl=0
      in=0
      nlast=5*no+nu
      do 20 i=1,nu
         noff=1
      do 21 j=1,nu
         kk=kk+1
         call rrv4(kk,nu2,ti(noff))
         noff=noff+nu2
 21   continue
      irec=nlast+i
      call wr(nall4,irec,nu3,ti)
 20   continue
      return
      end
      subroutine ro2hppsm(io,no,nu,ti,t2)
      implicit double precision(a-h,o-z)
      character*8 lb1,lb2
      integer dirprd
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      dimension ti(1),t2(nu,nu,no,no)
      if (io.eq.0) then
         if(nirrep.gt.1)then
            call zeroma(t2,1,no2u2)
            lb1='SVAVB2  '
            lb2='SOAOB2  '
            i1=1+nu2
            i2=i1+nu2
            itot=i2+no2
         call grdsym(nu,nu,no,no,nu2,no2,t2,ti,ti(i1),ti(i2),lb1,lb2,46)
            call insitu(nu,nu,no,no,ti,t2,13)
            call tranmd(t2,no,nu,nu,no,23)
         else
            CALL GETLST(t2,1,NO2,1,1,63)
            call insitu(nu,nu,no,no,ti,t2,13)
            call tranmd(t2,no,nu,nu,no,23)
         endif
      endif
      if (io.eq.1) then
         if(nirrep.gt.1)then
            call zeroma(t2,1,no2u2)
            lb1='SVAVB2  '
            lb2='SOAOB2  '
            i1=1+nu2
            i2=i1+nu2
            itot=i2+no2
        call grdsym(nu,nu,no,no,nu2,no2,t2,ti,ti(i1),ti(i2),lb1,lb2,16)
            call insitu(nu,nu,no,no,ti,t2,13)
            call tranmd(t2,no,nu,nu,no,23)
         else
            CALL GETLST(T2,1,no2,1,1,16)
            call insitu(nu,nu,no,no,ti,t2,13)
            call tranmd(t2,no,nu,nu,no,23)
         endif
      endif
      if (io.eq.2) then
      call zeroma(t2,1,no2u2)
         if(nirrep.gt.1)then
            call zeroma(t2,1,no2u2)
            lb1='SVBOA2  '
            i1=1+nou
            i2=i1+nou
            itot=i2+nou
         call grdsym(nu,no,nu,no,nou,nou,t2,ti,ti(i1),ti(i2),lb1,lb1,25)
            call insitu(nu,no,nu,no,ti,t2,12)
         else
            CALL GETLST(T2,1,nou,1,1,25)
            call insitu(nu,no,nu,no,ti,t2,12)
         endif
      endif
      return
      end
      SUBROUTINE rdsymv(nu,t2,t2s,SYVEC1,syvec2,ifr,inr,nlist)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SYVEC1(1),syvec2(1),t2(nu,nu),t2s(1),ifr(1),inr(1)
      integer syvec1,dirprd,a,b,c,d,pop,vrt,spop,svrt,syvec2
      character*8 lb
      common/countv/irec
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),N(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      common/sym/pop(8,2),vrt(8,2),ntaa(6)
      common/sumpop/spop(8),svrt(8)
      common/flags/iflags(100)
      equivalence(icllvl,iflags(2))
      lb='SVAVB2  '
      no=nocco(1)
      nu=nvrto(1)
      nu2=nu*nu
      call getrec(20,'JOBARC',lb,nu2,syvec1)
      nfl=34
         call veccop(nu2,syvec2,syvec1)
      kkk=0
      nrec=nu2/1024+1
      irec=2*nrec
      DO 10 IRREP=1,NIRREP
         NSYDIS=IRPDPD(IRREP,ISYTYP(2,nLIST))
         NSYDSZ=IRPDPD(IRREP,ISYTYP(1,nLIST))
         NSIZ=NSIZ+NSYDIS*NSYDSZ
         do 9 inum=1,nsydis
            call zeroma(t2,1,nu2)
            kij=kij+1
            call getlst(t2s,inum,1,1,irrep,nlist)
            do 8 is=1,nsydsz
               kkk=kkk+1
               call pckindsk(kkk,3,syvec2,syvec1,nu2,nu2,a,b,c,d,nlist)
               t2(c,d)=t2s(is)
 8          continue
            call vpakv(nfl,nu,a,b,ifr,inr,t2)
 9       CONTINUE
 10   CONTINUE
      RETURN
      END
      subroutine rdvem4sm(io,no,nu,ti,v)
      implicit double precision(a-h,o-z)
      integer dirprd
      character*8 lb1,lb2
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      dimension ti(1),v(1)
      nou=no*nu
      if(io.eq.0)then
         if(nirrep.gt.1)then
            call zeroma(v,1,nou3)
            lb1='SVAVB2  '
            lb2='SVBOA2  '
            i1=1+nu2
            i2=i1+nu2
            itot=i2+nou
          call grdsym(nu,nu,nu,no,nu2,nou,v,ti,ti(i1),ti(i2),lb1,lb2,30)
            call tranmd(v,nu,nu,nu,no,23)
         else
            call getlst(v,1,nou,1,1,30)
            call tranmd(v,nu,nu,nu,no,23)
         endif
      endif
      if (io.eq.1)then
         if(nirrep.gt.1)then
            call zeroma(v,1,no3u)
            lb1='SOAOB2  '
            lb2='SOAVB2  '
            i1=1+no2
            i2=i1+no2
            itot=i2+nou
          call grdsym(nu,nu,nu,no,no2,nou,v,ti,ti(i1),ti(i2),lb1,lb2,10)
            call tranmd(v,nu,nu,nu,no,23)
         else
            call getlst(v,1,nou,1,1,10)
            call tranmd(v,nu,nu,nu,no,23)
         endif
      endif
      return
      end
      SUBROUTINE PCKINDsk(IVALUE,ISPIN,SYVEC1,SYVEC2,NSMSZ1,NSMSZ2,
     &                  I,J,A,B,list)
      IMPLICIT INTEGER (A-Z)
      DIMENSION SYVEC1(NSMSZ1),SYVEC2(NSMSZ2)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),N(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      ILRG(IX)=INT(0.5D0*(1.D0+DSQRT(8.d0*IX-7)))+1
      NNM1O2(IX)=(IX*(IX-1))/2
      IGETI(LSTELM,NUMA)=1+(LSTELM-1)/NUMA
      IGETA(LSTELM,NUMA)=LSTELM-(IGETI(LSTELM,NUMA)-1)*NUMA
      IF(IVALUE.Eq.0)THEN
       I=0
       J=0
       A=0
       B=0
       RETURN
      ENDIF
      NSIZ=0
      IBOTL=0
      IBOTR=0
      DO 10 IRREP=1,NIRREP
       NSYDIS=IRPDPD(IRREP,ISYTYP(2,LIST))
       NSYDSZ=IRPDPD(IRREP,ISYTYP(1,LIST))
       NSIZ=NSIZ+NSYDIS*NSYDSZ
       IF(IVALUE.LE.NSIZ)GOTO 20
       IBOTR=IBOTR+NSYDIS
       IBOTL=IBOTL+NSYDSZ
10    CONTINUE
20    NOFF=NSIZ-NSYDIS*NSYDSZ
      IOFF=IVALUE-NOFF
      NSYDIS=IRPDPD(IRREP,ISYTYP(2,LIST))
      NSYDSZ=IRPDPD(IRREP,ISYTYP(1,LIST))
      if(ioff.eq.0)then
        idisp=0
        ielemp=0
       else
      IDISP=IGETI(IOFF,NSYDSZ)
      IELEMP=IGETA(IOFF,NSYDSZ)
      endif
      IELEMU=SYVEC1(IBOTL+IELEMP)
      IDISU=SYVEC2(IBOTR+IDISP)
      if(list.eq.13)then
      J=1+(IDISU-1)/NOCCO(1)    !ijkL
      I=IDISU-(J-1)*NOCCO(1)    !ijKl
      B=1+(IELEMU-1)/NOCCO(1)   !iJkl
      A=IELEMU-(B-1)*NOCCO(1)   !Ijkl
      return
      endif
      if(list.eq.10)then
      J=1+(IDISU-1)/NOCCO(1)    !ijkA
      I=IDISU-(J-1)*NOCCO(1)    !ijKa
      B=1+(IELEMU-1)/NOCCO(1)   !iJka
      A=IELEMU-(B-1)*NOCCO(1)   !Ijka
      return
      endif
      if(list.eq.16)then
      J=1+(IDISU-1)/NOCCO(1)    !abiJ
      I=IDISU-(J-1)*NOCCO(1)    !abIj
      B=1+(IELEMU-1)/NVRTO(1)   !aBij
      A=IELEMU-(B-1)*NVRTO(1)   !Abij
      return
      endif
      if(list.eq.25)then
      J=1+(IDISU-1)/NVRTO(1)    !biaJ
      I=IDISU-(J-1)*NVRTO(1)    !biAj
      B=1+(IELEMU-1)/NVRTO(1)   !bIaj
      A=IELEMU-(B-1)*NVRTO(1)   !Biaj
      return
      endif
      if(list.eq.30)then
      J=1+(IDISU-1)/NVRTO(1)    !abcI
      I=IDISU-(J-1)*NVRTO(1)    !abCi
      B=1+(IELEMU-1)/NVRTO(1)   !aBci
      A=IELEMU-(B-1)*NVRTO(1)   !Abci
      endif
      if(list.eq.46)then        
      J=1+(IDISU-1)/NOCCO(1)    !abiJ
      I=IDISU-(J-1)*NOCCO(1)    !abIj
      B=1+(IELEMU-1)/NVRTO(1)   !aBij
      A=IELEMU-(B-1)*NVRTO(1)   !Abij
      return
      endif
      if(list.eq.233)then
      J=1+(IDISU-1)/NVRTO(1)    !abcD
      I=IDISU-(J-1)*NVRTO(1)    !abCd
      B=1+(IELEMU-1)/NVRTO(1)   !aBcd
      A=IELEMU-(B-1)*NVRTO(1)   !Abcd
      endif
      RETURN
      END
      subroutine izo(ia,n)
      implicit double precision(a-h,o-z)
      dimension ia(n)
      do 10 i=1,n
         ia(i)=0
 10   continue
      return
      end
      subroutine rifr(nfile,ifr,inr)
      implicit double precision (a-h,o-z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension ifr(1),inr(1)
      nrec=nu2/1024+1
      mm=mod(nu2,1024)
      iof=1
      do 20 i=1,nrec-1
         call  ird(nfile,i,1024,ifr(iof))
         iof=iof+1024
 20   continue
      call  ird(nfile,nrec,mm,ifr(iof))
      iof=1
      do 30 i=1,nrec-1
         call  ird(nfile,i+nrec,1024,inr(iof))
         iof=iof+1024
 30   continue
      call  ird(nfile,2*nrec,mm,inr(iof))
      return
      end
      subroutine wwp(nfile,n,irc,v)
      implicit double precision(a-h,o-z)
      dimension v(n)
      common/buforv/buf(341),ibuf(341),iiiii
      data tr/0.1d-13/
      in=0
      nr=0
      irec=irc
      do 10 i=1,n
         x=v(i)
         if (dabs(x).gt.tr)then
            in=in+1
            buf(in)=x
            ibuf(in)=i
         endif
         if (in.eq.341)then
            nr=nr+1
            write(nfile,rec=irec)buf,ibuf,in
            in=0
            irec=irec+1
         endif
 10   continue
      write(nfile,rec=irec)buf,ibuf,in
      irc=nr+1
      return
      end
      subroutine wifr(nfile,ifr,inr)
      implicit double precision (a-h,o-z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension ifr(1),inr(1)
      nrec=nu2/1024+1
      mm=mod(nu2,1024)
      iof=1
      do 20 i=1,nrec-1
         call  iwrt(nfile,i,1024,ifr(iof))
         iof=iof+1024
 20   continue
      call  iwrt(nfile,nrec,mm,ifr(iof))
      iof=1
      do 30 i=1,nrec-1
         call  iwrt(nfile,i+nrec,1024,inr(iof))
         iof=iof+1024
 30   continue
      call  iwrt(nfile,2*nrec,mm,inr(iof))
      return
      end
      subroutine ird(npak,i,n,ia)
      implicit double precision(a-h,o-z)
      common/buforv/ibuf(1024)
      dimension ia(n)
      read(npak,rec=i)ibuf
      call icopymm(ibuf,ia,n)
      return
      end
      subroutine iwrt(npak,i,n,ia)
      implicit double precision(a-h,o-z)
      common/buforv/ibuf(1024)
      dimension ia(n)
      call icopymm(ia,ibuf,n)
      write(npak,rec=i)ibuf
      return
      end
      SUBROUTINE rdsym(n1,n2,n3,n4,nsq1,nsq2,t2,t2s,
     *SYVEC1,SYVEC2,nlist)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SYVEC1(1),SYVEC2(1),t2(n1,n2,n3,n4),t2s(1)
      integer syvec1,syvec2,dirprd
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),N(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      no=nocco(1)
      nu=nvrto(1)
      kkk=0
      kij=0
      DO 10 IRREP=1,NIRREP
         NSYDIS=IRPDPD(IRREP,ISYTYP(2,nLIST))
         NSYDSZ=IRPDPD(IRREP,ISYTYP(1,nLIST))
         NSIZ=NSIZ+NSYDIS*NSYDSZ
         do 9 inum=1,nsydis
         kij=kij+1
         call getlst(t2s,inum,1,1,irrep,nlist)
         do 8 is=1,nsydsz
            kkk=kkk+1
            call pckindsk(kkk,3,syvec1,syvec2,nsq1,nsq2,i,j,ia,ib,nlist)
            t2(ia,ib,i,j)=t2s(is)
 8       continue
 9    CONTINUE
 10   CONTINUE
      RETURN
      END
      subroutine rdhhhhsm(no,ti,v)
      implicit double precision(a-h,o-z)
      character*8 lb1,lb2
      dimension ti(1),v(no,no,no,no)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      integer a,b,dirprd
         if(nirrep.gt.1)then 
            call zeroma(v,1,no4)
            lb1='SOAOB2  '
            i1=1+no2
            i2=i1+no2
            itot=i2+nou
          call grdsym(no,no,no,no,no2,no2,v,ti,ti(i1),ti(i2),lb1,lb1,13)
         else
            call getlst(v,1,no2,1,1,13)
         endif
      return
      end
      subroutine rdvemsm(io,no,nu,ti,v)
      implicit double precision(a-h,o-z)
      integer dirprd
      character*8 lb1,lb2
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      dimension ti(1),v(1)
      nou=no*nu
      if(io.eq.0)then
         if(nirrep.gt.1)then
            call zeroma(v,1,nou3)
            lb1='SVAVB2  '
            lb2='SVBOA2  '
            i1=1+nu2
            i2=i1+nu2
            itot=i2+nou
          call grdsym(nu,nu,nu,no,nu2,nou,v,ti,ti(i1),ti(i2),lb1,lb2,30)
            call tranmd(v,nu,nu,nu,no,23)
         else
            call getlst(v,1,nou,1,1,30)
            call tranmd(v,nu,nu,nu,no,23)
         endif
      endif
      if (io.eq.1)then
         if(nirrep.gt.1)then
            call zeroma(v,1,no3u)
            lb1='SOAOB2  '
            lb2='SOAVB2  '
            i1=1+no2
            i2=i1+no2
            itot=i2+nou
          call grdsym(nu,nu,nu,no,no2,nou,v,ti,ti(i1),ti(i2),lb1,lb2,10)
            call tranmd(v,nu,nu,nu,no,23)
         else
            call getlst(v,1,nou,1,1,10)
            call tranmd(v,nu,nu,nu,no,23)
         endif
      endif
      return
      end
      subroutine readnu2(i,j,nu,v,ifr,inr)
      implicit double precision(a-h,o-z)
      dimension v(1),ifr(1),inr(1)
      common/buforv/buf(341),ibuf(341),iiiii
      nu2=nu*nu
      call zeroma(v,1,nu2)
         ikol=(i-1)*nu+j
         ifst=ifr(ikol)
         nrec=inr(ikol)
         do 40 ire=1,nrec
            irec=ifst+ire-1
            read(34,rec=irec)buf,ibuf,in
            do 41 ii=1,in
               k=ibuf(ii)
               v(k)=buf(ii)
 41         continue
 40      continue
         call transq(v,nu)
      return
      end
      subroutine sym2v(nu,v)
      implicit double precision(a-h,o-z)
      integer a,b,e,f
      dimension v(1)
      data zero/0.0d0/,two/2.0d+0/,half/0.5d+0/      
      iabef=0
      nu2=nu*nu
      nu3=nu*nu2
      do 10 a=1,nu
         do 11 b=1,nu
            do 12 e=1,nu
               do 13 f=e,nu
                  iabef=iabef+1
                  iabefor=(a-1)*nu3+(b-1)*nu2+(e-1)*nu+f
                  v(iabef)=v(iabefor)
 13            continue
 12         continue
 11      continue
 10   continue
      return
      end
      subroutine sym2o2(no,nu,o2)
      implicit double precision(a-h,o-z)
      integer e,f
      dimension o2(1)
      data  half/0.5d+0/      
      ii=0
      no2=no*no
      do 12 e=1,nu
         ie=(e-1)*nu
         do 13 f=e,nu
            ii=ii+1
            ief=(ii-1)*no2+1
            iefor=(ie+f-1)*no2+1
            if(e.eq.f)call vecmul(o2(iefor),no2,half)
            call veccop(no2,o2(ief),o2(iefor))
 13      continue
 12   continue
      return
      end
      subroutine getv4(nu,lrec,ti,v,vec,ivec,ifr,inr)
      implicit double precision(a-h,o-z)
      integer a,b
      dimension ti(1),v(1),vec(1),ivec(1),ifr(1),inr(1)
      common/buforv/buf(341),ibuf(341),iiiii
      common/wpak/nfr(30),nsz(30)
      common/pak/nv4,intr
      common/lnt/lrecb,lrecw
      common/newio/nvv,nhh,npp,nvo,nt3,no1,nu1
      data tr/0.1d-12/
      irec=0
      nfr(1)=1
      nu2=nu*nu
      ijkl=0
      in=0
      intotal=0
      call rifr(34,ifr,inr)
 450  format(18i4)
      do 10 i=1,nu
      do 11 j=1,nu
         call readnu2(i,j,nu,v,ifr,inr)
         do 25 l=1,nu2
            ijkl=ijkl+1
            x=v(l)
            if(dabs(x).gt.tr)then
               intotal=intotal+1
               in=in+1
               vec(in)=x
               ivec(in)=ijkl
               if(in.eq.lrec)then
                  irec=irec+1
                  call ww(irec,lrec,vec,ivec)
                  in=0
               endif
            endif
 25      continue
 11      continue
 10   continue
      if(in.ne.0)then
         irec=irec+1
         call ww(irec,in,vec,ivec)
         endif
      nsz(1)=intotal
      nfr(2)=irec+1
      return
      end
      subroutine grdsym(n1,n2,n3,n4,nsq1,nsq2,t2,t2s,s1,s2,lb1,lb2,nlst)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*8 lb1,lb2
      integer s1,s2
      dimension t2(1),t2s(1),s1(nsq1),s2(nsq2)
      call getrec(20,'JOBARC',lb1,nsq1,s1)
      call getrec(20,'JOBARC',lb2,nsq2,s2)
      call rdsym(n1,n2,n3,n4,nsq1,nsq2,t2,t2s,s1,s2,nlst)
      return
      end
