      subroutine gett3s(no,nu,o3,w1,w2,t3,isym,iadw2,iadblk,len)
      implicit double precision (a-h,o-z)
      integer a,b,c,pop,vrt,dirprd
      logical ijeql,noneql
      common/newio/ni1,ni2,ni3,ni4,ntt3,no1,nu1
      dimension o3(nu,nu,nu),t3(nu,nu,nu),w1(1),w2(1),iadw2(1),
     *iadblk(1),isym(50,2),len(8,8)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
      COMMON /T3IOOF/ IJKPOS(8,8,8,2),IJKLEN(36,8,4),IJKOFF(36,8,4),
     1                NCOMB(4)
      index(i)=i*(i-1)/2
      if(no.gt.50)then 
      write(6,*)'you have number of correlated electrons
     *           greater than 100'
      stop
      endif
      call zeroma(o3,1,nu3)
      nu3=nu*nu*nu
      do 8 irpijk=1,nirrep
         do 7 irpc=1,nirrep
            irpab=dirprd(irpijk,irpc)
            len(irpc,irpijk)=irpdpd(irpab,1)*vrt(irpc,1)
 7       continue
 8    continue
      kk=0
      do 5 i=1,nirrep
         do 6 j=1,pop(i,1)
            kk=kk+1
            isym(kk,1)=i
            isym(kk,2)=j
 6       continue
 5    continue
      do 100 inew=1,no
         do 99 jnew=1,inew
            do 98 knew=1,jnew
               if(inew.eq.knew)goto 98
               ijeql=.false.
               noneql=.false.
               irpi=isym(inew,1)
               irpj=isym(jnew,1)
               irpk=isym(knew,1)
               irpij =dirprd(irpi,irpj)
               irpijk=dirprd(irpij,irpk)
               do 2 irrep=1,nirrep
                  if (irrep.eq.1)then
                     iadblk(irrep)=1
                  else
                     iadblk(irrep)=iadblk(irrep-1)+len(irrep-1,irpijk)
                  endif
 2             continue
               irpk=isym(jnew,1)
               k   =isym(jnew,2)
               irpj=isym(inew,1)
               j   =isym(inew,2)
               irpi=isym(knew,1)
               i   =isym(knew,2)
               irpjk =dirprd(irpj,irpk)
               irpijk=dirprd(irpi,irpjk)
               ijk=ijkpos(irpi,irpj,irpk,2)
               ioff=ijkoff(ijk,irpijk,2)
               if(irpi.eq.irpj)ijeql=.true.
               noneql=.not.ijeql
               if(ijeql)nij=(pop(irpi,1)*(pop(irpj,1)-1))/2
               if(noneql)nij=pop(irpi,1)*pop(irpj,1)
               if(ijeql) ijkval=ioff+(k-1)*nij+index(j-1)+i
               if(noneql)ijkval=ioff+(k-1)*nij+(j-1)*pop(irpi,1)+i
               call zeroma(o3,1,nu3)
               call getlist(o3,ijkval,1,1,irpijk,2)
c               call getlist(o3,ijkval,1,1,irpijk,6)
c               write(6,*)'getlist,ijkval,irpijk',ijkval,irpijk
c               call chksumsk('getls 1 ',o3,nu3)
               call symtrw2m(o3,w1,w2,iadblk,iadw2,irpijk,1,2)
               call trant3(o3,nu,3)
               call zeroma(t3,1,nu3)
               if (inew.eq.jnew)then
                  call adaiij(nu,o3,t3)
                  goto 999
               endif
               if(jnew.eq.knew)then
                  call adaijj(nu,o3,t3)
                  goto 999
               endif
ckij
               irpk=isym(knew,1)
               k   =isym(knew,2)
               irpi=isym(jnew,1)
               i   =isym(jnew,2)
               irpjk =dirprd(irpj,irpk)
               irpijk=dirprd(irpi,irpjk)
               ijk=ijkpos(irpi,irpj,irpk,2)
               ioff=ijkoff(ijk,irpijk,2)
               if(irpi.eq.irpj)ijeql=.true.
               noneql=.not.ijeql
               if(ijeql)nij=(pop(irpi,1)*(pop(irpj,1)-1))/2
               if(noneql)nij=pop(irpi,1)*pop(irpj,1)
               if(ijeql) ijkval=ioff+(k-1)*nij+index(j-1)+i
               if(noneql)ijkval=ioff+(k-1)*nij+(j-1)*pop(irpi,1)+i
               call zeroma(o3,1,nu3)
               call getlist(o3,ijkval,1,1,irpijk,2)
               call symtrw2m(o3,w1,w2,iadblk,iadw2,irpijk,1,2)
               call trant3(o3,nu,3)
               call zeroma(t3,1,nu3)
               call adakij(nu,o3,t3)
c
               irpk=isym(jnew,1)
               k   =isym(jnew,2)
               irpj=isym(inew,1)
               j   =isym(inew,2)
               irpi=isym(knew,1)
               i   =isym(knew,2)
               irpjk =dirprd(irpj,irpk)
               irpijk=dirprd(irpi,irpjk)
               ijk=ijkpos(irpi,irpj,irpk,2)
               ioff=ijkoff(ijk,irpijk,2)
               ijeql=.false.
               noneql=.false.
               if(irpi.eq.irpj)ijeql=.true.
               noneql=.not.ijeql
               if(ijeql)nij=(pop(irpi,1)*(pop(irpj,1)-1))/2
               if(noneql)nij=pop(irpi,1)*pop(irpj,1)
               if(ijeql) ijkval=ioff+(k-1)*nij+index(j-1)+i
               if(noneql)ijkval=ioff+(k-1)*nij+(j-1)*pop(irpi,1)+i
               call zeroma(o3,1,nu3)
               call getlist(o3,ijkval,1,1,irpijk,2)
               call symtrw2m(o3,w1,w2,iadblk,iadw2,irpijk,1,2)
               call trant3(o3,nu,3)
               call adajik(nu,o3,t3)
               ijeql=.false.
               noneql=.false.
               irpk=isym(inew,1)
               k   =isym(inew,2)
               irpj=isym(jnew,1)
               j   =isym(jnew,2)
               irpi=isym(knew,1)
               i   =isym(knew,2)
               irpjk =dirprd(irpj,irpk)
               irpijk=dirprd(irpi,irpjk)
               ijk=ijkpos(irpi,irpj,irpk,2)
               ioff=ijkoff(ijk,irpijk,2)
               if(irpi.eq.irpj)ijeql=.true.
               noneql=.not.ijeql
               if(ijeql)nij=(pop(irpi,1)*(pop(irpj,1)-1))/2
               if(noneql)nij=pop(irpi,1)*pop(irpj,1)
               if(ijeql) ijkval=ioff+(k-1)*nij+index(j-1)+i
               if(noneql)ijkval=ioff+(k-1)*nij+(j-1)*pop(irpi,1)+i
               call zeroma(o3,1,nu3)
               call getlist(o3,ijkval,1,1,irpijk,2)
               call symtrw2m(o3,w1,w2,iadblk,iadw2,irpijk,1,2)
               call trant3(o3,nu,3)
               call adaijk(nu,o3,t3)
 999           continue
               ias=it3(inew,jnew,knew)
               call vpakt3new(ias,i,j,k,nu3,t3)
 98         continue
 99      continue
 100  continue
      return
      end
      subroutine adaiij(nu,o3,t3)
      implicit double precision (a-h,o-z)
      integer a,b,c
      dimension o3(nu,nu,nu),t3(nu,nu,nu)
      nu3=nu*nu*nu
      call zeroma(t3,1,nu3)
      do 100 a=1,nu
         do 99 b=1,a
            b1=b
            if(a.eq.b)b1=b-1
            do 98 c=1,b1
               t3(a,c,b)=-o3(a,b,c)
               if(b.ne.c)t3(c,a,b)=t3(c,a,b)+o3(c,a,b)
               t3(c,a,b)= o3(c,a,b)-o3(b,a,c)
               if(b.ne.c)t3(c,b,a)=-o3(b,a,c)
               if(a.ne.b)t3(b,c,a)=t3(b,c,a)-o3(b,a,c)
 98         continue
 99      continue
 100  continue
      return
      end
      subroutine adaijj(nu,o3,t3)
      implicit double precision (a-h,o-z)
      integer a,b,c
      dimension o3(nu,nu,nu),t3(nu,nu,nu)
      nu3=nu*nu*nu
      call zeroma(t3,1,nu3)
      do 100 a=1,nu
         do 99 b=1,a
            b1=b
            if(a.eq.b)b1=b-1
            do 98 c=1,b1
               if(a.ne.b)t3(c,b,a)=-o3(b,a,c)
               if(b.ne.c)then
               t3(b,a,c)=-o3(c,a,b)
               t3(b,c,a)=-o3(c,a,b)
               endif
               t3(c,a,b)=-o3(b,a,c)
 98         continue
 99      continue
 100  continue
      return
      end
      subroutine adaijk(nu,o3,t3)
      implicit double precision (a-h,o-z)
      integer a,b,c
      dimension o3(nu,nu,nu),t3(nu,nu,nu)
      nu3=nu*nu*nu
      do 100 a=1,nu
         do 99 b=1,a
            b1=b
            if(a.eq.b)b1=b-1
            do 98 c=1,b1
               t3(a,c,b)=-o3(a,b,c)
               if(b.ne.c)then
                  t3(c,a,b)=t3(c,a,b)+o3(c,a,b)
                  if(a.ne.b)t3(b,c,a)=t3(b,c,a)-o3(b,a,c)
               endif
 98         continue
 99      continue
 100  continue
      return
      end
c jik
      subroutine adajik(nu,o3,t3)
      implicit double precision (a-h,o-z)
      integer a,b,c
      dimension o3(nu,nu,nu),t3(nu,nu,nu)
      nu3=nu*nu*nu
      do 100 a=1,nu
         do 99 b=1,a
            b1=b
            if(a.eq.b)b1=b-1
            do 98 c=1,b1
               if(b.ne.c)t3(c,a,b)=t3(c,a,b)-o3(b,a,c)
               if(a.ne.b)t3(c,b,a)=-o3(b,a,c)
 98         continue
 99      continue
 100  continue
      return
      end
      subroutine adakij(nu,o3,t3)
      implicit double precision (a-h,o-z)
      integer a,b,c
      dimension o3(nu,nu,nu),t3(nu,nu,nu)
      nu3=nu*nu*nu
      do 100 a=1,nu
         do 99 b=1,a
            b1=b
            if(a.eq.b)b1=b-1
            do 98 c=1,b1
               t3(b,a,c)=-o3(c,a,b)
               t3(b,c,a)=t3(b,c,a)-o3(c,a,b)
 98         continue
 99      continue
 100  continue
      return
      end
      subroutine gettest(o3,nu)
      implicit double precision (a-h,o-z)
      dimension o3(1)
      nu3=nu*nu*nu
      call zeroma(o3,1,nu3)
      call getlist(o3,1,1,1,1,2)
      return
      end
