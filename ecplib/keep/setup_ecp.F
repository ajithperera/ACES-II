      subroutine set(namstr,iqmstr,jcostr,nucstr,nrcstr,NHARM)
c-----------------------------------------------------------------------      
c     subroutine doit fills the common-blocks for subroutine paloop
c     version for general contracted basissets
c     input:
c     namstr: names of the atoms (strange storage as integer!)
c     iqmstr: max. l-quantum-number+1 of a shell
c     jcostr: number of shells
c     nucstr: number of the uncontracted functions of the shells > NUCO
c             (exponents)
c     nrcstr: number of blocks(general contracted) > NRCO
c             (contractions)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
cch   symmetry values
      integer eor,and,or,dstrt
C-----------------------------------------------------------------------
C     parameter from param.h and lquant.h (TURBOMOLE)
C-----------------------------------------------------------------------
#include "baslims.par"
#include "ecplims.par"
C
C-----------------------------------------------------------------------
      parameter (nftmax = (ndilmx+1)*(ndilmx+2)*(ndilmx*3)/6)
C-----------------------------------------------------------------------
C     nftmax  total number of reducible gaussians up to {ndilmx+1}
C     SPECIAL parameter needed in gradient calculations
C-----------------------------------------------------------------------
      dimension iqmstr(ndi10),jcostr(ndi10,ndilmx)
C SG 8/25/98 namstr changed to character*4
      character*4 namstr(ndi10),namsml(ndi10)
      dimension nucstr(ndi10,ndilmx,maxjco),nrcstr(ndi10,ndilmx,maxjco) 
cch     dimension gout(20000,1)
cch      dimension ghelp(ndgout,1)
      dimension ipeoff(ndi10,ndilmx,maxjco),iccoff(ndi10,ndilmx,maxjco)
      dimension newoff(ndi10,ndilmx,maxjco),numpri(ndi10,ndilmx,maxjco)
      dimension imin(ndilmx),imax(ndilmx),nfcts(ndilmx)
      dimension ifeld(ndi10)
      dimension iadr(0:1,0:1,0:1)
cch   iadr is never used in this program cause mkadr=0!
      logical ja
      common /palcali/ nshdim
      common /palcal/ toll
c-----------------------------------------------------------------------
c     namsml              : sml = small
c     gout                : intbatch
c     ipeoff              : offset for primitive exponents
c     iccoff              : offset for contraction coefficients
c     newoff              : offsets for primitives
c     numpri              : number of primitive functions of a shell
c     imin,imax           : lower(upper) offset of a cartesian component
c     nfcts               : number of types of a cartesian component 
c     ifeld               : checkarray 
c-----------------------------------------------------------------------
      common /infoa / xyz(3,ndi10),charg(ndi10),wmass(ndi10),natoms
      common /symshe/ mulsh(ndi9)
      common /modez / zetm(ndi9)
      common /pairij/ ipq(ndi4+1)
      common /typinf/ ltmax,ltdrv,nftnft,nftdrv
      common /powers/ jx(nftmax),jy(nftmax),jz(nftmax),
     1                ix(nftmax),iy(nftmax),iz(nftmax)
c-----------------------------------------------------------------------
c     common from AcesII
c-----------------------------------------------------------------------
      COMMON /DAT/ ALPHA(MXTNPR),CONT(MXTNCC),CENT(3,MXTNSH),
     &             CORD(100,3),CHARGE(100),FMULT(8),TLA,TLC 
      common/indx/pc(512),dstrt(8,MXCBF),ntap,lu2,nrss,nucz,itag,
     &            maxlop,maxlot,kmax,nmax,khkt(7),mult(8),isytyp(3),
     &            itype(7,28),and(8,8),or(8,8),eor(8,8),nparsu(8),
     &            npar(8),mulnuc(100),nhkt(MXTNSH),mul(MXTNSH),
     &            nuco(MXTNSH),nrco(MXTNSH),jstrt(MXTNSH),
     &            nstrt(MXTNSH),mst(MXTNSH),jrs(MXTNSH)
      common/xa/xand(3,8,nh4)
CDEBUG
c-----------------------------------------------------------------------
c     write(*,*) 'Exponenten und Kontr.koeff.'
c     write(*,*)
c     do 2345 i=1,100
c       write(*,*) i,alpha(i)
c2345 continue
c     write(*,*) 
c     do 2346 i=1,100
c       write(*,*) i,cont(i)
c2346 continue
c     write(*,*)
c-----------------------------------------------------------------------
CDEBUG
c-----------------------------------------------------------------------
c fill xand in /xa/
c-----------------------------------------------------------------------
      do 2000 j=1,3
        do 2100 i=1,maxlop
          xand(j,i,1)=1.
          do 2200 k=2,nh4-1,2
            xand(j,i,k+1)=1.
            xand(j,i,k)=pc(and(isytyp(j),i))
 2200     continue
 2100   continue
 2000 continue
c-----------------------------------------------------------------------
c     first create some arrays for help
c-----------------------------------------------------------------------

c     detect the different types of atoms

      do i = 1, ndi10
         namsml(i) = '    '
      end do

c     find offsets of primitive exponents and contraction coefficients
c     --> ipeoff, iccoff

CDEBUG
c     write(*,*)
c     write(*,*) 'jrs'
c     write(*,*) jrs
c     write(*,*) 'jstrt'
c     write(*,*) jstrt
c     write(*,*)
CDEBUG
      icount=0
      call izero(ifeld,ndi10)
      call izero(ipeoff,ndi10*ndilmx*maxjco)  
      call izero(iccoff,ndi10*ndilmx*maxjco)  
      do 1200 iat=1,NMAX
        do 1250 ilq=1,iqmstr(iat)
          do 1300 ijco=1,jcostr(iat,ilq)
            icount=icount+1
            if (ifeld(iat).ne.1) then
              ipeoff(iat,ilq,ijco)=jstrt(icount)
              iccoff(iat,ilq,ijco)=jrs(icount)
CDEBUG
c             write(*,1) iat,ilq,ijco,ipeoff(iat,ilq,ijco)
c   1         format('iat: ',i3,'ilq: ',i3,' ijco: ',i3,' ipeoff: ',i3)
c             write(*,2) iat,ilq,ijco,iccoff(iat,ilq,ijco)
c   2         format('iat: ',i3,'ilq: ',i3,' ijco: ',i3,' iccoff: ',i3)
CDEBUG
            endif
 1300     continue
 1250   continue
        ifeld(iat)=1 
 1200 continue 

c     create shell-structure
c     and arrays expnt and contr

      call zero(expnt,ndi13)
c      call zero(contr,ndi13)
      call zero(contr,ndi13*ndico)
      call izero(numpri,ndi10*ndilmx*maxjco)
      icount=0
      kcount=0
c     lcount=1
      ncount=0
      do 1500 iat=1,NMAX
        do 1550 ilq=1,iqmstr(iat)
          do 1575 ijco=1,jcostr(iat,ilq)
            ncount=ncount+1
            do 1600 inuco=1,nucstr(iat,ilq,ijco)
              icount=icount+1
              do 1650 inrco=1,nrcstr(iat,ilq,ijco)
                contr(icount,inrco)=CONT(iccoff(iat,ilq,ijco)+inuco+   
     &          nucstr(iat,ilq,ijco)*(inrco-1))
CDEBUG
c                write(*,3) icount,inrco,contr(icount,inrco)
CDEBUG
   3            format('icount: ',i3,' inrco: ',i3,' contr: ',f12.8)
 1650         continue
              expnt(icount)=ALPHA(ipeoff(iat,ilq,ijco)+inuco)
CDEBUG
c             write(*,4) icount,expnt(icount)
CDEBUG
   4          format('icount: ',i3,' expnt: ',f12.8)
              kcount=kcount+1
 1600       continue
c           numpri(iat,ilq,ijco)=numpri(iat,ilq,ijco)+kcount
            numpri(iat,ilq,ijco)=kcount
CDEBUG
c           write(*,5) iat,ilq,ijco,numpri(iat,ilq,ijco)
CDEBUG
   5        format('iat: ',i3,' ilq: ',i3,' ijco: ',i3,' numpri: ',i3)
            kcount=0
            numcon(ncount)=nrcstr(iat,ilq,ijco)
CDEBUG
c           write(*,6) ncount,numcon(ncount)
CDEBUG
   6        format('ncount: ',i3,' numcon: ',i3)
 1575     continue
c         numpri(iat,ilq)=numpri(iat,ilq)+kcount
c         write(*,5) iat,ilq,numpri(iat,ilq)
c5        format('iat: ',i3,' ilq: ',i3,' numpri: ',i3)
c         lcount=lcount+1
c         kcount=0
 1550   continue
c       lcount=0
 1500 continue

CDEBUG
c     write(*,*)
c     stop 'in doit!'
CDEBUG
c     create new offsets of primitive exponents

      call izero(newoff,ndi10*ndilmx*maxjco)
      khelp=1
c     kcount=0
      do 3000 iat=1,NMAX
c       kcount=1
        do 3100 ilq=1,iqmstr(iat)
          do 3200 ijco=1,jcostr(iat,ilq)
c           newoff(iat,kcount)=khelp
            newoff(iat,ilq,ijco)=khelp
CDEBUG
c           write(*,99) iat,ilq,ijco,newoff(iat,ilq,ijco)
  99        format('iat: ',i3,' ilq: ',i3,' ijco: ',i3,' newoff: ',i3)
CDEBUG
c           khelp=khelp+numpri(iat,kcount)
            khelp=khelp+numpri(iat,ilq,ijco)
c           kcount=kcount+1
 3200     continue
 3100   continue
 3000 continue

c-----------------------------------------------------------------------
c     fill in common infoa
c-----------------------------------------------------------------------

      call zero(xyz,3*ndi10)
      natoms=NMAX
      do 7000 i=1,natoms
        do 7100 j=1,3
          xyz(j,i)=CORD(i,j)
 7100   continue
 7000 continue
      do 8000 i=1,natoms
        charg(i)=CHARGE(i)
 8000 continue

c-----------------------------------------------------------------------
c     fill in common /typinf/
c-----------------------------------------------------------------------

      ltmax=1
      do 10 iat=1,NMAX
        if (ltmax.le.iqmstr(iat)) then 
          ltmax=iqmstr(iat)
        endif  
   10 continue      
      nftdrv=ltmax*(ltmax+1)*(ltmax+2)/6

c-----------------------------------------------------------------------
c     fill in common /nshell/
c-----------------------------------------------------------------------

c Doing the k's from common nshel
      
      if (NHARM.eq.0) then

c spherical harmonics are used!

        do l=1,ndilmx
           imin(l)=(l-1)*l*(l+1)/6+1
           imax(l)=l*(l+1)*(l+2)/6 
           nfcts(l)=2*l-1
        enddo

      else 

        do l=1,ndilmx
           imin(l)=(l-1)*l*(l+1)/6+1
           imax(l)=l*(l+1)*(l+2)/6 
           nfcts(l)=l*(l+1)/2
        enddo
      endif

      nbf=0
      nshll=0
      indsh=0
C
c loop over atoms

      do iat=1,NMAX

c loop over l-quantum-numbers

        do ilq=1,iqmstr(iat)

c loop over contracted basisfunctions

          do ijco=1,jcostr(iat,ilq)

             indsh=indsh+1
             mulilq=mul(indsh)
             nshll=nshll+1
             katom(nshll)=iat
             ktype(nshll)=ilq
             kmini(nshll)=imin(ilq)
             kmaxi(nshll)=imax(ilq)
             ndegen(nshll)=nfcts(ilq)             
             kbfn(nshll)=nbf+1
             mulsh(nshll)=mulilq
             nprims(nshll)=numpri(iat,ilq,ijco)
             kprim(nshll)=newoff(iat,ilq,ijco)
             nbf=nbf+ndegen(nshll)*fmult(mulilq)

          Enddo
        Enddo
      Enddo

      nshell=nshll

CDEBUG
c      write(*,*)
c      do 3912 i=1,nshell
c       write(*,*) 'nshell=',i
c       write(*,3913) i,katom(i),i,ktype(i)
c3913   format('katom<',i2,'>  ',i2,' ktype<',i2,'>  ',i2)
c       write(*,3914) i,kprim(i),i,kbfn(i)
c3914   format('kprim<',i2,'>  ',i2,' kbfn<',i2,'>   ',i2)
c       write(*,3915) i,nprims(i),i,ndegen(i)
c3915   format('nprims<',i2,'> ',i2,' ndegen<',i2,'> ',i2)
c       write(*,3916) i,kmini(i),i,kmaxi(i)
c3916   format('kmini<',i2,'>  ',i2,' kmaxi<',i2,'>  ',i2)
c       write(*,*)
c3912  continue
c      write(*,*) 'nbf=',nbf
c      write(*,*)
cCDEBUG

c-----------------------------------------------------------------------
c     fill in common modez
c-----------------------------------------------------------------------
      do 4000 ii=1,nshell
        zetm(ii)=expnt(kprim(ii))
        do 4100 ipr=kprim(ii)+1,kprim(ii)+nprims(ii)-1
          zetm(ii)=min(expnt(ipr),zetm(ii))
 4100   continue
 4000 continue

CDEBUG
c     write(*,*)
c     do 7398 i=1,nshell
c       write(*,7399) i,zetm(i)
c7399   format('zetm<',i2,'> ',f13.8)
c7398 continue
c     write(*,*)
CDEBUG

c-----------------------------------------------------------------------
c     fill in common pairij
c-----------------------------------------------------------------------

      call setipq(ipq,ndi4)

c-----------------------------------------------------------------------
c     fill in common pairij
c-----------------------------------------------------------------------

      call monom(ltmax,nft,jx,jy,jz,iadr,nftmax,1,0)

c-----------------------------------------------------------------------
c     fill common /palcal/
c-----------------------------------------------------------------------
      nshdim=nshell
c      TLA=1.0d-12
c      toll=-log10(TLA)+8.0
       toll=25.0    
c-----------------------------------------------------------------------
c     Ende 
c-----------------------------------------------------------------------
      return
      end 
