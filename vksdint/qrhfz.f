      subroutine qrhfz(xova,x1ia,x1ib,
     &                 eval,icore,maxcor,iuhf)

c  this routine evaluates the following quantities for 
c  HF-KS     gradient calculations.

c   z(i,j) = x(i,j)/(eval(i,rhf)-eval(j,rhf)

c  in addition, the modified x intermediate (xtw) is formed

c  xtw(a,i) = x(a,i) - sum z(k,j) a(kj,ai) 
C                      k,J                   

      implicit double precision (a-h,o-z)
      integer dirprd,dissyw,disful,pop,vrt,poprhf,vrtrhf
      integer popdoc,vrtdoc,a,b,c,d
      dimension xova(1),xovb(1),eval(1),icore(maxcor)
      dimension x1ia(1),x1ib(1),xa2a(1),xa2b(1)

      common/machsp/iintln,ifltln,iintfp,ialone,ibitwd
      common/flags/iflags(100)
      common/sym/pop(8,2),vrt(8,2),nt(2),nf1(2),nf2(2) 
      common/syminf/nstart,nirrep,irrepa(255,2),dirprd(8,8)
      common/sympop/irpdpd(8,22),isytyp(2,500),ntot(18)
      common/info/nocco(2),nvrto(2)


      data zero,one,onem,two /0.d0,1.d0,-1.d0,2.d0/
      mxcor=maxcor
      norb=nocco(1)+nvrto(1) 
      write(*,*) 'norb=',norb
      call getrec(20,'JOBARC','SCFEVLA0',IINTFP*NORB,EVAL)

      ipos=1

      do 10 irrep=1,nirrep
       nocca  =pop(irrep,1)
       noccb  =pop(irrep,2)
       do 100 i=1,nocca
         eigi=eval(i)
         do 110 j=1,nocca
          eigj=eval(j)
          if( j .ne. i) then
          denom=one/(eigi-eigj)
          x1ia(ipos)=-x1ia(ipos)*denom
          else
             x1ia(ipos)=zero             
             x1ib(ipos)=zero               
          end if
          ipos=ipos+1
110      continue 
100     continue
c       call kkk(nocca*nocca,x1ia)
C  XTW(A,I)=X(A,I)
c          - sum Zkj(kM) [<kA||MI> + <kI||MA>] 
C            k,M
        write(*,*)'xov'
       call qrhfxtw1(xova,xovb,x1ia,
     &               irrep,icore,maxcor)
c       call kkk(9,xova)
10    continue
      return
      end
