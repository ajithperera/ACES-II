      subroutine qrhfxtw1(xova,xovb,zija,
     &                    irrepz,icore,maxcor)
c  this routine forms contributions to the modified X
c  intermediate required for hf-ks   gradients. 
C
c   XTW(A,I) = X(A,I) - SUM Z(k,M) [<kM||AI> + <AM||kI>]
C                       k,M
C
c  here a,i are the virtual and occupied orbitals in the rhF
c       reference,
c       K,M ARE THE  doubly occupied orbitals
c       in the rhf reference
      
      implicit double precision (a-h,o-z)
      integer dirprd,dissyw,disful,pop,vrt,poprhf,vrtrhf
      integer popdoc,vrtdoc
      dimension xova(1),xovb(1),icore(maxcor)
      dimension ZijA(1),ZijB(1)

      common/machsp/iintln,ifltln,iintfp,ialone,ibitwd
      common/sym/pop(8,2),vrt(8,2),nt(2),nf1(2),nf2(2) 
      common/syminf/nstart,nirrep,irrepa(255,2),dirprd(8,8)
      common/sympop/irpdpd(8,22),isytyp(2,500),ntot(18)
      common/info/nocco(2),nvrto(2)

      INDXF(I,J,N)=I+(J-1)*N

      data zero,one,onem,two /0.d0,1.d0,-1.d0,2.d0/

      mxcor=maxcor

C TERM I:  SUM Z(1,M) * [<MI||1A> + <1I||MA>]  Z - AA   W - AAAA
C          1,M

C      if(ispinp.eq.2) then
C
C THE OPEN SHELL ORBITAL IS OCCUPIED IN THE ALPHA SPIN CASE,
C HENCE THE <1A||MI> INTEGRALS ARE OF TYPE <IJ||KA>
C
        LISTW=7
c      else
C
C THE OPEN SHELL ORBITAL IS OCCUPIED IN THE BETA SPIN CASE,
C HENCE THE <1A||MI> INTEGRALS ARE OF TYPE <ij||ka>
C
c         LISTW=8
c      ENDIF
C
C LOOP OVER IRREPS OF MI AND 1A

       numsyw=irpdpd(1,isytyp(2,listw))
       dissyw=irpdpd(1,isytyp(1,listw))

       n=irpdpd(1,18)
       d=irpdpd(1,14)

       disful=irpdpd(1,21)
       
       i000=1
       i010=i000+iintfp*disful*numsyw
       i020=i010+iintfp*n*D  
       if(i020.ge.maxcor)call insmem('qrhfxtw1',i020,maxcor)

       call getlst(icore(i000),1,numsyw,2,1,listw)

       call symexp2(1,pop(1,1),disful,dissyw,numsyw,
     &              icore(i000),icore(i000))
c       call kkk(disful*9,icore(i000))
       call getlst(icore(i010),1,n,1,1,10)

       call zjkApxia(icore(i000),zija,xova,pop(1,1),
     &               vrt(1,1),icore(i010))
      return
      end
