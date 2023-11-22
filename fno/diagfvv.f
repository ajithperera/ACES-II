      subroutine diagfvv(fvv,fvv2,fvvevl,fvvevl2,fnodrop,cp,cp2,scr,
     &   maxcor,uhf,ncp)
c
c The reduced virtual virtual fock matrix is diagonalized
c
c       no-d        no-d        no-d          no-d
c  no-d Z~     no-d F'     no-d Z     =  no-d e'
c
c Z also contributes to C'
c matrix multiplication is done on the fly
c
c     no-d         no-d         no-d
c  so C'      no-d Z      =  so C'
c
c=================variable declarations and common blocks======================
c     implicit double precision (a-h,o-z)
      implicit none
C     Common Blocks
      integer nocco(2),nvrto(2)
      common/info/nocco,nvrto
      integer iintln,ifltln,iintfp
      common/machsp/iintln,ifltln,iintfp
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      integer pop(8,2),vrt(8,2),nt(2),nd1(2),nd2(2)
      common/sym/pop,vrt,nt,nd1,nd2
      integer naobasis,naobas(8)
      common/aoinfo/naobasis,naobas
      logical rohf,qrhf,semi
      common/fnoparam/rohf,qrhf,semi
C     Input Variables
      integer fnodrop(8,2),uhf,maxcor,ncp(2)
C     Output Variables
      double precision fvvevl(*),fvvevl2(*)
C     Input/Output Variables
      double precision fvv(*),fvv2(*),cp(*),cp2(*)
C     Pre-allocated Local variables
      double precision scr(maxcor)
C     Local variables
      integer vrtdrop,aos,drop,virt,irrep,i000,i010,ii,icp,icp2,spin,
     &   ifvvevl,ifvvevl2,ifvv,ifvv2,numdrop
      logical do_eig
      double precision one,zilch
      character*5 csptype(2)
      character*8 wlabel(2)
      data csptype /'alpha','beta '/
      data wlabel  /'WMATRXVA','WMATRXVB'/
      parameter(one=1.d0)
      parameter(zilch=0.d0)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      i000=1
      do_eig=(.not.rohf).or.semi
      do 10 spin=1,uhf+1
         icp=1+ncp(1)*(spin-1)
         icp2=1+ncp(1)*(spin-1)
         ifvv=1+nd2(1)*(spin-1)
         ifvv2=1+nd2(1)*(spin-1)
         ifvvevl=1+nvrto(1)*(spin-1)
         ifvvevl2=1+nvrto(1)*(spin-1)
         write(6,*)' ',csptype(spin),' spin'
         numdrop=0
         do 20 irrep=1,nirrep
            drop=fnodrop(irrep,spin)
            virt=vrt(irrep,spin)
            vrtdrop=virt-drop
            aos=naobas(irrep)
            numdrop=numdrop+vrtdrop**2
            if(vrtdrop.gt.0) then
               i010=i000+aos*vrtdrop
               if(i010.gt.maxcor) call insmem('diagfvv',i010,maxcor)
               if (do_eig) call eig(fvv(ifvv),scr(i000),1,vrtdrop,0)
               call dcopy(vrtdrop,fvv(ifvv),vrtdrop+1,fvvevl(ifvvevl),1)
               
               if (do_eig) then
                  call dcopy(vrtdrop**2,scr(i000),1,fvv(ifvv),1)
                  call xgemm('n','n',aos,vrtdrop,vrtdrop,one,cp(icp),aos
     &               ,fvv(ifvv),vrtdrop,zilch,scr(i000),aos)               
                  call dcopy(aos*vrtdrop,scr(i000),1,cp(icp),1)
               endif
            endif
            if(drop.gt.0)then
               i010=i000+aos*drop
               if(i010.gt.maxcor) call insmem('diagfvv',i010,maxcor)
               if (do_eig) call eig(fvv2(ifvv2),scr(i000),1,drop,0)
               call dcopy(drop,fvv2(ifvv2),drop+1,fvvevl2(ifvvevl2),1)
               
               if (do_eig) then
                  call dcopy(drop**2,scr(i000),1,fvv2(ifvv2),1)
                  call xgemm('n','n',aos,drop,drop,one,cp2(icp2),aos,
     &               fvv2(ifvv2),drop,zilch,scr(i000),aos)
                  call dcopy(aos*drop,scr(i000),1,cp2(icp2),1)
               endif
            endif
c//////////////////////////////////////////////////////////////////////////////
            if(virt.gt.0)then
               write(6,*)' new virtual orbital energies for irrep',irrep
               write(6,1000)
               write(6,2000)(fvvevl(ii),ii=ifvvevl,ifvvevl+vrtdrop-1)
               write(6,1000)
            endif
c////////////////////////////////////////////////////////////////////////////// 
            ifvv=ifvv+vrtdrop**2
            ifvv2=ifvv2+drop**2
            ifvvevl=ifvvevl+vrtdrop
            ifvvevl2=ifvvevl2+drop
            icp=icp+aos*vrtdrop
            icp2=icp2+aos*drop
 20      continue
         if (semi) call putrec(20,'JOBARC',wlabel(spin),numdrop*iintfp,
     &      fvv(1+nd2(1)*(spin-1)))
 10   continue
      return
 1000 format(t3,70('-'))
 2000 format((4(2x,f13.10)))
      end
