      subroutine trnsfockrohf(cp,cp2,fvv,fvv2,fnodrop,scr,maxcor,ncp)
c     For ROHF references
c         so      so    no-d           no-d
c    no-d C'~  so F  so C'     =  no-d F'
c
c==============variable declarations and common blocks======================
      implicit none
C     Common blocks
      integer iintln,ifltln,iintfp
      common/machsp/iintln,ifltln,iintfp
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      integer pop(8,2),vrt(8,2),nt(2),nd1(2),nd2(2)
      common/sym/pop,vrt,nt,nd1,nd2
      integer naobasis,naobas(8)
      common/aoinfo/naobasis,naobas
C     Input Variables
      integer maxcor,fnodrop(8,2),ncp(2)
      double precision cp(*),cp2(*)
C     Output Variables
      double precision fvv(*),fvv2(*)
C     Pre-allocated Local Variables
      double precision scr(maxcor)
C     Local variables
      integer virt,drop,vrtdrop,spin,i000,i010,i020,i030,irrep,icp,
     &   icp2,ii,ifvv,ifvv2,ifock,aos
      double precision one,zilch
      character*8 cfock(2)
      data cfock /'FOCKA   ','FOCKB   '/
      parameter(one=1.d0)
      parameter(zilch=0.d0)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do 10 spin=1,2
         icp=1+ncp(1)*(spin-1)
         icp2=1+ncp(1)*(spin-1)
         ifvv=1+nd2(1)*(spin-1)
         ifvv2=1+nd2(1)*(spin-1)
         i000=1
         i010=i000+naobasis**2
         if(i010.gt.maxcor) call insmem('trnsfock',i030,maxcor)
C     Get AO Fock matrix
         call getrec(20,'JOBARC',cfock(spin),naobasis**2*iintfp,
     &      scr(i000))
         ifock=1
         do 20 irrep=1,nirrep
            virt=vrt(irrep,spin)
            drop=fnodrop(irrep,spin)
            vrtdrop=virt-drop
            aos=naobas(irrep)
            i020=i010+aos**2
            i030=i020+virt*aos
            if(i030.gt.maxcor)call insmem('trnsfock',i030,maxcor)            
            call blkcpy2(scr(i000),naobasis,naobasis,scr(i010),aos,aos,
     &         ifock,ifock)
            if (vrtdrop.gt.0) then
               call xgemm('t','n',vrtdrop,aos,aos,one,cp(icp),aos,
     &            scr(i010),aos,zilch,scr(i020),vrtdrop)               
               call xgemm('n','n',vrtdrop,vrtdrop,aos,one,scr(i020),
     &            vrtdrop,cp(icp),aos,zilch,fvv(ifvv),vrtdrop)
            endif
            if (drop.gt.0) then
               call xgemm('t','n',drop,aos,aos,one,cp2(icp2),aos,
     &            scr(i010),aos,zilch,scr(i020),drop)
               call xgemm('n','n',drop,drop,aos,one,scr(i020),drop,
     &            cp2(icp2),aos,zilch,fvv2(ifvv2),drop)
            endif
            ifock=ifock+aos
            icp=icp+aos*vrtdrop
            icp2=icp2+aos*drop
            ifvv=ifvv+vrtdrop**2
            ifvv2=ifvv2+drop**2
 20      continue
 10   continue
      return
1000  format((8(2x,I9)))
2000  format(t3,70('-'))
      end
