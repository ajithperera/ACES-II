      subroutine trnsfock(u,u2,fvv,fvv2,fnodrop,scr,maxcor,uhf)
c
c         mo      mo    no-d           no-d
c    no-d U'~  mo F  mo U'     =  no-d F'
c
c==============variable declarations and common blocks======================
      implicit none
C     Common blocks
      integer nocco(2),nvrto(2)
      common/info/nocco,nvrto
      integer iintln,ifltln,iintfp
      common/machsp/iintln,ifltln,iintfp
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      integer pop(8,2),vrt(8,2),nt(2),nd1(2),nd2(2)
      common/sym/pop,vrt,nt,nd1,nd2
C     Input Variables
      integer maxcor,uhf,fnodrop(8,2)
      double precision u(*),u2(*)
C     Output Variables
      double precision fvv(*),fvv2(*)
C     Pre-allocated Local Variables
      double precision scr(maxcor)
C     Local variables
      integer virt,drop,vrtdrop,spin,i000,i010,i020,i030,irrep,nmo,iu,
     &   iu2,ii,ifvv,ifvv2,ifock
      double precision one,zilch
      character*8 ceval(2)
      data ceval    /'SCFEVLA0','SCFEVLB0'/
      parameter(one=1.d0)
      parameter(zilch=0.d0)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      iu=1
      iu2=1
      do 10 spin=1,uhf+1
         ifvv=1+nd2(1)*(spin-1)
         ifvv2=1+nd2(1)*(spin-1)
         nmo=nocco(spin)+nvrto(spin)
         i000=1
         i010=i000+nmo**2
         if(i010.gt.maxcor) call insmem('trnsfock',i010,maxcor)
         call zero(scr(i000),nmo**2)
         call getrec(20,'JOBARC',ceval(spin),nmo*iintfp,scr(i000))
         do ii=1,nmo-1
            scr(i000+ii+ii*nmo)=scr(i000+ii)
            scr(i000+ii)=0.0d0
         end do
         ifock=1+nocco(spin)
         do 20 irrep=1,nirrep
            virt=vrt(irrep,spin)
            drop=fnodrop(irrep,spin)
            vrtdrop=virt-drop
            i020=i010+virt**2
            i030=i020+virt**2
            if(i030.gt.maxcor)call insmem('trnsfock',i030,maxcor)            
            call blkcpy2(scr(i000),nmo,nmo,scr(i010),virt,virt,ifock,
     &         ifock)
            if (vrtdrop.gt.0) then
               call xgemm('t','n',vrtdrop,virt,virt,one,u(iu),virt,
     &            scr(i010),virt,zilch,scr(i020),vrtdrop)               
               call xgemm('n','n',vrtdrop,vrtdrop,virt,one,scr(i020),
     &            vrtdrop,u(iu),virt,zilch,fvv(ifvv),vrtdrop)
            endif
            if (drop.gt.0) then
               call xgemm('t','n',drop,virt,virt,one,u2(iu2),virt,
     &            scr(i010),virt,zilch,scr(i020),drop)
               call xgemm('n','n',drop,drop,virt,one,scr(i020),drop,
     &            u2(iu2),virt,zilch,fvv2(ifvv2),drop)
            endif
            ifock=ifock+virt
            iu=iu+virt*vrtdrop
            iu2=iu2+virt*drop
            ifvv=ifvv+vrtdrop**2
            ifvv2=ifvv2+drop**2
 20      continue
 10   continue
      return
1000  format((8(2x,I9)))
2000  format(t3,70('-'))
      end
