      subroutine fullfock(cp,cp2,fock,fnodrop,scr,maxcor,uhf,ncp)
c
c         ao      ao    d              d
c    no-d U'~  ao F  ao U'     =  no-d F'
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
      integer naobasis,naobas(8)
      common/aoinfo/naobasis,naobas
C     Input Variables
      integer maxcor,uhf,fnodrop(8,2),ncp(2)
      double precision cp(*),cp2(*)
C     Output Variables
      double precision fock(*)
C     Pre-allocated Local Variables
      double precision scr(maxcor)
C     Local variables
      integer virt,drop,vrtdrop,spin,i000,i010,i020,i030,irrep,nao,icp,
     &   icp2,iloc,ifock
      double precision one,zilch
      character*8 cfock(2)
      data cfock    /'FOCKA   ','FOCKB   '/
      parameter(one=1.d0)
      parameter(zilch=0.d0)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      icp=1
      icp2=1
      do 10 spin=1,uhf+1
         ifock=1+nd2(1)*(spin-1)
         i000=1
         i010=i000+naobasis**2
         if(i010.gt.maxcor) call insmem('trnsfock',i010,maxcor)
         call zero(scr(i000),naobasis**2)
         call getrec(20,'JOBARC',cfock(spin),naobasis**2*iintfp,
     &      scr(i000))
         iloc=1
         do 20 irrep=1,nirrep
            nao=naobas(irrep)
            virt=vrt(irrep,spin)
            drop=fnodrop(irrep,spin)
            vrtdrop=virt-drop
            if ((vrtdrop.gt.0).and.(drop.gt.0)) then
               i020=i010+nao**2
               i030=i020+nao*vrtdrop
               if(i030.gt.maxcor)call insmem('trnsfock',i030,maxcor)
               call blkcpy2(scr(i000),naobasis,naobasis,scr(i010),nao,
     &            nao,iloc,iloc)
               call xgemm('t','n',vrtdrop,nao,nao,one,cp(icp),nao,
     &            scr(i010),nao,zilch,scr(i020),vrtdrop)
               call xgemm('n','n',vrtdrop,drop,nao,one,scr(i020),
     &            vrtdrop,cp2(icp2),nao,zilch,fock(ifock),vrtdrop)
            endif
            ifock=ifock+vrtdrop*drop
            icp=icp+nao*vrtdrop
            icp2=icp2+nao*drop
            iloc=iloc+nao
 20      continue
 10   continue
      return
1000  format((8(2x,I9)))
2000  format(t3,70('-'))
      end
