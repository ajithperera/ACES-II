      subroutine unsemi(dvv,bcktrns,scr,maxcor)
      implicit none
C     Back transform density matrix from semi-canonical basis
C
C     Common Blocks
      integer iintln,ifltln,iintfp
      common/machsp/iintln,ifltln,iintfp
      integer pop(8,2),vrt(8,2),nt(2),nd1(2),nd2(2)
      common/sym/pop,vrt,nt,nd1,nd2
      integer nstart,nirrep
      common/syminf/nstart,nirrep
C     Input Variables
      integer maxcor
C     Input/Output Variables
      double precision dvv(*)
C     Pre-allocated Local Variables
      double precision bcktrns(*),scr(maxcor)
C     Local Variables
      integer spin,irrep,idvv,ibck,i000,i010,virt
      double precision one,zilch
      character*8 wlabel(2)
      data wlabel /'WMATRXVA','WMATRXVB'/
      parameter(one=1.d0)
      parameter(zilch=0.d0)
C----------------------------------------------------------------
      idvv=1
      i000=1
      do 1 spin=1,2
         ibck=1
         call zero(bcktrns,nd2(spin))
         call getrec(20,'JOBARC',wlabel(spin),nd2(spin)*iintfp,
     &      bcktrns)
         do 11 irrep=1,nirrep
            virt=vrt(irrep,spin)
            i010=i000+virt**2
            if(i010.gt.maxcor)
     &      call insmem('unsemi',i010,maxcor)
            call xgemm('n','n',virt,virt,virt,one,bcktrns(ibck),virt,
     &         dvv(idvv),virt,zilch,scr(i000),virt)
            call xgemm('n','t',virt,virt,virt,one,scr(i000),virt,
     &         bcktrns(ibck),virt,zilch,dvv(idvv),virt)
            idvv=idvv+virt**2
            ibck=ibck+virt**2
 11      continue
 1    continue

      return
      end
            
            
