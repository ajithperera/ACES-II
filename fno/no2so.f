      subroutine no2so(u,u2,cp,cp2,fnodrop,scr,maxcor,uhf,ncp)
c
c       vrt     no-d        no-d
c    so C   vrt U      =  so C'
c
c================variable declarations and common blocks=====================
      implicit none
C     Common blocks
      integer nocco(2),nvrto(2)
      common/info/nocco,nvrto
      integer iintln,ifltln,iintfp
      common/machsp/iintln,ifltln,iintfp
      integer nstart,nirrep,irreps(255,2),dirprd(8,8)      
      common/syminf/nstart,nirrep,irreps,dirprd
      integer pop(8,2),vrt(8,2),nt(2),nd1(2),nd2(2)
      common/sym/pop,vrt,nt,nd1,nd2
      integer naobasis,naobas(8)
      common/aoinfo/naobasis,naobas
      logical rohf,qrhf,semi
      common/fnoparam/rohf,qrhf,semi
C     Input variables
      integer fnodrop(8,2),uhf,maxcor,ncp(2)
      double precision u(*),u2(*)
C     Input/Output variables
      double precision cp(*),cp2(*)
C     Pre-allocated Local variables
      double precision scr(maxcor)
C     Local variables
      integer spin,iu,iu2,icp,icp2,nmo,i000,i010,i020,iscfc,iscfr,irrep,
     &   virt,drop,vrtdrop,aos,ii,jj,basis2
      double precision one,zilch
      character*8 cscfevc(2)
      parameter(one=1.d0)
      parameter(zilch=0.d0)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      iu=1
      iu2=1
      if (semi) then
         cscfevc(1)='ROHFEVCA'
         cscfevc(2)='ROHFEVCB'    
      else
         cscfevc(1)='SCFEVCA0'
         cscfevc(2)='SCFEVCB0'
      endif
      do 10 spin=1,uhf+1
         icp=1+ncp(1)*(spin-1)
         icp2=1+ncp(1)*(spin-1)
         nmo=nocco(spin)+nvrto(spin)
         basis2=nmo*naobasis
         i000=1
         i010=i000+basis2
         if(i010.gt.maxcor)
     &      call insmem('no2so',i010,maxcor)
         call zero(scr(i000),basis2)
         call getrec(20,'JOBARC',cscfevc(spin),basis2*iintfp,scr(i000))
         iscfc=nocco(spin)+1
         iscfr=1
         do 20 irrep=1,nirrep
            virt=vrt(irrep,spin)
            drop=fnodrop(irrep,spin)
            vrtdrop=virt-drop
            aos=naobas(irrep)
            i020=i010+virt*aos
            if (i020.gt.maxcor)
     &         call insmem('no2so',i020,maxcor)
            call zero(scr(i010),virt*aos)
            call blkcpy2(scr(i000),naobasis,nmo,scr(i010),aos,virt,
     &         iscfr,iscfc)
            if (vrtdrop.gt.0) call xgemm('n','n',aos,vrtdrop,virt,one,
     &         scr(i010),aos,u(iu),virt,zilch,cp(icp),aos)
            if (drop.gt.0) call xgemm('n','n',aos,drop,virt,one,
     &         scr(i010),aos,u2(iu2),virt,zilch,cp2(icp2),aos)
            iscfr=iscfr+aos
            iscfc=iscfc+virt
            iu=iu+virt*vrtdrop
            iu2=iu2+virt*drop
            icp=icp+aos*vrtdrop
            icp2=icp2+aos*drop
 20      continue
 10   continue
      return 
      end
