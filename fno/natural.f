      subroutine natural(u,cp,scr,maxcor,uhf,ncp)
c
c       vrt     no           no
c    so C   vrt U      =  so C'
c
c================variable declarations and common blocks=====================
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
      integer uhf,maxcor,ncp(2)
      double precision u(*)
C     Output Variables
      double precision cp(*)
C     Pre-allocated Local Variables
      double precision scr(maxcor)
C     Local variables      
      double precision one,zilch
      integer iscfr,iscfc,i000,i010,i020,mos,virt,irrep,icp,iu,nmo,spin
      character*8 cscfevc(2)
      data cscfevc  /'SCFEVCA0','SCFEVCB0'/
      parameter(one=1.d0)
      parameter(zilch=0.d0)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      iu=1
      i000=1
      do 10 spin=1,uhf+1
	 icp=1+ncp(1)*(spin-1)
         nmo=nocco(spin)+nvrto(spin)
         i010=i000+nmo**2
         i020=i010+nvrto(spin)**2
         if(i020.gt.maxcor) call insmem('natural',i020,maxcor)
         call getrec(20,'JOBARC',cscfevc(spin),nmo**2*iintfp,scr(i000))
         iscfc=nocco(spin)+1 
         iscfr=1
         do 20 irrep=1,nirrep
            virt=vrt(irrep,spin)
            mos=pop(irrep,spin)+virt
            call blkcpy2(scr(i000),nmo,nmo,scr(i010),mos,virt,iscfr,
     &         iscfc)
            call xgemm('n','n',mos,virt,virt,one,scr(i010),mos,u(iu),
     &         virt,zilch,cp(icp),mos)
            iscfr=iscfr+mos
            iscfc=iscfc+virt
            iu=iu+virt**2
            icp=icp+mos*virt
 20      continue
 10   continue
      return 
      end

