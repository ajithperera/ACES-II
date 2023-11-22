      subroutine dmprohf(cp,fvv,fnodrop,scfevl,scfevl2,scfevc,scfevc2,
     &   nmo,ncp)
c
c this routine writes the following to jobarc for semicanonical ROHF
C     ROHFEVC and ROHFEVL the PRE-semicanonicalized ROHF orbs
C
c
c
c================variable declarations and common blocks======================
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
      integer nmo,fnodrop(8,2),ncp(2)
      double precision cp(*),fvv(*)
C     Pre-allocated Local Variables
      double precision scfevl(nmo),scfevl2(nmo),scfevc(naobasis,nmo),
     &   scfevc2(naobasis,nmo)
C     Local Variables
      integer nmodrop,vrtdrop,orbs,ii,occ,spin,ifvv,icp,aos,iscfevl,
     &   irrep,iscfevl2,iscfcol,iscfrow2,iscfcol2
      character*8 cscfevl(2),cscfevc(2)
      data cscfevl /'ROHFEVLA','ROHFEVLB'/
      data cscfevc /'ROHFEVCA','ROHFEVCB'/
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ifvv=1
      do 10 spin=1,2
         icp=1+ncp(1)*(spin-1)
c----------------write the new orbital energies to jobarc---------------------
         call getrec(20,'JOBARC',cscfevl(spin),nmo*iintfp,scfevl)
         call zero(scfevc2,nmo*naobasis)
         call zero(scfevl2,nmo)
         iscfevl=1
         iscfevl2=1
         nmodrop=nmo
         do irrep=1,nirrep
            nmodrop=nmodrop-fnodrop(irrep,spin)
            occ=pop(irrep,spin)
            vrtdrop=vrt(irrep,spin)-fnodrop(irrep,spin)
            call dcopy(occ,scfevl(iscfevl),1,scfevl2(iscfevl2),1)
            iscfevl=iscfevl+occ
            iscfevl2=iscfevl2+occ
            call dcopy(vrtdrop,fvv(ifvv),vrtdrop+1,scfevl2(iscfevl2),1)
            ifvv=ifvv+vrtdrop**2
            iscfevl2=iscfevl2+vrtdrop
         end do
         call putrec(20,'JOBARC',cscfevl(spin),nmodrop*iintfp,scfevl2)
c--------------write the new orbital coefficents to jobarc--------------------
         call getrec(20,'JOBARC',cscfevc(spin),nmo*naobasis*iintfp,
     &      scfevc)
         iscfcol=1
         iscfcol2=1
         iscfrow2=1
         do irrep=1,nirrep
            aos=naobas(irrep)
            occ=pop(irrep,spin)
            vrtdrop=vrt(irrep,spin)-fnodrop(irrep,spin)
            call dcopy(occ*naobasis,scfevc(1,iscfcol),1,
     &         scfevc2(1,iscfcol2),1)
            iscfcol=iscfcol+occ
            iscfcol2=iscfcol2+occ
            do orbs=1,vrtdrop
               call dcopy(aos,cp(icp),1,scfevc2(iscfrow2,iscfcol2),1)
               icp=icp+aos
               iscfcol2=iscfcol2+1
            end do
            iscfrow2=iscfrow2+aos
         end do

         call putrec(20,'JOBARC',cscfevc(spin),naobasis*nmodrop*iintfp,
     &      scfevc2)
 10   continue

      return
      end
