      subroutine dmpgrad(cp,cp2,fvvevl,fvvevl2,fnodrop,scfevl,scfevl2,
     &   scfevc,scfevc2,orbirrep,nmo,uhf,ncp)
c
c this routine writes the following to jobarc
c	new molecular orbitals and energies 
c 	dropmo records consisting of numdropa/numdropb and fnodropa/fnodropb
c       iflags(2) calculation level in reduced space
c
c
c================variable declarations and common blocks=====================
      implicit none
#include "flags.h"
C     Common blocks
      integer nocco(2),nvrto(2)
      common/info/nocco,nvrto
      integer iflags(100)
      common/flags/iflags
      integer iintln,ifltln,iintfp
      common/machsp/iintln,ifltln,iintfp
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      integer pop(8,2),vrt(8,2),nt(2),nd1(2),nd2(2)
      common/sym/pop,vrt,nt,nd1,nd2
      integer naobasis,naobas(8)
      common/aoinfo/naobasis,naobas
C     Input Variables
      integer nmo,uhf,fnodrop(8,2),ncp(2)
      double precision cp(*),cp2(*),fvvevl(*),fvvevl2(*)
C     Pre-allocated Local Variables
      integer orbirrep(nmo)
      double precision scfevl(nmo),scfevl2(nmo),scfevc(naobasis,nmo),
     &   scfevc2(naobasis,nmo)
C     Local Variables
      integer nmodrop(2),vrtdrop(nirrep),orbs,ii,occ,spin,ifvvevl,icp,
     &   iscfevl,irrep,iorb,iscfevl2,nvrttmp,iscfcol,iscfrow2,iscfcol2,
     &   iunit,icp2,ifvvevl2,mos,dropocc,droplen
      double precision mbptenergy
      character*8 cscfevl(2),cscfevc(2),symvrt(2),irrepl(2)
      data cscfevl  /'SCFEVLA0','SCFEVLB0'/
      data cscfevc  /'SCFEVCA0','SCFEVCB0'/
      data symvrt   /'SYMPOPVA','SYMPOPVB'/
      data irrepl   /'IRREPALP','IRREPBET'/
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      nvrttmp=nvrto(1)
      do 10 spin=1,uhf+1
         icp=1+ncp(1)*(spin-1)
         icp2=1+ncp(1)*(spin-1)
         ifvvevl=1+nvrttmp*(spin-1)
         ifvvevl2=1+nvrttmp*(spin-1)
c----------------write the new orbital energies to jobarc---------------------
         call getrec(20,'JOBARC',cscfevl(spin),nmo*iintfp,scfevl)
         call zero(scfevc2,nmo*naobasis)
         call zero(scfevl2,nmo)
         iscfevl=1
         iscfevl2=1
         nmodrop(spin)=nmo
         iorb=1
         do irrep=1,nirrep
            nmodrop(spin)=nmodrop(spin)-fnodrop(irrep,spin)
            occ=pop(irrep,spin)
            vrtdrop(irrep)=vrt(irrep,spin)-fnodrop(irrep,spin)
            do ii=1,occ+vrt(irrep,spin)
               orbirrep(iorb+ii-1)=irrep
            end do
            iorb=iorb+occ+vrt(irrep,spin)
            call dcopy(occ,scfevl(iscfevl),1,scfevl2(iscfevl2),1)
            iscfevl=iscfevl+occ
            iscfevl2=iscfevl2+occ
            call dcopy(vrtdrop(irrep),fvvevl(ifvvevl),1,
     &         scfevl2(iscfevl2),1)
            ifvvevl=ifvvevl+vrtdrop(irrep)
            iscfevl2=iscfevl2+vrtdrop(irrep)
            do ii=1,fnodrop(irrep,spin)
               scfevl2(iscfevl2)=fvvevl2(ifvvevl2)+999.0d0
               iscfevl2=iscfevl2+1
               ifvvevl2=ifvvevl2+1
            end do
         end do
         call putrec(20,'JOBARC',irrepl(spin),nmo,orbirrep)
         call putrec(20,'JOBARC',symvrt(spin),nirrep,vrt(1,spin))
         call putrec(20,'JOBARC',cscfevl(spin),nmo*iintfp,scfevl2)
c--------------write the new orbital coefficents to jobarc--------------------
         call getrec(20,'JOBARC',cscfevc(spin),nmo*naobasis*iintfp,
     &      scfevc)
         iscfcol=1
         iscfcol2=1
         iscfrow2=1
         do irrep=1,nirrep
            occ=pop(irrep,spin)
            mos=naobas(irrep)
            call dcopy(occ*naobasis,scfevc(1,iscfcol),1,
     &         scfevc2(1,iscfcol2),1)
            iscfcol=iscfcol+occ
            iscfcol2=iscfcol2+occ
            do orbs=1,vrtdrop(irrep)
               call dcopy(mos,cp(icp),1,scfevc2(iscfrow2,iscfcol2),1)
               icp=icp+mos
               iscfcol2=iscfcol2+1
            end do
            do orbs=1,fnodrop(irrep,spin)
               call dcopy(mos,cp2(icp2),1,scfevc2(iscfrow2,iscfcol2),1)
               icp2=icp2+mos
               iscfcol2=iscfcol2+1
            end do
            iscfrow2=iscfrow2+mos
         end do

         call putrec(20,'JOBARC',cscfevc(spin),naobasis*nmo*iintfp,
     &      scfevc2)
 10   continue
      if (uhf.eq.0) nmodrop(2)=nmodrop(1)
      call getrec(20,'JOBARC','FNOFREEZ',1,dropocc)
      call putrec(20,'JOBARC','NVRTORB',2,nvrto)
      call putrec(20,'JOBARC','NUMDROPA',1,nmo-nmodrop(1)+dropocc)
      call getrec(0,'JOBARC','MODROPA ',droplen,orbirrep)
      if (droplen.gt.0)
     &   call getrec(20,'JOBARC','MODROPA ',droplen,orbirrep)
      do ii=nmodrop(1)+1,nmo
         orbirrep(ii-nmodrop(1)+dropocc)=ii
      end do
      call putrec(20,'JOBARC','MODROPA ',nmo-nmodrop(1)+dropocc,
     &   orbirrep)
      if (uhf.ne.0) then
         call getrec(0,'JOBARC','MODROPB ',droplen,orbirrep)
         if (droplen.gt.0)
     &      call getrec(20,'JOBARC','MODROPB ',droplen,orbirrep)
         do ii=nmodrop(2)+1,nmo
            orbirrep(ii-nmodrop(2)+dropocc)=ii
         end do
         call putrec(20,'JOBARC','NUMDROPB',1,nmo-nmodrop(2)+dropocc)         
         call putrec(20,'JOBARC','MODROPB ',nmo-nmodrop(2)+dropocc,
     &      orbirrep)
      endif
      
      if (uhf.eq.0) then
         nt(1)=0
         nd2(1)=0
         do irrep=1,nirrep
            vrt(irrep,2)=vrt(irrep,1)
            nt(1)=nt(1)+pop(irrep,1)*vrt(irrep,1)
            nd2(1)=nd2(1)+vrt(irrep,1)**2
         end do
         nd2(2)=nd2(1)
         nt(2)=nt(1)
      else
         nt(1)=0
         nt(2)=0
         nd2(1)=0
         nd2(2)=0
         do irrep=1,nirrep
            nt(1)=nt(1)+pop(irrep,1)*vrt(irrep,1)
            nt(2)=nt(2)+pop(irrep,2)*vrt(irrep,2)
            nd2(1)=nd2(1)+vrt(irrep,1)**2
            nd2(2)=nd2(2)+vrt(irrep,2)**2
         end do
      endif
      call putrec(20,'JOBARC','NDROPGEO',1,1)
      iflags(h_IFLAGS_nonhf)=1
      call putrec(20,'JOBARC','IFLAGS  ',100,iflags)
      call getrec(20,'JOBARC','TOTENERG',iintfp,mbptenergy)
      call putrec(20,'JOBARC','PREFNOEN',iintfp,mbptenergy)

      return
 666  call errex
2000  format(t3,70('-'))
2001  format((10(2x,I4)))
      end
