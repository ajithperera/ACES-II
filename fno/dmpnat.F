      subroutine dmpnat(cp,scfevl,scfevl2,ireorder,iscr,scfevc,scfevc2,
     &   scr,nmo,maxcor,uhf,ncp)
c combine orbitals 
c    occupied are hartree fock
c    virtual are natural orbitals 
c
c new orbital energies are diagonal fock matrix elements
c these appear to need to be sorted
c          so     so    hf/no          hf/no
c    hf/no Z~  so F  so Z      = hf/no f
c
c this routine writes the following to jobarc
c       new molecular orbitals and energies
c       iflags(2) calculation level for natural orbitals 
c 	turn non-hf flag on
c       turn property flag off
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
C     Input Variables
      integer uhf,maxcor,nmo,ncp(2)
      double precision cp(*)
C     Pre-allocated Local Variables
      integer ireorder(nmo),iscr(nmo)
      double precision scfevl(nmo),scfevl2(nmo),scfevc(nmo,nmo),
     &   scfevc2(nmo,nmo),scr(maxcor)
C     Local variables
      double precision tolerance,one,zilch
      integer iscf,ifock,occd,j,ievl,i000,i010,i020,i030,ivrt,iscfr,
     &   iscfc,icp,spin,i,iocc,mos,irrep,virt
      character*8 cscfevl(2),cscfevc(2),cfock(2)
      data cscfevl  /'SCFEVLA0','SCFEVLB0'/
      data cscfevc  /'SCFEVCA0','SCFEVCB0'/
      data cfock    /'FOCKA   ','FOCKB   '/
      parameter(tolerance=1.0d-9)
      parameter(one=1.d0)
      parameter(zilch=0.d0)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do 10 spin=1,uhf+1
         icp=1+ncp(1)*(spin-1)
         call getrec(20,'JOBARC',cscfevc(spin),nmo**2*iintfp,scfevc)
         iscfr=1
         iscfc=nocco(spin)+1
         do 20 irrep=1,nirrep
            virt=vrt(irrep,spin)
            mos=pop(irrep,spin)+virt
            call blkcpy(cp(icp),mos,virt,scfevc,nmo,nmo,iscfr,iscfc)
            icp=icp+mos*virt
            iscfr=iscfr+mos
            iscfc=iscfc+virt
 20      continue
c fill ireorder and convert from post vtran to prevtran ordering
         i=1
         iocc=1
         ivrt=nocco(spin)+1
         do 40 irrep=1,nirrep
            occd=pop(irrep,spin)
            virt=vrt(irrep,spin)
            do 41 j=iocc,iocc+occd-1
               ireorder(i)=j
               i=i+1
 41         continue
            do 42 j=ivrt,ivrt+virt-1
               ireorder(i)=j
               i=i+1
 42         continue
            iocc=iocc+occd
            ivrt=ivrt+virt
 40      continue
         do 50 i=1,nmo
            do 51 j=1,nmo
               scfevc2(j,i)=scfevc(j,ireorder(i))
 51         continue
 50      continue
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c do matrix multiplication
c          so     so    hf/no          hf/no
c    hf/no Z~  so F  so Z      = hf/no f

         call getrec(20,'JOBARC',cfock(spin),nmo**2*iintfp,scfevc)
         iscf=1
         ifock=1
         ievl=1
         i000=1
         do 60 irrep=1,nirrep
            mos=pop(irrep,spin)+vrt(irrep,spin)
            i010=i000+mos**2
            i020=i010+mos**2
            i030=i020+mos**2
            if(i030.gt.maxcor) call insmem('dmpnat',i030,maxcor)  
c copy orbitals
            call blkcpy2(scfevc2,nmo,nmo,scr(i000),mos,mos,iscf,iscf)
c copy so based fock matrix
            call blkcpy2(scfevc,nmo,nmo,scr(i010),mos,mos,ifock,ifock)
            
            call xgemm('t','n',mos,mos,mos,one,scr(i000),mos,scr(i010),
     &         mos,zilch,scr(i020),mos)

            call xgemm('n','n',mos,mos,mos,one,scr(i020),mos,scr(i000),
     &         mos,zilch,scr(i010),mos)

            call scopy(mos,scr(i010),mos+1,scr(i020),1)
 
c sort orbital energies
            do 65 i=1,mos
               iscr(i)=i
 65         continue 
            call piksr2(mos,scr(i020),iscr)
c write diagonal fock matrix elements to scfevl
            call scopy(mos,scr(i020),1,scfevl(ievl),1)

            write(6,*)' diagonal fock matrix elements'
            write(6,1000)(scfevl(i),i=ievl,ievl+mos-1)
 1000       format((5(2x,f13.6)))

c reorder orbitals
            do 66 i=1,mos
               call scopy(mos,scr(i000+(i-1)*mos),1,
     &            scr(i010+(iscr(i)-1)*mos),1)
 66         continue 

            call blkcpy(scr(i000),mos,mos,scfevc2,nmo,nmo,iscf,iscf) 

            iscf=iscf+mos
            ifock=ifock+mos
            ievl=ievl+mos
 60      continue
         call putrec(20,'JOBARC',cscfevl(spin),nmo*iintfp,scfevl)
         call putrec(20,'JOBARC',cscfevc(spin),nmo**2*iintfp,scfevc2)
 10   continue
c------------------------change the flags------------------------------
c turn non-hf on
      iflags(h_IFLAGS_nonhf)=1

      call putrec(20,'JOBARC','IFLAGS  ',100,iflags)
      return
      end


