
c subroutines from the vdens program are used to calculate the relaxed 
c density matrix for the virtual virtual block
c
c D(a,b) is then diagonalized block by block to yield natural orbital
c occupation numbers and a set of natural orbitals in the mo basis which 
c are notated as:
c                     no 
c                  mo U
c
c U is truncated based on natural orbital occupation numbers
c
c                     no-d
c                 vrt U
c
c Next, the natural orbitals are expressed in the SO basis by the 
c following matrix multiplication:
c
c       vrt     no-d         no-d
c    so C   vrt U      =  so C'
c
c C is read from either the scfevca0 or scfevcb0 record.  
c
c The fock matrix is formed in the basis of C'
c
c         so      so    no-d           no-d
c    no-d C'~  so F  so C'     =  no-d F' 
c
c The reduced virtual virtual fock matrix is diagonalized
c
c       no-d        no-d        no-d          no-d
c  no-d Z~     no-d F'     no-d Z     =  no-d e'
c
c Z also contributes to C' 
c
c     no-d         no-d         no-d
c  so C'      no-d Z      =  so C'
c
c The following records are written to jobarc
c
c       no-d        no-d
c    so C'     no-d e'       
c
c    SCFEVCA0  SCFEVLA0  NUMDROPA MODROPA
c
c The dropmo variables are used to prevent the integral transformation 
c of the dropped MOs.  This involves the following records on jobarc:
c
c   NUMDROPA NUMDROPB  integer value for number of MOs to be dropped 
c   MODROPA  MODROPB  integer array of which MOs will be dropped
c
c it is worthy to note that the drop mo records are based on the ascending 
c ordering of orbital energies
c-------------------------------------------------------------------------
c
      program fno
c============variable declarations and common blocks======================
      implicit none
#include <flags.h>
#include <flags.com>
#include <flags2.com>
C     Commn blocks
      integer icore
      common//icore(1)
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
      logical rohf,qrhf,semi,natorb,anaderiv
      common/fnoparam/rohf,qrhf,semi,natorb,anaderiv
C     Local variables
      integer evenflag,ifvv,ifvv2,icp,icp2,ifnodrp,ilist,ncp(2),i000,
     &   i010,i020,i030,i040,i050,i060,iu,iu2,iocc,drop,nvlen,ifvvevl,
     &   ifvvevl2,nmo,irrep,uhf,maxcor,spin,inext,idvv,nd2len,i0,
     &   idrop,idrop2,ireord
      character*1 svar(2)

      data svar /'A','B'/
c------------------------------------------------------------------------
      call aces_init(icore,i0,maxcor,uhf,.true.)
c
c this part of the code uses common blocks and lists of dropmo dimensions
c      
      call setmet(uhf)

      anaderiv = (iFlags2(h_IFLAGS2_geom_opt).gt.0 .OR. 
     &            Iflags(h_IFLAGS_vib) .GT. 1)
      anaderiv = .TRUE.
      idvv=i0
      nd2len=nd2(1)+uhf*nd2(2)
      inext=idvv+nd2len*iintfp
      if (nd2len*iintfp.gt.maxcor)
     &   call insmem('densvv',nd2len*iintfp,maxcor)
      maxcor=maxcor-nd2len*iintfp
      call zero(icore(idvv),nd2len)
      call densvv(icore(idvv),icore(inext),maxcor,uhf)
      if (rohf) then

         if (semi) then
            idrop=inext
            inext=idrop+nd2(2)*iintfp
            if (inext-i0.gt.maxcor)
     &         call insmem('unsemi',nd2(2)*iintfp,maxcor)
            maxcor=maxcor-nd2(2)*iintfp
            evenflag=mod(maxcor,iintfp)
            maxcor=(maxcor-evenflag)/iintfp            
            call unsemi(icore(idvv),icore(idrop),icore(inext),maxcor)
            inext=idrop
            maxcor=maxcor*iintfp+evenflag+nd2(2)*iintfp
         endif
         call freezelem(icore(idvv),icore(idvv+nd2(1)*iintfp))
         idrop=inext
         call zero(icore(idrop),(nvrto(1)+nocco(1))**2)
      endif
      call getrec(20,'JOBARC','FNOFREEZ',1,drop)
      if (drop.gt.0) then
C     Handle DROPMO by explicity dropping virtual orbitals here
C      and allowing vtran to drop occupieds
         idrop=inext
         idrop2=idrop+drop
         ireord=idrop2+drop
         ifvvevl=ireord+nocco(1)+nvrto(1)
         inext=ifvvevl+(nocco(1)+nvrto(1))*iintfp
         if (inext-i0.gt.maxcor) call insmem('dropvrt',inext,maxcor)
         call dropvrt(icore(idvv),drop,icore(idrop),icore(idrop2),
     &      icore(ireord),icore(ifvvevl),nocco(1)+nvrto(1),uhf)
         inext=idrop
      endif

      iu=inext
      iu2=iu+nd2len*iintfp
      iocc=iu2+nd2len*iintfp
      nvlen=nvrto(1)+nvrto(2)*uhf
      ilist=iocc+nvlen*iintfp
      icp=ilist+nvrto(1)+mod(nvrto(1),2)
      ncp(1)=0
      ncp(2)=0
      do irrep=1,nirrep
         ncp(1)=ncp(1)+naobas(irrep)*vrt(irrep,1)
         ncp(2)=ncp(2)+naobas(irrep)*vrt(irrep,2)
      end do
      icp2=icp+(ncp(1)+uhf*ncp(2))*iintfp
      ifvv=icp2+(ncp(1)+uhf*ncp(2))*iintfp
      ifvv2=ifvv+nd2len*iintfp
      ifnodrp=ifvv2+nd2len*iintfp
      ifvvevl=ifnodrp+16
      ifvvevl2=ifvvevl+nvlen*iintfp
      inext=ifvvevl2+nvlen*iintfp
      maxcor=maxcor-(inext-iu)
      if (maxcor.lt.0) call insmem('FNO',inext-iu,maxcor)
      evenflag=mod(maxcor,iintfp)
      maxcor=(maxcor-evenflag)/iintfp
      call diagdvv(icore(idvv),icore(iu),icore(iocc),icore(inext),
     &      maxcor,uhf)

c break off here for non-hf natural orbitals 
      Write(6,*) "THE NATORB", NATORB
      if(natorb)then
         call natural(icore(iu),icore(icp),icore(inext),maxcor,uhf,
     &      ncp)
         maxcor=maxcor*iintfp+evenflag
         nmo=nocco(1)+nvrto(1)
         i000=inext
         i010=i000+nmo*iintfp
         i020=i010+nmo*iintfp
         i030=i020+nmo+mod(nmo,2)
         i040=i030+nmo+mod(nmo,2)
         i050=i040+nmo**2*iintfp
         i060=i050+nmo**2*iintfp
         if((i060-i000).gt.maxcor)
     &      call insmem('dmpnat',i060-i000,maxcor)
         maxcor=maxcor-(i060-i000)
         evenflag=mod(maxcor,iintfp)
         maxcor=(maxcor-evenflag)/iintfp
         call dmpnat(icore(icp),icore(i000),icore(i010),icore(i020),
     &      icore(i030),icore(i040),icore(i050),icore(i060),nmo,maxcor,
     &      uhf,ncp)
         call aces_fin
         stop       
      endif

      maxcor=maxcor*iintfp+evenflag
      if (maxcor.lt.16) call insmem('FNO',16,maxcor)
      call redvrt(icore(iocc),icore(ilist),icore(ifnodrp),icore(inext),
     &   uhf,drop)
      
      evenflag=mod(maxcor,iintfp)
      maxcor=(maxcor-evenflag)/iintfp
      if (anaderiv) then
         call putrec(20,'JOBARC','FNODROP ',16,icore(ifnodrp))
      endif
      call shku(icore(iu),icore(iu2),icore(ifnodrp),icore(inext),maxcor,
     &   uhf)
      call no2so(icore(iu),icore(iu2),icore(icp),icore(icp2),
     &   icore(ifnodrp),icore(inext),maxcor,uhf,ncp)
      if (rohf) then
         call trnsfockrohf(icore(icp),icore(icp2),icore(ifvv),
     &      icore(ifvv2),icore(ifnodrp),icore(inext),maxcor,ncp)
      else
         call trnsfock(icore(iu),icore(iu2),icore(ifvv),icore(ifvv2),
     &      icore(ifnodrp),icore(inext),maxcor,uhf)
      endif
      if (semi) then
         nmo=nocco(1)+nvrto(1)
         i000=inext
         i010=i000+nmo
         i020=i010+nmo
         i030=i020+nmo*naobasis
         i040=i030+nmo*naobasis
         if((i040-i000).gt.maxcor) 
     &      call insmem('dmprohf',i050-i000,maxcor)
         maxcor=maxcor-(i040-i000)
C     Dump data for semicanonical ROHF calcs
         call dmprohf(icore(icp),icore(ifvv),icore(ifnodrp),icore(i000),
     &      icore(i010),icore(i020),icore(i030),nmo,ncp)
      endif
      call diagfvv(icore(ifvv),icore(ifvv2),icore(ifvvevl),
     &   icore(ifvvevl2),icore(ifnodrp),icore(icp),icore(icp2),
     &   icore(inext),maxcor,uhf,ncp)
      if (anaderiv) then
C     For FNO gradients need the entire fock matrix
         call fullfock(icore(icp),icore(icp2),icore(ifvv),
     &      icore(ifnodrp),icore(inext),maxcor,uhf,ncp)
         maxcor=maxcor*iintfp+evenflag
         nmo=nocco(1)+nvrto(1)
         i000=inext
         i010=i000+nmo*iintfp
         i020=i010+nmo*iintfp         
         i030=i020+nmo*naobasis*iintfp
         i040=i030+nmo*naobasis*iintfp
         i050=i040+nmo
         if((i050-i000).gt.maxcor) 
     &      call insmem('dmpgrad',i050-i000,maxcor)
         maxcor=maxcor-(i050-i000)
         evenflag=mod(maxcor,iintfp)
         maxcor=(maxcor-evenflag)/iintfp
C     Dump data for FNO gradients
         call dmpgrad(icore(icp),icore(icp2),icore(ifvvevl),
     &      icore(ifvvevl2),icore(ifnodrp),icore(i000),icore(i010),
     &      icore(i020),icore(i030),icore(i040),nmo,uhf,ncp)
         call aces_fin
         stop
      endif      
      maxcor=maxcor*iintfp+evenflag
c-----------------------------------------------------------------------------
C       dump reduced virtual space
c
      nmo=nocco(1)+nvrto(1)
      i000=inext
      i010=i000+nmo*iintfp
      i020=i010+nmo*iintfp
      i030=i020+nmo*naobasis*iintfp
      i040=i030+nmo*naobasis*iintfp
      i050=i040+nmo
      if((i050-i000).gt.maxcor) call insmem('dmprecs',i050-i000,maxcor)
      maxcor=maxcor-(i050-i000)
      evenflag=mod(maxcor,iintfp)
      maxcor=(maxcor-evenflag)/iintfp
      Write(6,*) "Engtering DMPRECS"
      call dmprecs(icore(icp),icore(ifvvevl),icore(ifnodrp),icore(i000),
     &   icore(i010),icore(i020),icore(i030),icore(i040),nmo,uhf,ncp)
      call aces_fin
      stop
      end
