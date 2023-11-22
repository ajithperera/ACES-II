        subroutine vhfksdint(icore,icrsiz,iuhf)

       implicit none

#include "aces.h"
#include "flags.h"
#include "flags.com"
#include "flags2.com"
#include "machsp.com" 

       integer icore(icrsiz),icrsiz,icrsiz2,iuhf
       integer ntpert,npert,ipert,ixpert,iypert,izpert,ixypert,
     &        ixzpert,iyzpert,itransx,itransy,itransz,
     &        nucind
      common/pert/ntpert,npert(8),ipert(8),ixpert,iypert,izpert,
     &            ixypert,ixzpert,iyzpert,itransx,itransy,itransz,
     &            nucind 
      logical c1
      integer pop,vrt,nt,nfmi,nfea
      common /sym/ pop(8,2),vrt(8,2),nt(2),nfmi(2),nfea(2)
      integer nstart,nirrep,irrepa,dirprd
      common/syminf/nstart,nirrep,irrepa(255,2),dirprd(8,8)
      integer irpdpd,isytyp,njunk
      common/sympop/irpdpd(8,22),isytyp(2,500),njunk(18) 
      integer istart
      integer ifiaaa,iend
      integer irep,ip,np,ixaia,ixaib
      integer ifijaa,ixija,ixijb,ieval
      double precision gradient(9),ddot
      integer ifiabb,natom,ntpert1
      logical trainv
      integer grad,igrad,istR,ioffseT,cgrad,idfaia
      istart=1
      call setbashfks
      call getrec(20,'JOBARC','NSYMPERT',27,ntpert)
      call getrec(20,'JOBARC','NREALATM',iintfp,natom)
      call formsd(icore(istart),icrsiz,iuhf)       
      call formfd(icore(istart),icrsiz,iuhf)
       
      ifiaaa=istart
      ifiabb=ifiaaa+nt(1)*iintfp
      ixaia =ifiabb+nt(2)*iintfp
      ixaib =ixaia +nt(1)*iintfp
      iend  =ixaib +nt(2)*iintfp
      icrsiz2=icrsiz-2*(nt(1)+nt(2))*iintfp

      
      call getrec(1,'JOBARC','HFKSFIAA',nt(1)*iintfp,
     &            icore(iFiaaa))
      call getrec(-1,'JOBARC','HFKSFIAB',nt(2)*iintfp,
     &            icore(ifiabb))

      call makdia1(icore(iFiaaa),
     &             icore(iend),icrsiz2,iuhf,.false.) 

#ifdef _DEBUG_LVL0
       call ppp(8,npert) 
#endif
       trainv=.false.
       ntpert1=ntpert
       if(iflags(73).eq.0 .and. nirrep .eq. 1 ) then
          trainv=.true.
          ntpert1=ntpert+3
        end if

       igrad=iend
       istr=igrad+ntpert*iintfp
       cgrad=istr+ntpert1*ntpert1*iintfp
       grad=cgrad+ntpert1*iintfp

      call dzero(icore(igrad),ntpert)
      call dzero(icore(cgrad),ntpert1)
        
      ioffset=1
       nirrep=1
 
      do 100 irep=1,nirrep
         np=npert(irep)
         if(np.eq.0) goto 100
         do 150 ip=1,np
           call getlst(icore(iXaia),ip,1,1,irep,280)
           if(iuhf .eq.0 ) then
           call calcgrad(icore(ixaia),icore(ifiaaa),
     &         icore(igrad),nt(1),ntpert,ioffset)
           else                
            call getlst(icore(iXaib),ip,1,1,irep,281)
            call calcgradBB(icore(ixaia),icore(ixaib),icore(ifiaaa),
     &       icore(ifiabb),icore(igrad),nt(1),nt(2),ntpert,ioffset)
           end if
          ioffset=ioffset+1
 150     continue
 100     continue
        if(.not. trainv) then
       call getrec(20,'JOBARC','SYMCOORD',iintfp*ntpert*ntpert,
     &              icore(istr))
CSSS       write(*,*) 'cartesian coordinate'
CSSS       call kkk(ntpert,icore(igrad))
       call traSGC(icore(igrad),icore(cgrad),icore(istr),
     &            ntpert)
       else
        call gen_mis_grad(icore(igrad),icore(cgrad),natom,ntpert1)
        call dcopy(ntpert,icore(igrad),1,icore(cgrad+3),1)
       end if
CSSS        write(*,*) 'symmtery coordinate'
CSSS       call kkk(ntpert1,icore(cgrad))
       call getrec(1,'JOBARC','GRADIENT',iintfp*ntpert1,
     &          icore(grad))
       call daxpy(ntpert1,1.0d0,icore(grad),1,icore(cgrad),1)
       call putrec(20,'JOBARC','GRADIENT',ntpert1*iintfp,
     &            icore(cgrad))
       return
       end 
 
