      subroutine gaussprak(flag,grid,type,grdangpts,gridxyz,gridwt,
     &         ntheta,nphi,glwght,glroot)

      implicit none
#include<aces.h>
#include<grid.com>

      integer grid,type,grdangpts(numgrid),flag

      M_REAL
     &    gridxyz(3,maxangpts,numgrid),gridwt(maxangpts,numgrid)

      integer maxlebgridP

      integer ntheta,nphi,itheta,iphi,ifact,ithea,itheb,iphia,iphib,
     &    i,nexppts(type)

      M_REAL
     &  pi,wthe,cost,sint,dphi,phii,wphi,cosp,sinp
 
      M_REAL
     &  glwght(ntheta,ntheta),glroot(ntheta,ntheta)

      M_REAL
     &    one,two,half
 
       parameter(one=1.0d0,two=2.0d0,half=0.5d0)

       call callstack_push('GAUSSPRAK')

      PI=ATAN(1.0D+00)*4.0D+00

      nexppts(type)=ntheta*nphi
      if(flag.eq.0) then
       ntotrad=ntotrad+nexppts(type)*gridlist(grid,3)
       if (nexppts(type).gt.maxangpts) then
          maxangpts=nexppts(type)
          maxanggrd=grid
      endif
        goto 999
      endif
      grdangpts(grid)=nexppts(type)

      do itheta=1,ntheta
      call findwt(-one,one,glroot(1,itheta),glwght(1,itheta),itheta)
      end do
  
      ntheta=max(ntheta,2)
      ifact=1
      ithea=1
      itheb=ntheta
      iphia=1
      iphib=nphi

      do itheta=ithea,itheb
        wthe=glwght(itheta,ntheta)
        cost=glroot(itheta,ntheta)
        sint=SQRT(one-cost**2)
        do iphi=iphia,iphib
           i=(itheta-ithea)*(iphib-iphia+1)+iphi
           dphi=two*pi/nphi
           phii=dphi*(iphi-half)
           wphi=dphi
           cosp=cos(phii)
           sinp=sin(phii)
           gridxyz(1,i,numgrid)=sint*cosp
           gridxyz(2,i,numgrid)=sint*sinp
           gridxyz(3,i,numgrid)=cost
           gridwt(i,numgrid)=wthe*wphi
        end do
      end do

 999  continue
      call callstack_pop
      return
      end

