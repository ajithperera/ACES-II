      subroutine calcgrad(xai,fia,grad,ov,ncoord,ioffset)
      implicit none
      integer ov,ncoord,ioffset
      double precision xai(ov),fia(ov),grad(ncoord),ddot
 
      grad(ioffset)=4.0d0*ddot(ov,fia,1,xai,1)
CSSS       write(*,*) 'ioffset=',ioffset,grad(ioffset)
      return
      end 

