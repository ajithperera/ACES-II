      subroutine calcgradBB(xaia,xaib,fia,fib,grad,ova,
     &                     ovb,ncoord,ioffset)
      implicit none
      integer ova,ovb,ncoord,ioffset
      double precision xaia(ova),xaib(ovb),fia(ova),
     &                 fib(ovb),grad(ncoord),ddot
 
      grad(ioffset)=2.0d0*ddot(ova,fia,1,xaia,1)+
     &               2.0d0*ddot(ovb,fib,1,xaib,1)

CSSS       write(*,*) 'ioffset=',ioffset,grad(ioffset)
      return
      end 

