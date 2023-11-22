      subroutine getblk2(wfull,wblock,nrowf,ncolf,nrowb,ncolb,
     &                   irowoff,icoloff)
      implicit none
      double precision wfull,wblock
      integer nrowf,ncolf,nrowb,ncolb,irowoff,icoloff
      integer i,j
C
      dimension wfull(nrowf,ncolf),wblock(nrowb,ncolb)
C
      do 20 j=1,ncolb
      do 10 i=1,nrowb
      wblock(i,j) = wfull(irowoff-1+i,icoloff-1+j)
   10 continue
   20 continue
      return
      end
