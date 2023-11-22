      subroutine scaled(d,ntot,n)
      implicit double precision(a-h,o-z)
      integer dirprd
      dimension d(ntot),n(8)
      common/syminf/nstart,nirrep,irrepa(255,2),dirprd(8,8)
C
      data two,half /2.d0,0.5d0/
C
      call sscal(ntot,two,d,1)
      ioff=0
      do 100 irrep=1,nirrep
       do 110 i=1,n(irrep)
        d(ioff+i*(i+1)/2) = d(ioff+i*(i+1)/2)*half
110    continue
       ioff=ioff+n(irrep)*(n(irrep)+1)/2
100   continue
      return 
      end
