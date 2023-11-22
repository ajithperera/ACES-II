      subroutine mn_invert(a, ainv, n, ipiv)
c
c the matrix a is inverted. The inverted matrix is obtained in ainv
c
      implicit none
c
      integer n, ipiv(n), info, i
      double precision
     &   a(n,n), ainv(n,n)
c
cmn      write(6,*) ' @ invert'
cmn      call output(a, 1,n,1,n,n,n,1)
      call zero(ainv, n*n)
      do i = 1, n
        ainv(i,i) = 1.0d0
      enddo
cmn      call output(ainv, 1,n,1,n,n,n,1)
c
      info = 0
      call DGESV(n, n, a, n, ipiv, ainv, n, info)
c
      if (info .ne. 0) then
        write(6,*) ' @invert : something wrong after dgesv', info
      call output(a, 1,n,1,n,n,n,1)
      call output(ainv, 1,n,1,n,n,n,1)
c        call errex
      endif
c
      return
      end
c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End: c
