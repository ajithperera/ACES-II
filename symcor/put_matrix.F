      subroutine put_matrix(A, nrow, ncol, perline, iunit, eformat)
c     
c     This subroutine puts a matrix on a file: sequential, formatted.
c     The matrix is written in  'perline' columns per line
c     an empty line after each set of columns
c     The format is given in label 999
c     
      implicit none
      integer nrow, ncol, perline, iunit
      logical eformat
      double precision A(nrow, ncol)
      integer nleft, jbegin, jend, i, j, nput
c
      if (perline .gt. 8) then
         write(6,*) ' @put_matrix:  maximum perline = 8 '
         call aces_exit(1)
      endif
c     
      nleft = ncol
      jbegin= 1

 100  nput = min(perline, nleft)
      jend = jbegin + nput -1
      do i = 1, nrow
         if (eformat) then
         write(iunit, 998) (A(i, j), j=jbegin,jend)
         else
         write(iunit, 999) (A(i, j), j=jbegin,jend)
         endif
      enddo
      write(iunit, *)
c
      jbegin = jbegin + nput
      nleft = nleft - nput
      if (nleft .gt. 0) goto 100
c
 998  format(8E14.6)
 999  format(8F12.6)
c
      return
      end
