      subroutine get_matrix(A, nrow, ncol, perline, iunit)
c     
c     This subroutine gets a matrix from a file: sequential, formatted.
c     The matrix is written in  'perline' columns per line
c     an empty line after each set of columns
c     The format is given in label 999
c     
      implicit none
      integer nrow, ncol, perline, iunit
      double precision A(nrow, ncol)
      integer nleft, jbegin, jend, i, j, nput
c
      if (perline .gt. 8) then
         write(6,*) ' @put_matrix:  maximum perline = 8 '
         call errex
      endif
c     
      nleft = ncol
      jbegin= 1

 100  nput = min(perline, nleft)
      jend = jbegin + nput -1
      do i = 1, nrow
         read(iunit, 999) (A(i, j), j=jbegin,jend)
      enddo
      read(iunit, *)
c
      jbegin = jbegin + nput
      nleft = nleft - nput
      if (nleft .gt. 0) goto 100
c
 999  format(8E20.14)
c
      return
      end
