      subroutine mn_fndmaxd(n, vec, diff, iloc, tol)
c
c  locate 'first' maximum absolute element in vec, and its positon.
c This means that values within tol are considered the same.
c
      implicit none
c
      integer n, iloc, i
      double precision
     &   vec(n), diff, amax, tol
c
      iloc=1
      if (n .lt. 1) then
        diff = 0.0d0
        goto 999
      endif
c
      amax = -100.0
      do i = 1, n
        if (abs(vec(i)) .gt. (amax + tol)) then
          amax = abs(vec(i))
          iloc = i
        endif
      enddo
c
      diff = abs(vec(iloc))
c
 999  return
      end
c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End: c
