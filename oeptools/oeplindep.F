      subroutine oeplindep (
     &     n, nwork, thres, eval, work, a)
c
c This routine finds the eigenvalues of the symmetric matrix A (only the
c upper half is referenced) using LAPACK routine DSYEV, and reports the
c smallest eigenvalue and number of eigenvalues below a certain
c threshold.
c
c In: n, nwork, thres, a
c Scr: work
c Out: eval, a is overwritten
c
c Igor Schweigert, Oct 2003
c $Id: oeplindep.FPP,v 1.1.1.1 2008/05/29 19:35:40 taube Exp $
c
      implicit none
c
c     Arguments
c
      integer
     &     n, nwork
c
      double precision
     &     thres, eval (n), work (nwork), a (n, n)
c
c     Local variables
c
      integer
     &     ierr, i, nmin
c
      double precision
     &     dmin, dabs
c
c     Messages to the user
c
 5905 format (
     &     t3,'@OEPLINDEP-F: DSYEV failed to converge an eigenvalue',/)
 5910 format (
     &     t3,'@OEPLINDEP-F: Illegal argument to DSYEV, argument #',i2,
     &     /)
 5930 format (
     &     t3,'@OEPLINDEP-I: Among ',i5,' eigenvalues, ',i4,' are',/,
     &     t3,'are below ',d7.1,' and the smallest one is ',d20.10/)
c
c     Find the linear dependencies in matrix A: Diagonilize the matrix.
c     <> Check if the diagonalization procedure converged, if not
c     complain and die. <> Find and print out the smalled eigenvalue,
c     and the number of eigenvalues below the threshold.
c
      call dsyev (
     &     'n', 'u', n, 
     &     a, n,
     &     eval, work, nwork, ierr)
c
      if (ierr.ne.0) then
         if (ierr.gt.0) then
            write (6, 5905)
         else
            write (6, 5910) abs (ierr)
         endif
         call errex
      endif
c
      dmin = 1.d+10
      nmin = 0
      do i=1, n
         if (dabs (eval (i)) .lt. dmin) dmin = dabs (eval (i))
         if (dabs (eval (i)) .lt. thres) nmin = nmin + 1
      enddo
c
      write (6, 5930) n, nmin, thres, dmin
c
      return
      end
