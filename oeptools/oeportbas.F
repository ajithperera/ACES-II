      subroutine oeportbas (
     &     nbas, nwrk, nind, thres, eval, work, ovrl)
c
c This routine constructs the matrix that transforms a given basis to
c the corresponding orhtonormal, linearly independent basis. (Such a
c matrix constist of the overlap matrix eigenvectors with non-zero
c eigenvalues.
c
c Only the upper diagonal of the symmetric overlap matrix is referenced.
c
c Note that, if the original basis has linear dependencies then the
c transformation matrix is rectangular with the second index
c corresponding to the new basis.
c
c In: nbas, nwork, thres
c Scr: eval, work
c Out: nind, ovrl is rewritten with the transformation matrix
c
c Igor Schweigert, Oct 2003
c $Id: oeportbas.FPP,v 1.1.1.1 2008/05/29 19:35:40 taube Exp $
c 
      implicit none
c
c     Arguments
c 
      integer
     &     nbas, nwrk, nind
c     
      double precision
     &     thres, eval (nbas), work (nwrk), ovrl (*)
c
c     Local variables
c     
      integer
     &     ierr, n
c
c     Messages to the user
c     
 5905 format (
     &     t3,'@OEPORTBAS-F: DSYEV failed to converge an eigenvalue',/)
 5910 format (
     &     t3,'@OEPORTBAS-F: Illegal argument to DSYEV, argument #',i2,
     &     /)
 5930 format (
     &     t3,'@OEPORTBAS-I: ',i4,' linear dependent functions (with',/,
     &     t3,'eigenvalues below ',d7.1,' ) have been removed from',/,
     &     t3,'the basis.',/)
c
c     Calculate the matrix that transform the original basis to the
c     orthogonal one: <> Diagonalize the overlap matrix. <> Check if
c     the diagonalization procedure converged, if not complain and die.
c     
      call dsyev (
     &     'v', 'u', nbas, 
     &     ovrl, nbas,
     &     eval, work, nwrk, ierr)
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
c     Calculate the number of linearly independent functions and collect
c     all the corresponding eigenvectors (ones with eigenvalues beyond
c     the threshold) in the beginning of the transformation matrix: <>
c     Loop over eigenvectors: <<>> If the current eigenvector has
c     eigenvalue below the threshold swap it with the last non-zero
c     eigenvalue, swap the eigenvalues correspondingly, decrement the
c     number of lin. indep. functions and redo the loop. <<>> If not,
c     increment the loop index and cycle. <> Check if the last vector
c     has eigenvalue below the threshold and decrement the number of
c     lin. indep. functions if it does (because the last vector has
c     never got checked.) <> Report the number of lin. dep. functions
c     discarded.
c
c     Note using DSYEV assumes NWRK > 3 * NBAS, so we can use NWRK to
c     swap eigenvectors. Also it is convenient to treat eigenvectors as
c     columns (to easy copy them from one matrix to another), so once a
c     certain number of eigenvectors is disregarded, the second index of
c     the transformation matrix correponds to the new, linearly
c     independent basis.
c     
      n = 1
      nind = nbas
      do while (n.lt.nind)
         if (dabs (eval (n)) .lt. thres) then
            work (1) = eval (n)
            eval (n) = eval (nind)
            eval (nind) = work (1)
            call dcopy (
     &           nbas,
     &           ovrl (1+(n-1)*nbas), 1,
     &           work, 1)
            call dcopy (
     &           nbas,
     &           ovrl (1+(nind-1)*nbas), 1,
     &           ovrl (1+(n-1)*nbas), 1)
            call dcopy (
     &           nbas,
     &           work, 1,
     &           ovrl (1+(nind-1)*nbas), 1)
            nind = nind - 1
         else
            n = n + 1
         endif
      enddo
c     
      if (dabs (eval (nind)) .lt. thres) nind = nind - 1
c     
      if (nind .lt. nbas)
     &     write (6, 5930) nbas - nind
c
c     Normalize the basis: <> Loop over the eigenvectors and scale them
c     with the square root of the correponding eigenvalue.
c     
      do n=1, nind
         call dscal (
     &        nbas, dsqrt (dabs (eval (n))),
     &        ovrl (1+(n-1)*nbas), 1)
      enddo
c     
      return
      end
