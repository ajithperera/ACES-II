      subroutine oepresinv (
     $     nbas, nocc, naux, svdthr, scfevl, auxph, 
     $     resevl, resfun, resinv)
c
c This routine constructs the response function in the auxiliary basis 
c set representation, and then finds its "singular" inverse using the
c Singular Value Decomposition.
c
c input: nbas, nocc, naux, scfevl, auxph, svdthr
c scratch: resfun, resevl
c output: resinv
c
c Igor Schweigert, Apr 2003
c $Id: oepresinv.FPP,v 1.1.1.1 2008/05/29 19:35:39 taube Exp $
c
      implicit none
c
c     .. the input parameters                                         ..
c
      integer
     $     nbas, nocc, naux, svdthr
c
      double precision
     $     scfevl (nbas), auxph (nbas-nocc, nocc, naux),
     $     resevl (naux), resfun (naux, naux), 
     $     resinv (naux, naux)
c
c     .. local variables                                              ..
c
      integer
     $     i, a, x, y, z, ndep, ierr
c     
      double precision
     &     threshold
c
c     .. messages to the user                                         ..
c
 1005 format (
     $     t3, '@OEPRESINV-F: DSYEV failed to converge an eigenvalue')
 1010 format (
     $     t3, '@OEPRESINV-F: Illegal argument to DSYEV, argument #',
     $     i2)
 1020 format (
     $     t3, '@OEPRESINV-I: ', i2, ' singularities of the response ',
     $     'function have been removed.')
c
c     .. build the response function (the upper triangle)             ..
c
      call zero (resfun, naux*naux)
      do y=1, naux
         do x=1, y
            do i=1, nocc
               do a=nocc+1, nbas
                  resfun (x, y) = resfun (x, y) +
     $                 auxph (a-nocc, i, x) * auxph (a-nocc, i, y) /
     $                 (scfevl (i) - scfevl (a))
               enddo
            enddo
         enddo
      enddo
c
      do y=1, naux
         do x=1, y
            resfun (y, x) = resfun (x, y)
         enddo
      enddo
c
      call dscal (naux**2, 4.d0, resfun, 1)
c
c     .. and diagonalize it. 
c     .. RESINV will be used as the WORK array                        ..
c
      call dsyev (
     $     'v', 'u', naux,
     $     resfun, naux,
     $     resevl, resinv, naux*naux, ierr)
c
c     .. check if eigenvalues have been found                         ..
c
      if (ierr.ne.0) then
         if (ierr.gt.0) then
            write (6, 1005)
         else
            write (6, 1010) abs (ierr)
         endif
         call errex
      endif
c
c     .. invert all the eigenvalues beyond the threshold              ..
c
      threshold = 10.d0 ** (-svdthr)
      ndep = 0
      do x=1, naux
         if ( dabs (resevl (x)) .gt. threshold ) then
            resevl (x) = 1.d0 / resevl (x)
         else
            resevl (x) = 0.d0
            ndep = ndep + 1
         endif
      enddo
c
c     .. let the user know how many singularities will be removed     ..
c     
      if (ndep.gt.0) write (6, 1020) ndep
c
c     .. construct the SVD inverse of the reponse function            ..
c
      call zero (resinv, naux*naux)
      do x=1, naux
         do y=1, naux
            do z=1, naux
               resinv (x, y) = resinv (x, y) +
     $              resfun (x, z) * resevl (z) * resfun (y, z)
            enddo
         enddo
      enddo
c
c     .. Note that if you want to debug this SVD procedure, you will  ..
c     .. need one extra matrix (naux, naux), and then if SVD is done  ..
c     .. right then A = A^(SVD) x A x A = A x A^(SVD) x A = etc. and  ..
c     .. 1 = A^(SVD) A = A x A^(SVD) if and only if A is not singular ..
c
c     .. let's check if we did the SVD properly                       ..
c
c$$$
c$$$      call xgemm (
c$$$     $     'n', 'n', naux, naux, naux,
c$$$     $     1.d0, resinv, naux,
c$$$     $     resfun, naux,
c$$$     $     0.d0, try, naux)
c$$$c
c$$$      call dcopy (
c$$$     $     naux * naux,
c$$$     $     resfun, 1,
c$$$     $     resinv, 1)
c$$$c
c$$$      call xgemm (
c$$$     $     'n', 'n', naux, naux, naux,
c$$$     $     -1.d0, resfun, naux,
c$$$     $     try, naux,
c$$$     $     1.d0, resinv, naux)
c$$$c
c$$$      write (6,*) 'Matrix that follows must be zero'
c$$$      call output (resinv, 1, naux, 1, naux, naux, naux, 1)
c
c
c     .. done ..
c
      return
      end
