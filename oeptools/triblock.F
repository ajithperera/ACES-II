      subroutine triblock (
     &     na, nb, a, b, c, d)
c
c This routine given two symmetric matrices A, B, and a rectungular
c matrix C constructs the upper diagonal of a symmetric matrix D with A
c and B as the diagonal blocks and C as off-diagonal. Note that only
c upper diagonals of matrices A,B, and D are referenced.
c
c D = A  C
c        B 
c
c The matrix dimensions should be as follows:
c A (na, na), B (nb, nb), C (na, nb), D (na+nb, na+nb)
c
c Igor Schweigert, Oct 2003
c $Id: triblock.FPP,v 1.1.1.1 2008/05/29 19:35:40 taube Exp $
c
      implicit none
c
c     Arguments
c
      integer
     &     na, nb
c
      double precision
     &     a (*), b (*), c (*), d (*)
c
c     Local variables
c
      integer
     &     n, nd
c
c     Construct D: <> Loop over first NA row of D, copy rows of A and C
c     into D <> Loop over remaining NB rows of D, copy B into D.
c
      nd = na + nb
      do n=1, na
         call dcopy (
     &        na-n+1,
     &        a (n+(n-1)*na), na,
     &        d (n+(n-1)*nd), nd)
         call dcopy (
     &        nb,
     &        c (n), na,
     &        d (na*nd+n), nd)
      enddo
c
      do n=1, nb
         call dcopy (
     &        nb,
     &        b (n+(n-1)*nb), nb,
     &        d (na*nd+na+n+(n-1)*nd), nd)
      enddo
c
      return
      end
