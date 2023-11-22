      subroutine oep1rhs (
     &     nbas, nocc, naux,
     &     eval, auxph, hfxph, rhs)
c
c This routine calculate the right hand side of the first-order OEP
c equation.
c
c In: nbas, nocc, naux, eval, hfxph, auxph
c Out: rhs
c
c Igor Schweigert, Jan 2004
c $Id: oep1rhs.FPP,v 1.1.1.1 2008/05/29 19:35:39 taube Exp $
c
      implicit none
c
c     Arguments
c
      integer
     &     nbas, nocc, naux
c
      double precision
     &     eval (nbas),
     &     auxph (nbas-nocc, nocc, naux),
     &     hfxph (nbas-nocc, nocc), rhs (naux)
c
c     Local variables
c
      integer
     &     x, i, a
c
      call zero (rhs, naux)
      do x=1, naux
         do i=1, nocc
            do a=1, nbas-nocc
               rhs (x) = rhs (x) +
     &              4.d0 * auxph (a, i, x) * hfxph (a, i)/
     &              (eval (i) - eval (a+nocc))
            enddo
         enddo
      enddo
c 
      return
      end
