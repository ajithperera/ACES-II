      subroutine oepmkener (
     &     nbas, ntri, nocc, scfitr,
     &     dmax, nucrep, etot, noconv,
     &     scfevl, scfden, scfham, h1eao, hfxao)
c
c This routine computes the OEP SCF (sum of the eigenvalues) or OEP full
c energy (expectation value of the full Hamiltonian with respect to the
c OEP single determinant plus the correlation energy if necessary), and
c print them. Which energy to print (and format) depedends on whether
c the OEP SCF iterations have converged or not (in other words if DMAX >
c TOL).
c
c Igor Schweigert, Sep 2003, $Id: oepmkener.FPP,v 1.1.1.1 2008/05/29 19:35:39 taube Exp $
c 
      implicit none
c
c     Arguments
c 
      integer
     &     nbas, ntri, nocc, scfitr
c
      logical noconv
c 
      double precision
     &     dmax, nucrep, etot,
     &     scfevl (*), scfden (*), scfham (*),
     &     h1eao (*), hfxao (*)
c
c     Local variables
c 
      integer
     &     i
c
      double precision
     &     ddot, oep0, oep1, anarmyofone
c 
      data anarmyofone /1.d0/
c 
c     Messages to the user
c     
 5000 format(t3,i6,5x,f20.10,10x,d20.10)
 5005 format(t3,1x,a10,f20.10,10x,d20.10,/)
 5010 format(t3,1x,a10,f20.10,/)
c
c     Calculate and print the expectation value of the OEP Hamiltonian
c     with respect to the OEP single determinant (just the sum of the
c     eigenvalues plus nuclear repulsion.)
c
      oep0 = nucrep +
     &     ddot (nocc, scfevl, 1, anarmyofone, 0)
c
c     If the OEP SCF iterations have not converged, the SCF energy is
c     all we need, so print it together with the density difference and
c     return. It seems SCF energy is used in DMPSCF, so pass it back
c     through ETOT.
c     
      if (noconv) then
         write (6, 5000) scfitr, oep0, dmax
         etot = oep0
         return
      endif
c
c
c     Calculate the first-order OEP energy (expectation value of the
c     full Hamiltonian with repect to the OEP single determinant). The
c     formula is E = 1/2 sum (mu,nu) D (mu,nu) x (H1e (mu,nu) + F (mu,
c     nu)). Since all matrices are in the triangular format, we will
c     scale the diagonal AOMEs of density by a factor 1/2, perform the
c     summation without prefactor 1/2, and then scale density back:
c     
c       1/ scale the density diagonal elements by a factor of 1/2
c       2/ contract the density with one-particle Hamiltonian, Hartree
c          Hamiltonian, and non-local exchange, and print
c       3/ scale the density diagonal elements by a factor of 2 
c
      do i=1, nbas
         scfden (i*(i+1)/2) = 0.5d0 * scfden (i*(i+1)/2)
      enddo
c
      oep1 = nucrep
      oep1 = oep1 +
     &     ddot (ntri, scfden, 1, h1eao, 1)  +
     &     ddot (ntri, scfden, 1, scfham, 1) +
     &     ddot (ntri, scfden, 1, hfxao,  1)
c
      do i=1, nbas
         scfden (i*(i+1)/2) = 2.0d0 * scfden (i*(i+1)/2)
      enddo
c
c     Print the SCF and total OEP energy
c     
      write (6, 5005) 'E0 = ', oep0, dmax
      write (6, 5010) 'E (OEP) = ', oep1
c
c     The total energy will be used later on in the post-SCF block of
c     VSCF, pass it back in ETOT.
c     
      etot = oep1
c     
      return
      end
