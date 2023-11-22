      subroutine oep1enr (
     &     nbas, dens, oneh, hart, hfx, oepxe)
c
c This routine calculates the first-order OEP energy, which is just the
c Hartree-Fock energy expressions with the OEP orbitals.
c
c Igor Schweigert, Feb 2004
c $Id: oep1enr.FPP,v 1.1.1.1 2008/05/29 19:35:39 taube Exp $
c 
      implicit none
c
c     Arguments
c     
      integer nbas
c     
      double precision
     &     dens (*), oneh (*), hart (*), hfx (*), oepxe
c
c     Local variables
c     
      integer n
      double precision ddot
c
c     Calculate the zeroth- and first-order contribution to the total
c     energy: <> Scale the diagonal part of the density matrix by a
c     factor of one half. <> Contract the density matrix with the
c     core-Hamiltonian, core-Hamiltonian plus Coulomb, and exchange
c     AOMEs. <> Scale the diagonal part of the density matrix by a
c     factor of two.
c     
      do n = 1, nbas
         dens (n*(n+1)/2) = .5d0 * dens (n*(n+1)/2)
      enddo
c     
      oepxe = oepxe +
     &     ddot (nbas*(nbas+1)/2, dens, 1, oneh, 1) +
     &     ddot (nbas*(nbas+1)/2, dens, 1, hart, 1) +
     &     ddot (nbas*(nbas+1)/2, dens, 1, hfx,  1)
c     
      do n = 1, nbas
         dens (n*(n+1)/2) = 2.d0 * dens (n*(n+1)/2)
      enddo
c     
      return
      end

