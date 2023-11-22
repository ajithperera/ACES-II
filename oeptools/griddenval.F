      subroutine griddenval (
     &     nprm, nbas,
     &     prmval, prm2ao, dens, scrax1, denval)
c
c This routine calculates the value of the density at the given point in
c the real space. 
c
c Igor Schweigert, Nov 2002
c $Id: griddenval.FPP,v 1.1.1.1 2008/05/29 19:35:40 taube Exp $
c
      implicit none
c
c     Arguments.
c
      integer
     &     nprm, nbas
c
      double precision
     &     prmval (nprm), prm2ao (nprm, nbas), dens (nbas, nbas),
     &     scrax1 (nbas), denval
c
c     Local variables.
c
      integer i,j
c
c     Calculate the value of the density at a given point in space: <>
c     Contract the primitive values with contraction coefficients to get
c     the AO values. <> Contract the AO values with the density matrix
c     to get the density value.
c
      call dgemm(
     &     'n', 'n', 1, nbas, nprm,
     &     1.d0, prmval, 1, 
     &     prm2ao, nprm,
     &     0.d0, scrax1, 1)
c
      denval = 0.d0
      do i=1, nbas
         do j=1, nbas
            denval = denval + dens (i,j) * scrax1 (i) * scrax1 (j)
         enddo
      enddo
c
      return
      end
