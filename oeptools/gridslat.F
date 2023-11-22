      subroutine gridslat (
     &     nprm, naos, prm2ao, valao, dens, denval, nucatt,
     &     scr1xa, scraxp, scraxa, sltval)
c
c This routine calculates the AO matrix elements of 
c the orbital-dependent Slater potential at a given grid point.
c
c Igor Schweigert, 2003 
c $Id: gridslat.FPP,v 1.1.1.1 2008/05/29 19:35:40 taube Exp $
c
      implicit none
c
c     .. input parameters ..
c
      integer
     &     naos, nprm
c
      double precision
     &     denval, dens (naos, naos), valao (naos),
     &     prm2ao (nprm, naos), nucatt (nprm, nprm),
     &     scr1xa (naos), scraxp (naos, nprm), scraxa (naos, naos),
     &     sltval
c
c     .. local variables ..
c
      integer i,j
c
c     .. transform the nuc.att. integrals to the AO representation ..
c
      call xgemm(
     &     't', 'n', naos, nprm, nprm,
     &     1.d0, prm2ao, nprm,
     &     nucatt, nprm,
     &     0.d0, scraxp, naos)
c
      call xgemm(
     &     'n', 'n', naos, naos, nprm,
     &     1.d0, scraxp, naos,
     &     prm2ao, nprm,
     &     0.d0, scraxa, naos)
c
c     .. contract the density matrix with the values of the AO's      ..
c     .. note the factor 1/2 to "undo" the spin summation for RHF     ..
c
      call xgemm(
     &     'n', 'n', 1, naos, naos,
     &     .5d0, valao, 1, 
     &     dens, naos,
     &     0.d0, scr1xa, 1)
c
c     .. contract the "density-contracted" AO values with nuc.att. integ.'s..
c     .. to get the value of the Slater potential at the grid point   ..
c     .. factor 2 comes from the RHF spin summation for               ..
c
      sltval = 0.d0
      do i=1, naos
         do j=1, naos
            sltval = sltval -
     &           scr1xa (i) * scraxa (i,j) * scr1xa (j)
         enddo
      enddo
      sltval = 2.d0 * sltval / denval
c
c     .. done ..
c
      return
      end
