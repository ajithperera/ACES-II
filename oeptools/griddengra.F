      subroutine griddengra (
     &     nprm, naos, aosval, prmgra, prm2ao, 
     &     dens, aosgra, rhogra)
c
c This routine calculates the value of the density gradient at 
c the given grid point. It calculates the values of AO gradients
c from the primitive gradients and then contracts the former ones
c with the density matrix. 
c
c Igor Schweigert, Dec 2002
c $Id: griddengra.FPP,v 1.1.1.1 2008/05/29 19:35:40 taube Exp $
c
      implicit none
c
c     .. input parameters ..
c
      integer
     &     nprm, naos
c
      double precision
     &     aosval (naos), prmgra (nprm, 3), prm2ao (naos, nprm),
     &     dens (naos, naos), aosgra (naos, 3), rhogra (3)
c
c     .. local variables ..
c
      integer i,j,k
c
c     .. loop over the 3 components of the gradient ..
c
      do k=1, 3
c
c        .. calculate the values of AO gradients at the grid point ..
c
         call xgemm(
     &        'n', 'n', 1, naos, nprm,
     &        1.d0, prmgra (1,k), 1, 
     &        prm2ao, nprm,
     &        0.d0, aosgra (1,k), 1)
c
c        .. calculate the value of the density gradient ..
c
         rhogra (k) = 0.d0
         do i=1, naos
            do j=1, naos
               rhogra (k) = rhogra (k) + 
     &              dens (i,j) * (
     &              aosgra (i ,k) * aosval (j) +
     &              aosgra (j ,k) * aosval (i) )
            enddo
         enddo

      enddo
c
c     .. done ..
c
      return
      end
