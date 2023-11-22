      subroutine pack_coeffs_oed(alpha,ixalpha,pcoeff,ixpcoeff, 
     *                           ncfps,npfps,m,n,alpha_pack,
     *                           nalpha_pack,pcoeff_pack,
     *                           npcoef_pack,ccbeg,ccend,
     *                           indx_cc,ccbeg_pack,ccend_pack, 
     *                           Nprims,Nshells,Nalpha,Npcoef,
     *                           Max_cfps,Max_pfps)
c---------------------------------------------------------------------------
c   Formats the integral exponents and contraction coefficients for use
c   in the ERD integral package.
c---------------------------------------------------------------------------
      implicit none
      integer ncfps(Nshells),npfps(Nshells),ixalpha(Nshells)
      integer ixpcoeff(Nshells)
      integer ccbeg(Nprims),ccend(Nprims),indx_cc(Nshells)
      integer ccbeg_pack(4*Max_cfps),ccend_pack(4*Max_cfps)
      integer m, n, r, s
      integer Max_cfps,Max_pfps
      integer nalpha,npcoef,Nshells,Nprims
      integer nalpha_pack,npcoef_pack

      integer ialpha
      integer ipcoeff
      integer quad(2)
      double precision alpha(Nalpha),pcoeff(Npcoef)
      double precision alpha_pack(4*Max_pfps)
      double precision pcoeff_pack(4*Max_pfps*Max_cfps)
      integer i, j, num, ishell
      integer k, l, icc

      quad(1) = m
      quad(2) = n

      ialpha  = 1
      ipcoeff = 1
      nalpha_pack = 0
      npcoef_pack = 0

      do j = 1, 2 
         ishell = quad(j)
         do i = 1, npfps(ishell)
            alpha_pack(ialpha+i-1) = alpha(ixalpha(ishell)+i-1)
         enddo
         ialpha = ialpha + npfps(ishell)
         nalpha_pack = nalpha_pack + npfps(ishell)
      enddo
  
      icc = 1
      do j = 1, 2
         ishell = quad(j)
         num    = npfps(ishell)*ncfps(ishell)
         do i = 1, num
            pcoeff_pack(ipcoeff+i-1) = pcoeff(ixpcoeff(ishell)+i-1)
         enddo

         ipcoeff = ipcoeff + num
         npcoef_pack = npcoef_pack + num
      enddo

      icc = 1
      do j = 1, 2
         ishell = quad(j)

            do k = 1, ncfps(ishell)
               ccbeg_pack(icc) = ccbeg(indx_cc(ishell)+k-1)
               ccend_pack(icc) = ccend(indx_cc(ishell)+k-1)
               icc = icc + 1
            enddo
      enddo

      return
      end  
