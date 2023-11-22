      Subroutine Erd_scratch_mem_calc(Nshells,Nprims,Nalpha,Npcoef,
     *                                shtype,ncfps,npfps, 
     *                                atom,coords,alpha,alpha_pack,
     *                                pcoeff,pcoeff_pack,
     *                                ixalpha,ixpcoeff,ccbeg,
     *                                ccbeg_pack,ccend,ccend_pack,
     *                                indx_cc,spherical,
     *                                nuclei,Max_cfps,Max_pfps,
     *                                imax,zmax)
c-------------------------------------------------------------------------
c   Calculates the maximum amount (in words) of scratch memory required
c   by the Flocke integral package.  The calculation is carried out in 
c   parallel over the number of processors in the job. 
c
c   Arguments:
c	nshells			Number of shells 
c       shtype			Array of shell types (angular momentum)
c	ncfps			Number of contracted functions per shell.
c	npfps			Number of primitive functions per shell.
c       atom                    Array of atomic labels (one for each shell).
c       coords                  Array of coords per shell.
c	imax			Returned number of integer words required.
c	zmax 			Returned number of fl. pt. words required.
c       me                      Rank of current processor.
c--------------------------------------------------------------------------

      implicit none

      integer nshells 
      integer shtype(nshells),ncfps(nshells),npfps(nshells)
      integer atom(nshells)
      integer ixalpha(nshells), ixpcoeff(nshells)
      integer ccbeg(Nprims),ccend(Nprims)
      integer ccbeg_pack(4*Max_cfps),ccend_pack(4*Max_cfps)
      integer indx_cc(Nshells)
      integer imax,zmax,Max_cfps,Max_pfps
      integer nuclei,nalpha,ncoeff,Nprims,Npcoef
      integer nalpha_pack,Npcoef_pack

      logical calc_2der

      double precision coords(3,nshells)
      double precision alpha(Nalpha),alpha_pack(4*Max_pfps)
      double precision pcoeff(Npcoef)
      double precision pcoeff_pack(4*Max_pfps*Max_cfps)

      integer i,j
      integer m, n, r, s
      integer iopt, zopt 
      integer imin, zmin
      integer ixderc
      integer der1x, der1y, der1z
      integer der2x, der2y, der2z
      integer dercx, dercy, dercz
      integer der_flags(12), iflag, jflag 
      integer lmax, ncmax, npmax
      integer m1, m2, mrange, mleft

      logical spherical

      imax = 0
      zmax = 0

      do m = 1, nshells
      do n = 1, nshells
      do r = 1, nshells
      do s = 1, nshells
         call pack_coeffs(alpha,ixalpha,pcoeff,ixpcoeff,
     *                    ncfps,npfps,m,n,r,s,alpha_pack,
     *                    nalpha_pack,pcoeff_pack,npcoef_pack,
     *                    ccbeg,ccend,indx_cc,ccbeg_pack,
     *                    ccend_pack,Nprims,Nshells,Nalpha,
     *                    Npcoef,Max_cfps,Max_pfps)

c-------------------------------------------------------------------------
c   ERD version 2
c-------------------------------------------------------------------------

         call ERD__MEMORY_ERI_BATCH(nalpha,ncoeff,ncfps(m),ncfps(n),
     *                              ncfps(r),ncfps(s),npfps(m),
     *                              npfps(n),npfps(r),npfps(s),
     *                              shtype(m),shtype(n),shtype(r),
     *                              shtype(s),coords(1,m),coords(2,m),
     *                              coords(3,m),coords(1,n),
     *                              coords(2,n),coords(3,n), 
     *                              coords(1,r),coords(2,r),
     *                              coords(3,r),coords(1,s),
     *                              coords(2,s),coords(3,s), 
     *                              alpha_pack,pcoeff_pack,spherical,
     *                              imin,iopt,zmin,zopt)
         imax = max0(imax, imin)
         zmax = max0(zmax, zmin)
      enddo
      enddo
      enddo
      enddo   

c------------------------------------------------------------------------
c   Now calculate the memory for the one-electron integrals.
c------------------------------------------------------------------------

      do m = 1, nshells
      do n = 1, nshells
        call pack_coeffs_oed(alpha,ixalpha,pcoeff,ixpcoeff,
     *                       ncfps,npfps,m,n,alpha_pack,nalpha,
     *                       pcoeff_pack,ncoeff,ccbeg,ccend,indx_cc,
     *                       ccbeg_pack,ccend_pack,Nprims,Nshells,
     *                       Nalpha,Npcoef,Max_cfps,Max_pfps)

        call oed__memory_nai_batch(nalpha,ncoeff,
     *                             ncfps(m),ncfps(n), 
     *                             npfps(m),npfps(n), 
     *                             shtype(m),shtype(n),
     *                             coords(1,m),coords(2,m),
     *                             coords(3,m),coords(1,n),
     *                             coords(2,n),coords(3,n),
     *                             nuclei,alpha_pack,pcoeff_pack,
     *                             spherical,imin,iopt,zmin,zopt)

        imax = max0(iopt, imax)
        zmax = max0(zopt, zmax)

        call oed__memory_kin_batch(nalpha,ncoeff,
     *                             ncfps(m),ncfps(n), 
     *                             npfps(m),npfps(n), 
     *                             shtype(m),shtype(n),
     *                             coords(1,m),coords(2,m),
     *                             coords(3,m),coords(1,n),
     *                             coords(2,n),coords(3,n),
     *                             alpha_pack,pcoeff_pack,
     *                             spherical,imin,iopt,zmin,zopt)
         imax = max0(iopt, imax)
         zmax = max0(zopt, zmax)

         call oed__memory_ovl_batch(nalpha,ncoeff,
     *                              ncfps(m),ncfps(n), 
     *                              npfps(m),npfps(n), 
     *                              shtype(m),shtype(n),
     *                              coords(1,m),coords(2,m),
     *                              coords(3,m),coords(1,n),
     *                              coords(2,n),coords(3,n),
     *                              alpha_pack,pcoeff_pack,
     *                              spherical,imin,iopt,zmin,zopt)

         imax = max0(iopt, imax)
         zmax = max0(zopt, zmax)

      enddo
      enddo

      return
      end
