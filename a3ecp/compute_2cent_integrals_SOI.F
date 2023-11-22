      subroutine compute_2cent_integrals_SOI
     &                              (a1,a2,b1,b2,scr,
     &                              maxblk,iscr,coords,coeffs,
     &                              alphas, ccbeg, 
     &                              ccend, end_nfps, nshells, nbasis,
     &                              zmax, intmax, npfps, ncfps, 
     &                              ivangmom,
     &                              ixalpha, indx_cc, Atm_4shell,
     &                              ixpcoef,ispherical, outl0_Im,
     *                              outl1_R,outl1_Im, nsend, 
     &                              Nalpha, Npcoef,ncenters,charge,
     *                               max_centers)
c---------------------------------------------------------------------------
c   "Work" routine for the integral worker task.
c---------------------------------------------------------------------------
      implicit none
      integer a1, a2, b1, b2,max_centers 
      integer aa1,aa2,bb1,bb2
      integer adim, bdim, cdim, ddim  
      integer m1, m2, n1, n2
      integer i, j, n, m, r, s
      integer a,b,c,d
      integer nshells, nbasis, zmax, intmax, ispherical

      integer num_to_do, nsend
      integer nints, maxblk
      integer nalpha_pack, npcoeff_pack,npcoeff
      integer ncsum, next, nfirst
      integer me, ierr, nalpha, npcoef, nprim_quad, ncoef_quad
      integer matm, natm, ratm, satm

      integer imin, zmin, iblk, zblk, iatm
      integer ivangmom(*), ixalpha(*), ixpcoef(*), npfps(*)
      integer ncfps(*), indx_cc(*), Atm_4shell(*), ERD_index

      logical skip
      logical mn_symmetry
      logical rs_symmetry
      logical mn_rs_symmetry
      logical*8 l8true, l8spherical
      logical spherical, erfints

      double precision x1,y1,z1
      double precision x2,y2,z2

cPV we need this for NAI Lo, L-1, L+1
      integer ncenters
      double precision xn(max_centers), yn(max_centers),
     *                 zn(max_centers),charge(max_centers),
     *                 ERD_Scale
      integer iatom
      integer DERCX,DERCY,DERCZ
      integer DER1X, DER1Y,DER1Z,IXDERC,DER2X,DER2Y,DER2Z
      integer Ea1,Ea2,Eb1,Eb2
      logical prakash
c end PV
      double precision coords(3,*), coeffs(*), alphas(*)
      double precision outl0_Im(a1:a2,b1:b2),
     *                 outl1_R(a1:a2,b1:b2),outl1_Im(a1:a2,b1:b2)
      double precision scr(*)   
      double precision erfcoef
      integer iscr(*)

      integer ccbeg(*), ccend(*), end_nfps(*)

      integer max_dim_coeff
      parameter (max_dim_coeff = 5000)
      integer ccbeg_pack(max_dim_coeff), ccend_pack(max_dim_coeff)
      double precision alpha_pack(max_dim_coeff), 
     *                 pcoeff_pack(max_dim_coeff)
      save me,alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack

      adim = Nbasis
      bdim = Nbasis

      l8true = .true.
      spherical = (ispherical .eq. 1)
      l8spherical = spherical
          
      nsend = adim*bdim
      if (nsend .lt. 0) then
         print *,'ERROR IN INTEGRAL WORKER ',me,' nsend = ',nsend
         print *,'adim,bdim = ',adim,bdim
         call Errex
      endif

c----------------------------------------------------------------------------
c   Clear the output array.
c----------------------------------------------------------------------------

C
      do b = 1, Nbasis
      do a = 1, Nbasis
         outl0_Im(a,b) = 0.d0
      end do
      end do
  
      do b = 1, Nbasis
      do a = 1, Nbasis
         outl1_R(a,b) = 0.d0
      end do
      end do

      do b = 1, Nbasis
      do a = 1, Nbasis
         outl1_Im(a,b) = 0.d0
      end do
      end do

C
      Write(6,*)
      Write(6,*) "The atom label for each shell"
      Write(6,"(4(1x,I4))") (Atm_4shell(i), i=1, Nshells)
      Write(6,*)

c-----------------------------------------------------------------------
c   Find the shell blocks for which we shall loop through.
c-----------------------------------------------------------------------

         call lookup_shell(end_nfps, nshells, a1, m1)
         call lookup_shell(end_nfps, nshells, a2, m2)
         call lookup_shell(end_nfps, nshells, b1, n1)
         call lookup_shell(end_nfps, nshells, b2, n2)
      write(*,*) "I am here"

       do m=1,nshells
         iatom=Atm_4shell(m)
         xn(iatom) = coords(1,iatom)
         yn(iatom) = coords(2,iatom)
         zn(iatom) = coords(3,iatom)
       end do


      IXDERC = 0
      DERCX  = 0
      DERCY  = 0
      DERCZ  = 0
      write(*,*) "I am here"

      do m = m1, m2
         matm = Atm_4shell(m)
         x1 = coords(1,matm)
         y1 = coords(2,matm)
         z1 = coords(3,matm)
         do n = n1, n2
            natm = Atm_4shell(n)
            x2 = coords(1,natm)
            y2 = coords(2,natm)
            z2 = coords(3,natm)
CSSS            write(*,*) matm,natm
            call pack_coeffs_oed(alphas, ixalpha, coeffs, ixpcoef, 
     *                           ncfps, npfps, m, n, 
     *                           alpha_pack, nalpha, 
     *                           pcoeff_pack, npcoeff, 
     *                           ccbeg, ccend, indx_cc,
     *                           ccbeg_pack, ccend_pack,max_dim_coeff)

      ncsum = ncfps(m) + ncfps(n) 

      Write(6,*) "@-entry to OED_batch"
      Write(6,"((a,2(1x,i4)))") "m, n :", m, n

      DER1X = 1
      DER1Y = 0
      DER1Z = 0
      DER2X = 0
      DER2Y = 1
      DER2Z = 0
  
      write(*,*) "I am here"
      call oed__soi_batchXY(intmax,zmax,
     *                 nalpha,npcoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m),ivangmom(n),x1,y1,z1,x2,y2,z2,
     *  ncenters, xn, yn,zn, charge, ixderc, !ncenters,charge,ixderc,
     *                 der1x,der1y,der1z,        
     *                 der2x,der2y,der2z,
     * dercx,dercy,dercz, alpha_pack,!dercx,dercy,derz,ixderc,
     * pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .true., iscr,
     +                 nints, ! nints no. of integrals in scr
     +                 nfirst,! first point 
     +                 scr) ! integral are coming in scr
      if(nints .gt. 0) then
         Ea2 = end_nfps(m)
         if(m .eq. 1) then
            Ea1=1
         else
            Ea1 = end_nfps(m-1)+1
         end if
         Eb2 = end_nfps(n)
         if(n .eq. 1) then
           Eb1 = 1
         else
           Eb1 = end_nfps(n-1)+1
         end if
      call add_integrals2(outl0_Im,a1,a2,b1,b2,
     *                   scr(nfirst),Ea1,Ea2,Eb1,Eb2,1.0d0)
      endif ! nints 
       write(*,*) "I am here"

      DER1X = 0
      DER1Y = 1
      DER1Z = 0
      DER2X = 1
      DER2Y = 0
      DER2Z = 0

      call oed__SOI_batchYX(intmax,zmax,
     *                 nalpha,npcoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m),ivangmom(n),x1,y1,z1,x2,y2,z2,
     *  ncenters, xn, yn,zn, charge, ixderc, !ncenters,charge,ixderc,
     *                 der1x,der1y,der1z,
     *                 der2x,der2y,der2z,
     * dercx,dercy,dercz, alpha_pack,!dercx,dercy,derz,ixderc see above 
     * pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .true., iscr,
     +                 nints, ! nints no. of integrals in scr
     +                 nfirst,! first point 
     +                 scr) ! integral are coming in scr
      if(nints .gt. 0) then
         Ea2 = end_nfps(m)
         if(m .eq. 1) then
            Ea1=1
         else
            Ea1 = end_nfps(m-1)+1
         end if
         Eb2 = end_nfps(n)
         if(n .eq. 1) then
           Eb1 = 1
         else
           Eb1 = end_nfps(n-1)+1
         end if
      call add_integrals2(outl0_Im,a1,a2,b1,b2,
     *                   scr(nfirst),Ea1,Ea2,Eb1,Eb2,-1.0d0)
       endif ! nints
       end do !n
       end do !m
#ifdef _DEBUG_LVL0
      Write(6,*) "The <mu| SOIZ |nu > integral matrix"
      do b = 1, Nbasis
      Write(6,*)
         Write(6, "(5(1x,F15.10))") (outl0_Im(a,b), a=1, nbasis)
      enddo
#endif

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c Lets for Lx 
c Equation 34 simplified form
c Lx =f/i { I^i_y I^j_z - I^i_z I^j_y}*Ix

c First construct the imaginary part. of L+1 & copy to L-1

      do m = m1, m2
         matm = Atm_4shell(m)
         x1 = coords(1,matm)
         y1 = coords(2,matm)
         z1 = coords(3,matm)
         do n = n1, n2
            natm = Atm_4shell(n)

            x2 = coords(1,natm)
            y2 = coords(2,natm)
            z2 = coords(3,natm)
            call pack_coeffs_oed(alphas, ixalpha, coeffs, ixpcoef,
     *                           ncfps, npfps, m, n,
     *                           alpha_pack, nalpha,
     *                           pcoeff_pack, npcoeff,
     *                           ccbeg, ccend, indx_cc,
     *                           ccbeg_pack, ccend_pack,max_dim_coeff)
      ncsum = ncfps(m) + ncfps(n)

c YZ part
      DER1X = 0
      DER1Y = 1
      DER1Z = 0
      DER2X = 0
      DER2Y = 0
      DER2Z = 1

      call oed__soi_batchYZ(intmax,zmax,
     *                 nalpha,npcoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m),ivangmom(n),x1,y1,z1,x2,y2,z2,
     *  ncenters, xn, yn,zn, charge, ixderc, !ncenters,charge,ixderc
     *                 der1x,der1y,der1z,
     *                 der2x,der2y,der2z,
     * dercx,dercy,dercz, alpha_pack,!dercx,dercy,derz,ixderc see above 
     * pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .true., iscr,
     +                 nints, ! nints no. of integrals in scr
     +                 nfirst,! first point 
     +                 scr) ! integral are coming in scr
      if(nints .gt. 0) then
         Ea2 = end_nfps(m)
         if(m .eq. 1) then
            Ea1=1
         else
            Ea1 = end_nfps(m-1)+1
         end if
         Eb2 = end_nfps(n)
         if(n .eq. 1) then
           Eb1 = 1
         else
           Eb1 = end_nfps(n-1)+1
         end if
      call add_integrals2(outl1_Im,a1,a2,b1,b2,
     *                   scr(nfirst),Ea1,Ea2,Eb1,Eb2,1.0d0)
      endif ! nints
c ZY= (mu^z||nu^y)
      DER1X = 0
      DER1Y = 0
      DER1Z = 1
      DER2X = 0
      DER2Y = 1
      DER2Z = 0
      call oed__SOI_batchZY(intmax,zmax,
     *                 nalpha,npcoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m),ivangmom(n),x1,y1,z1,x2,y2,z2,
     *  ncenters, xn, yn,zn, charge, ixderc, !ncenters,charge,ixderc
     *                 der1x,der1y,der1z,
     *                 der2x,der2y,der2z,
     * dercx,dercy,dercz, alpha_pack,!dercx,dercy,derz,ixderc see above 
     * pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .true., iscr,
     +                 nints, ! nints no. of integrals in scr
     +                 nfirst,! first point 
     +                 scr) ! integral are coming in scr
      if(nints .gt. 0) then
         Ea2 = end_nfps(m)
         if(m .eq. 1) then
            Ea1=1
         else
            Ea1 = end_nfps(m-1)+1
         end if
         Eb2 = end_nfps(n)
         if(n .eq. 1) then
           Eb1 = 1
         else
           Eb1 = end_nfps(n-1)+1
         end if

      call add_integrals2(outl1_Im,a1,a2,b1,b2,
     *                   scr(nfirst),Ea1,Ea2,Eb1,Eb2,-1.0d0)
       endif ! nints
        end do ! n
        end do ! m
#ifdef _DEBUG_LVL0
      Write(6,*) "The <mu| SOIX |nu > integral matrix"
      do b = 1, Nbasis
      Write(6,*)
         Write(6, "(5(1x,F15.10))") (outl1_Im(a,b), a=1, nbasis)
      enddo
#endif

c***************************************************************************
c***************************************************************************
c Ly = { (mu^z || nu^ x) - ( mu^x | nu^z) }
c first term
c (mu^z || nu^ x)
      do m = m1, m2
         matm = Atm_4shell(m)
         x1 = coords(1,matm)
         y1 = coords(2,matm)
         z1 = coords(3,matm)
         do n = n1, n2
            natm = Atm_4shell(n)
            x2 = coords(1,natm)
            y2 = coords(2,natm)
            z2 = coords(3,natm)
            call pack_coeffs_oed(alphas, ixalpha, coeffs, ixpcoef,
     *                           ncfps, npfps, m, n,
     *                           alpha_pack, nalpha,
     *                           pcoeff_pack, npcoeff,
     *                           ccbeg, ccend, indx_cc,
     *                           ccbeg_pack, ccend_pack,max_dim_coeff)
      ncsum = ncfps(m) + ncfps(n)
      DER1X = 0
      DER1Y = 0
      DER1Z = 1
      DER2X = 1
      DER2Y = 0
      DER2Z = 0
      call oed__soi_batchZX(intmax,zmax,
     *                 nalpha,npcoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m),ivangmom(n),x1,y1,z1,x2,y2,z2,
     *  ncenters, xn, yn,zn, charge, ixderc, !ncenters,charge,ixderc
     *                 der1x,der1y,der1z,
     *                 der2x,der2y,der2z,
     * dercx,dercy,dercz, alpha_pack,!dercx,dercy,derz,ixderc see above 
     * pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .true., iscr,
     +                 nints, ! nints no. of integrals in scr
     +                 nfirst,! first point 
     +                 scr) ! integral are coming in scr
      if(nints .gt. 0) then
         Ea2 = end_nfps(m)
         if(m .eq. 1) then
            Ea1=1
         else
            Ea1 = end_nfps(m-1)+1
         end if
         Eb2 = end_nfps(n)
         if(n .eq. 1) then
           Eb1 = 1
         else
           Eb1 = end_nfps(n-1)+1
         end if
      call add_integrals2(outl1_R,a1,a2,b1,b2,
     *                   scr(nfirst),Ea1,Ea2,Eb1,Eb2,1.0d0)
      endif ! nints
c - ( mu^x | nu^z)
      DER1X = 1
      DER1Y = 0
      DER1Z = 0
      DER2X = 0
      DER2Y = 0
      DER2Z = 1
      call oed__soi_batchXZ(intmax,zmax,
     *                 nalpha,npcoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m),ivangmom(n),x1,y1,z1,x2,y2,z2,
     *  ncenters, xn, yn,zn, charge, ixderc, !ncenters,charge,ixderc
     *                 der1x,der1y,der1z,
     *                 der2x,der2y,der2z,
     * dercx,dercy,dercz, alpha_pack,!dercx,dercy,derz,ixderc see above 
     * pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .true., iscr,
     +                 nints, ! nints no. of integrals in scr
     +                 nfirst,! first point 
     +                 scr) ! integral are coming in scr
      if(nints .gt. 0) then
         Ea2 = end_nfps(m)
         if(m .eq. 1) then
            Ea1=1
         else
            Ea1 = end_nfps(m-1)+1
         end if
         Eb2 = end_nfps(n)
         if(n .eq. 1) then
           Eb1 = 1
         else
           Eb1 = end_nfps(n-1)+1
         end if
c Add {(mu^z||nu^x) - (mu^x||nu^z)} Imaginary part
      call add_integrals2(outl1_R,a1,a2,b1,b2,
     *                   scr(nfirst),Ea1,Ea2,Eb1,Eb2,-1.0d0)
      endif ! nints
      end do ! n
      end do ! m
#ifdef _DEBUG_LVL0
      Write(6,*) "The <mu| SOIY |nu > integral matrix"
      do b = 1, Nbasis
      Write(6,*) 
         Write(6, "(5(1x,F15.10))") (outl1_R(a,b), a=1, nbasis)
      enddo
#endif
      return
      end

