      subroutine compute_angmom_integrals
     &                              (a1,a2,b1,b2,scr,
     &                              maxblk,iscr,coords,coeffs,
     &                              alphas, ccbeg, 
     &                              ccend, end_nfps, nshells, nbasis,
     &                              zmax, intmax, npfps, ncfps, 
     &                              ivangmom,
     &                              ixalpha, indx_cc, Atm_4shell,
     &                              ixpcoef,ispherical, out_lx,
     &                              out_ly,out_lz, nsend, 
     &                              Nalpha, Npcoef,ncenters,
     &                              max_centers)
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
      integer ncfps(*), indx_cc(*), Atm_4shell(*)

      logical skip
      logical mn_symmetry
      logical rs_symmetry
      logical mn_rs_symmetry
      logical*8 l8true, l8spherical
      logical spherical, erfints

      double precision x1,y1,z1
      double precision x2,y2,z2
      double precision xc,yc,zc

cPV we need this for NAI Lo, L-1, L+1
      integer ncenters
      integer iatom
      integer DERCX,DERCY,DERCZ
      integer DER1X, DER1Y,DER1Z,IXDERC,DER2X,DER2Y,DER2Z
      integer XMOMDER,YMOMDER,ZMOMDER
      integer Ea1,Ea2,Eb1,Eb2
c end PV
      double precision coords(3,*), coeffs(*), alphas(*)
      double precision out_lx(a1:a2,b1:b2),
     *                 out_ly(a1:a2,b1:b2),out_lz(a1:a2,b1:b2)
      double precision scr(*)   
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
         out_lx(a,b) = 0.d0
         out_ly(a,b) = 0.d0
         out_lz(a,b) = 0.d0
      end do
      end do
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "The atom label for each shell in angmom_int"
      Write(6,"(4(1x,I4))") (Atm_4shell(i), i=1, Nshells)
      Write(6,*)
#endif

c-----------------------------------------------------------------------
c   Find the shell blocks for which we shall loop through.
c-----------------------------------------------------------------------

      call lookup_shell(end_nfps, nshells, a1, m1)
      call lookup_shell(end_nfps, nshells, a2, m2)
      call lookup_shell(end_nfps, nshells, b1, n1)
      call lookup_shell(end_nfps, nshells, b2, n2)
C  
C Set the center ar which the angular momentum integrals are 
C evaluated. Most cases it is at the origin (0.0,0.0,0.0), 
C but it is not necessary.

      XC = 0.0D0
      YC = 0.0D0
      ZC = 0.0D0

C
C The X component of the angular momentum; -{<mu|(y-C)|derz(nu)> - 
C <mu|(z-C)|dery(nu)>}
C
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

CSSS      Write(6,*)
CSSS      Write(6,*) "@-entry to OED_batch"
CSSS      Write(6,"((a,2(1x,i4)))") "m, n :", m, n
      

            DER1X = 0
            DER1Y = 0
            DER1Z = 0
            DER2X = 0
            DER2Y = 0
            DER2Z = 1

            XMOMDER = 0
            YMOMDER = 1
            ZMOMDER = 0
  
            call oed__gener_xyz_derv_batch(intmax,zmax,
     *                 nalpha,npcoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 der1x, der1y, der1z,
     *                 der2x, der2y, der2z,
     *                 xmomder, ymomder, zmomder,
     *                 xc,yc,zc,
     *                 alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .false., iscr, nints, nfirst,
     *                 scr)
CSSS            Write(6,"(a,1x,i4)") "Nints :", Nints

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
C
              call add_integrals2(out_lx,a1,a2,b1,b2,
     *                            scr(nfirst),Ea1,Ea2,Eb1,Eb2,-1.0d0)
            endif ! nints 

           DER1X = 0
           DER1Y = 0
           DER1Z = 0
           DER2X = 0
           DER2Y = 1
           DER2Z = 0

           XMOMDER = 0
           YMOMDER = 0
           ZMOMDER = 1

            call oed__gener_xyz_derv_batch(intmax,zmax,
     *                 nalpha,npcoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 der1x, der1y, der1z,
     *                 der2x, der2y, der2z,
     *                 xmomder, ymomder, zmomder,
     *                 xc,yc,zc,
     *                 alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .false., iscr, nints, nfirst,
     *                 scr)

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
C
              call add_integrals2(out_lx,a1,a2,b1,b2,
     *                            scr(nfirst),Ea1,Ea2,Eb1,Eb2,1.0d0)
            endif ! nints 

       end do !n
       end do !m
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The <mu| L_x |nu > integral matrix"
      call output(out_lx, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
CSSS      do b = 1, Nbasis
CSSS      Write(6,*)
CSSS         Write(6, "(5(1x,F15.10))") (out_lx(a,b), a=1, b)
CSSS      enddo
#endif
C
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C The Y component of the angular momentum; -{<mu|(y-C)|derz(nu)> -
C <mu|(z-C)|derx(nu)>}
C
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
            DER1Z = 0
            DER2X = 0
            DER2Y = 0
            DER2Z = 1

            XMOMDER = 1
            YMOMDER = 0
            ZMOMDER = 0

            call oed__gener_xyz_derv_batch(intmax,zmax,
     *                 nalpha,npcoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 der1x, der1y, der1z,
     *                 der2x, der2y, der2z,
     *                 xmomder, ymomder, zmomder,
     *                 xc,yc,zc,
     *                 alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .false., iscr, nints, nfirst,
     *                 scr)

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
C
              call add_integrals2(out_ly,a1,a2,b1,b2,
     *                            scr(nfirst),Ea1,Ea2,Eb1,Eb2,-1.0d0)
            endif ! nints 

           DER1X = 0
           DER1Y = 0
           DER1Z = 0
           DER2X = 1
           DER2Y = 0
           DER2Z = 0

           XMOMDER = 0
           YMOMDER = 0
           ZMOMDER = 1

            call oed__gener_xyz_derv_batch(intmax,zmax,
     *                 nalpha,npcoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 der1x, der1y, der1z,
     *                 der2x, der2y, der2z,
     *                 xmomder, ymomder, zmomder,
     *                 xc,yc,zc,
     *                 alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .false., iscr, nints, nfirst,
     *                 scr)

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
C
              call add_integrals2(out_ly,a1,a2,b1,b2,
     *                            scr(nfirst),Ea1,Ea2,Eb1,Eb2,1.0d0)
            endif ! nints 

       end do !n
       end do !m
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The <mu| L_y |nu > integral matrix"
      call output(out_ly, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
C      do b = 1, Nbasis
C      Write(6,*)
C         Write(6, "(5(1x,F15.10))") (out_ly(a,b), a=1, nbasis)
C      enddo
#endif
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C The Z component of the angular momentum; -{<mu|(y-C)|derz(nu)> -
C <mu|(z-C)|derx(nu)>}
C
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
            DER1Z = 0
            DER2X = 0
            DER2Y = 1
            DER2Z = 0

            XMOMDER = 1
            YMOMDER = 0
            ZMOMDER = 0

            call oed__gener_xyz_derv_batch(intmax,zmax,
     *                 nalpha,npcoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 der1x, der1y, der1z,
     *                 der2x, der2y, der2z,
     *                 xmomder, ymomder, zmomder,
     *                 xc,yc,zc,
     *                 alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .false., iscr, nints, nfirst,
     *                 scr)

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
C
              call add_integrals2(out_lz,a1,a2,b1,b2,
     *                            scr(nfirst),Ea1,Ea2,Eb1,Eb2,-1.0d0)
            endif ! nints 

           DER1X = 0
           DER1Y = 0
           DER1Z = 0
           DER2X = 1
           DER2Y = 0
           DER2Z = 0

           XMOMDER = 0
           YMOMDER = 1
           ZMOMDER = 0

            call oed__gener_xyz_derv_batch(intmax,zmax,
     *                 nalpha,npcoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 der1x, der1y, der1z,
     *                 der2x, der2y, der2z,
     *                 xmomder, ymomder, zmomder,
     *                 xc,yc,zc,
     *                 alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .false., iscr, nints, nfirst,
     *                 scr)

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
C
              call add_integrals2(out_lz,a1,a2,b1,b2,
     *                            scr(nfirst),Ea1,Ea2,Eb1,Eb2,1.0d0)
            endif ! nints 

       end do !n
       end do !m
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The <mu| L_z |nu > integral matrix"
      call output(out_lz, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
C      do b = 1, Nbasis
C      Write(6,*) 
C         Write(6, "(5(1x,F15.10))") (out_lz(a,b), a=1, nbasis)
C      enddo
#endif
      return
      end

