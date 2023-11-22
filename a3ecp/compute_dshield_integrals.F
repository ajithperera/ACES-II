      subroutine compute_dshield_integrals
     &                              (a1,a2,b1,b2,scr,
     &                              maxblk,iscr,coords,coeffs,
     &                              alphas, ccbeg, 
     &                              ccend, end_nfps, nshells, nbasis,
     &                              zmax, intmax, npfps, ncfps, 
     &                              ivangmom,
     &                              ixalpha, indx_cc, Atm_4shell,
     &                              ixpcoef,ispherical, 
     &                              out_dsh_xx,out_dsh_yx,
     &                              out_dsh_zx,
     &                              out_dsh_xy,out_dsh_yy,
     &                              out_dsh_zy,
     &                              out_dsh_xz,out_dsh_yz,
     &                              out_dsh_zz,
     &                              nsend, 
     &                              Nalpha, Npcoef,ncenters,charge,
     &                              max_centers)
c---------------------------------------------------------------------------
c   "Work" routine for the integral worker task.
c---------------------------------------------------------------------------
      implicit none
c
      integer a1, a2, b1, b2, max_centers 
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
      double precision xg,yg,zg
      double precision goc_x,goc_y,goc_z
      double precision xn(max_centers), yn(max_centers), 
     *                 zn(max_centers)
      double precision Ncharge

cPV we need this for NAI Lo, L-1, L+1
      integer ncenters, component 
      integer iatom
      integer DERCX,DERCY,DERCZ
      integer DER1X, DER1Y,DER1Z,IXDERC,DER2X,DER2Y,DER2Z
      integer Ea1,Ea2,Eb1,Eb2
c end PV
      double precision coords(3,*), coeffs(*), alphas(*),
     *                 charge(*)
      double precision out_dsh_xx(a1:a2,b1:b2),
     *                 out_dsh_yx(a1:a2,b1:b2),
     *                 out_dsh_zx(a1:a2,b1:b2),
     *                 out_dsh_xy(a1:a2,b1:b2),
     *                 out_dsh_yy(a1:a2,b1:b2),
     *                 out_dsh_zy(a1:a2,b1:b2),
     *                 out_dsh_xz(a1:a2,b1:b2),
     *                 out_dsh_yz(a1:a2,b1:b2),
     *                 out_dsh_zz(a1:a2,b1:b2)
   
      double precision scr(*), sum
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
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "The atom label for each shell and coordinates"
      Write(6,"(4(1x,I4))") (Atm_4shell(i), i=1, Nshells)
      Write(6,*)
         Do Iatm =1, Ncenters
         Write(6, "(6(1x,F12.8))") (Coords(i,Iatm), i=1,3)
      Enddo
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
C Xg, Yg, and Zg are the gauge origin (currently set to the origin of
C Cartesian coordinate system.
C
      xg = 0.0D0
      yg = 0.0D0
      zg = 0.0D0

      do b = 1, Nbasis
         do a = 1, Nbasis
            out_dsh_xx(a,b) = 0.0d0
            out_dsh_yx(a,b) = 0.0d0
            out_dsh_zx(a,b) = 0.0d0
            out_dsh_xy(a,b) = 0.0d0
            out_dsh_yy(a,b) = 0.0d0
            out_dsh_zy(a,b) = 0.0d0
            out_dsh_xz(a,b) = 0.0d0
            out_dsh_yz(a,b) = 0.0d0
            out_dsh_zz(a,b) = 0.0d0
         end do
      end do
c-------------------------------------------------------------------------
c   Build an array of coordinates per atom.
c-------------------------------------------------------------------------

      do m = 1, nshells
         iatom = Atm_4shell(m)
         xn(iatom) = coords(1,iatom)
         yn(iatom) = coords(2,iatom)
         zn(iatom) = coords(3,iatom)

#ifdef _DEBUG_LVL0
         Write(6,*)
         Write(6,"(a)") "The center coordinates"
         write(6,"(1x,I3,3F12.6)") iatom, xn(iatom),yn(iatom),
     &                             zn(iatom)
#endif
      enddo
C
C Loop over the centers and compute the electric field integral x, y and z
C component for each center.
C
      do iatm =1, Ncenters
CSSS      do iatm =1, 1

C         do b = 1, Nbasis
C         do a = 1, Nbasis
C            out_dsh_xx(a,b) = 0.0d0
C            out_dsh_xy(a,b) = 0.0d0
C            out_dsh_xz(a,b) = 0.0d0
C            out_dsh_yx(a,b) = 0.0d0
C            out_dsh_yy(a,b) = 0.0d0
C            out_dsh_yx(a,b) = 0.0d0
C            out_dsh_zx(a,b) = 0.0d0
C            out_dsh_zy(a,b) = 0.0d0
C            out_dsh_zz(a,b) = 0.0d0
C         end do
C         end do

         ixderc = iatm
         ncharge = charge(iatm)
C  
C The X component of the electric field; {<mu|derx_C(y-C)|nu> 
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
CSS#ifdef _NOSKIP
           DERCX = 1
           DERCY = 0
           DERCZ = 0
           DER1X = 0
           DER1Y = 0
           DER1Z = 0
           DER2X = 0
           DER2Y = 0
           DER2Z = 0

           call oed__gener_efield_batch(intmax,zmax,
     *                 nalpha,npcoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 ncenters,iatm,xn,yn,zn,charge,ixderc,
     *                 der1x, der1y, der1z,
     *                 der2x, der2y, der2z,
     *                 dercx, dercy, dercz,
     *                 alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .false., iscr, nints, nfirst,
     *                 scr)
            
           If(nints .gt. 0) then
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
#ifdef _DEBUG_LVL0
      Write(6,*) "Efield block"
      Write(6,"(4(1x,i4))") a1,a2,b1,b2
      Write(6,"(4(1x,i4))") Ea1,Ea2,Eb1,Eb2
#endif
C
                 goc_x = (x1 - xg)*Ncharge
                 goc_y = (y1 - yg)*Ncharge
                 goc_z = (z1 - zg)*Ncharge

C                 goc_x = (x1 - xg)
C                 goc_y = (y1 - yg)
C                 goc_z = (z1 - zg)

                 call add_integrals2(out_dsh_xx,a1,a2,b1,b2,
     *                         scr(nfirst),Ea1,Ea2,Eb1,Eb2,
     *                         -goc_x)

                 call add_integrals2(out_dsh_yx,a1,a2,b1,b2,
     *                         scr(nfirst),Ea1,Ea2,Eb1,Eb2,
     *                          -goc_y)

                 call add_integrals2(out_dsh_zx,a1,a2,b1,b2,
     *                         scr(nfirst),Ea1,Ea2,Eb1,Eb2,
     *                         -goc_z)
C
           endif ! nints 
CSSS       Write(6,*) "intermxx", matm, natm
CSSS       call output(out_dsh_xx, 1, nbasis, 1, nbasis, nbasis, nbasis,1)
C
C
CSSS#endif
C
C loop over x, y and z component of center A (ket) derivatives 
C
CSSS#ifdef _NOSKIP
            Do component = 1, 3
              
CSSS               Write(6,"(a,3(1x,i3))") "Component,SHAB:", 
CSSS     &                               component,ivangmom(m),ivangmom(n)
               DERCX = 1
               DERCY = 0
               DERCZ = 0
               DER2X = 0
               DER2Y = 0
               DER2Z = 0

               if (component .EQ. 1) Then
                  DER1X = 1
                  DER1Y = 0
                  DER1Z = 0
               else if (component .EQ. 2) then
                  DER1X = 0
                  DER1Y = 1
                  DER1Z = 0
               else if (component .EQ. 3) then
                  DER1X = 0
                  DER1Y = 0
                  DER1Z = 1
               endif

               call oed__gener_efield_batch(intmax,zmax,
     *                 nalpha,npcoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 ncenters,iatm,xn,yn,zn,charge,ixderc,
     *                 der1x, der1y, der1z,
     *                 der2x, der2y, der2z,
     *                 dercx, dercy, dercz,
     *                 alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .false., iscr, nints, nfirst,
     *                 scr)

               if (nints .gt. 0) then
                  Ea2 = end_nfps(m)
                  if (m .eq. 1) then
                     Ea1 = 1
                  else
                     Ea1 = end_nfps(m-1)+1
                  endif

                  Eb2 = end_nfps(n)
                  if (n .eq. 1) then
                     Eb1 = 1
                  else
                     Eb1 = end_nfps(n-1)+1
                  endif
#ifdef _DEBUG_LVL0
      Write(6,*) "dshield block"
      Write(6,"(4(1x,i4))") a1,a2,b1,b2
      Write(6,"(4(1x,i4))") Ea1,Ea2,Eb1,Eb2
#endif
C
                  If (component .EQ. 1) call add_integrals2(
     *                                    out_dsh_xx,
     *                                    a1,a2,b1,b2,
     *                                    scr(nfirst),Ea1,Ea2,Eb1,Eb2,
     *                                    -Ncharge)

                  If (component .EQ. 2) call add_integrals2(
     *                                    out_dsh_yx,
     *                                    a1,a2,b1,b2,
     *                                    scr(nfirst),Ea1,Ea2,Eb1,Eb2,
     *                                    -Ncharge)

                  If (component .EQ. 3) call add_integrals2(
     *                                    out_dsh_zx,
     *                                    a1,a2,b1,b2,
     *                                    scr(nfirst),Ea1,Ea2,Eb1,Eb2,
     *                                    -Ncharge)
               endif

            end do 
CSSS#endif 
         end do !n
         end do !m

#ifdef _DEBUG_LVL0
      Write(6,*) "The <mu| dsh_xx |nu > integral matrix"
      call output(out_dsh_xx, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
      Write(6,*) "The <mu| dsh_yx |nu > integral matrix"
      call output(out_dsh_yx, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
      Write(6,*) "The <mu| dsh_zx |nu > integral matrix"
      call output(out_dsh_zx, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
C      do b = 1, Nbasis
C         Write(6, "(5(1x,F15.10))") (out_dsh_xx(a,b), a=1, nbasis)
C      enddo
#endif
CSSS#ifdef _NOSKIP
C
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C The Y component of the electric field; {<mu|dery_C(y-C)|nu>
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

            DERCX = 0
            DERCY = 1
            DERCZ = 0
            DER1X = 0
            DER1Y = 0
            DER1Z = 0
            DER2X = 0
            DER2Y = 0
            DER2Z = 0
 
            call oed__gener_efield_batch(intmax,zmax,
     *                 nalpha,npcoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 ncenters,iatm,xn,yn,zn,charge,ixderc,
     *                 der1x, der1y, der1z,
     *                 der2x, der2y, der2z,
     *                 dercx, dercy, dercz,
     *                 dercx, dercy, dercz,
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
            endif
#ifdef _DEBUG_LVLM
      Write(6,*) "yy dshield block-1"
      Write(6,"(4(1x,i4))") a1,a2,b1,b2
      Write(6,"(4(1x,i4))") Ea1,Ea2,Eb1,Eb2
#endif
                 goc_x = (x1 - xg)*Ncharge
                 goc_y = (y1 - yg)*Ncharge
                 goc_z = (z1 - zg)*Ncharge

C                 goc_x = (x1 - xg)
C                 goc_y = (y1 - yg)
C                 goc_z = (z1 - zg)
C
                 call add_integrals2(out_dsh_xy,a1,a2,b1,b2,
     *                               scr(nfirst),Ea1,Ea2,Eb1,Eb2,
     *                               -goc_x)

                 call add_integrals2(out_dsh_yy,a1,a2,b1,b2,
     *                               scr(nfirst),Ea1,Ea2,Eb1,Eb2,
     *                               -goc_y)

                 call add_integrals2(out_dsh_zy,a1,a2,b1,b2,
     *                               scr(nfirst),Ea1,Ea2,Eb1,Eb2,
     *                               -goc_z)
#ifdef _DEBUG_LVLM
      Write(6,*) "yy dshield block-2"
      Write(6,"(4(1x,i4))") a1,a2,b1,b2
      Write(6,"(4(1x,i4))") Ea1,Ea2,Eb1,Eb2
#endif
C
C loop over x, y and z component of center A (ket) derivatives
C
            Do component = 1, 3
    
               DERCX = 0
               DERCY = 1
               DERCZ = 0
               DER2X = 0
               DER2Y = 0
               DER2Z = 0

               if (component .EQ. 1) Then
                  DER1X = 1
                  DER1Y = 0
                  DER1Z = 0
               else if (component .EQ. 2) then
                  DER1X = 0
                  DER1Y = 1
                  DER1Z = 0
               else if (component .EQ. 3) then
                  DER1X = 0
                  DER1Y = 0
                  DER1Z = 1
               endif

               call oed__gener_efield_batch(intmax,zmax,
     *                 nalpha,npcoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 ncenters,iatm,xn,yn,zn,charge,ixderc,
     *                 der1x, der1y, der1z,
     *                 der2x, der2y, der2z,
     *                 dercx, dercy, dercz,
     *                 alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .false., iscr, nints, nfirst,
     *                 scr)

               if (nints .gt. 0) then
                  Ea2 = end_nfps(m)
                  if (m .eq. 1) then
                     Ea1 = 1
                  else
                     Ea1 = end_nfps(m-1)+1
                  endif

                  Eb2 = end_nfps(n)
                  if (n .eq. 1) then
                     Eb1 = 1
                  else
                     Eb1 = end_nfps(n-1)+1
                  endif
               endif
#ifdef _DEBUG_LVLM
      Write(6,*) "yy dshield block-3"
      Write(6,"(4(1x,i4))") a1,a2,b1,b2
      Write(6,"(4(1x,i4))") Ea1,Ea2,Eb1,Eb2
#endif
C
               If (component .EQ. 1) call add_integrals2(out_dsh_xy,
     *                                    a1,a2,b1,b2,
     *                                    scr(nfirst),Ea1,Ea2,Eb1,Eb2,
     *                                    -Ncharge)

               If (component .EQ. 2) call add_integrals2(out_dsh_yy,
     *                                    a1,a2,b1,b2,
     *                                    scr(nfirst),Ea1,Ea2,Eb1,Eb2,
     *                                    -Ncharge)

               If (component .EQ. 3) call add_integrals2(out_dsh_zy,
     *                                    a1,a2,b1,b2,
     *                                    scr(nfirst),Ea1,Ea2,Eb1,Eb2,
     *                                    -Ncharge)
            enddo

         end do !n
         end do !m
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The <mu| dsh_xy |nu > integral matrix"
      call output(out_dsh_xy, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
      Write(6,*) "The <mu| dsh_yy |nu > integral matrix"
      call output(out_dsh_yy, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
      Write(6,*) "The <mu| dsh_zy |nu > integral matrix"
      call output(out_dsh_zy, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
C      do b = 1, Nbasis
C      Write(6,*)
C         CWrite(6, "(5(1x,F15.10))") (out_dsh_xy(a,b), a=1, nbasis)
C      enddo
#endif
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C The Z component of the electric field integral; -{<mu|derz_C(y-C)|nu> -
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

            DERCX = 0
            DERCY = 0
            DERCZ = 1
            DER1X = 0
            DER1Y = 0
            DER1Z = 0
            DER2X = 0
            DER2Y = 0
            DER2Z = 0

            call oed__gener_efield_batch(intmax,zmax,
     *                 nalpha,npcoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 ncenters,iatm,xn,yn,zn,charge,ixderc,
     *                 der1x, der1y, der1z,
     *                 der2x, der2y, der2z,
     *                 dercx, dercy, dercz,
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
            endif ! nints 

            goc_x = (x1 - xg)*Ncharge
            goc_y = (y1 - yg)*Ncharge
            goc_z = (z1 - zg)*Ncharge

C            goc_x = (x1 - xg)
C            goc_y = (y1 - yg)
C            goc_z = (z1 - zg)

C
c       Write(6,*) "intermzz", m ,n , x1, y1, z1, scr(nfirst)
C       call output(out_dsh_zz, 1, nbasis, 1, nbasis, nbasis, nbasis,1)

                 call add_integrals2(out_dsh_xz,a1,a2,b1,b2,
     *                               scr(nfirst),Ea1,Ea2,Eb1,Eb2,
     *                               -goc_x)

                 call add_integrals2(out_dsh_yz,a1,a2,b1,b2,
     *                               scr(nfirst),Ea1,Ea2,Eb1,Eb2,
     *                               -goc_y)

                 call add_integrals2(out_dsh_zz,a1,a2,b1,b2,
     *                               scr(nfirst),Ea1,Ea2,Eb1,Eb2,
     *                              -goc_z)
c       Write(6,*) "intermzz", m ,n 
c       call output(out_dsh_zz, 1, nbasis, 1, nbasis, nbasis, nbasis,1)
C
C loop over x, y and z component of center A (ket) derivatives
C
            Do component = 3, 3

               DERCX = 0
               DERCY = 0
               DERCZ = 1
               DER2X = 0
               DER2Y = 0
               DER2Z = 0
 
               if (component .EQ. 1) Then
                  DER1X = 1
                  DER1Y = 0
                  DER1Z = 0
               else if (component .EQ. 2) then
                  DER1X = 0
                  DER1Y = 1
                  DER1Z = 0
               else if (component .EQ. 3) then
                  DER1X = 0
                  DER1Y = 0
                  DER1Z = 1
               endif

               call oed__gener_efield_batch(intmax,zmax,
     *                 nalpha,npcoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 ncenters,iatm,xn,yn,zn,charge,ixderc,
     *                 der1x, der1y, der1z,
     *                 der2x, der2y, der2z,
     *                 dercx, dercy, dercz,
     *                 alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .false., iscr, nints, nfirst,
     *                 scr)

               if (nints .gt. 0) then
                  Ea2 = end_nfps(m)
                  if (m .eq. 1) then
                     Ea1 = 1
                  else
                     Ea1 = end_nfps(m-1)+1
                  endif

                  Eb2 = end_nfps(n)
                  if (n .eq. 1) then
                     Eb1 = 1
                  else
                     Eb1 = end_nfps(n-1)+1
                  endif
               endif
CSSS       Write(6,*) "intermzz", m ,n , x1, y1, z1, scr(nfirst)
C
               If (component .EQ. 1) call add_integrals2(out_dsh_xz,
     *                                    a1,a2,b1,b2,
     *                                    scr(nfirst),Ea1,Ea2,Eb1,Eb2,
     *                                    -Ncharge)

               If (component .EQ. 2) call add_integrals2(out_dsh_yz,
     *                                    a1,a2,b1,b2,
     *                                    scr(nfirst),Ea1,Ea2,Eb1,Eb2,
     *                                    -Ncharge)

               If (component .EQ. 3) call add_integrals2(out_dsh_zz,
     *                                    a1,a2,b1,b2,
     *                                    scr(nfirst),Ea1,Ea2,Eb1,Eb2,
     *                                    -Ncharge)
            enddo

         end do !n
         end do !m
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The <mu| dsh_xz |nu > integral matrix"
      call output(out_dsh_zx, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
      Write(6,*) "The <mu| dsh_yz |nu > integral matrix"
      call output(out_dsh_zy, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
      Write(6,*) "The <mu| dsh_zz |nu > integral matrix"
      call output(out_dsh_zz, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
C      Write(6,*)
C      Write(6, "(5(1x,F15.10))") ((out_dsh_zz(a,b),a=1,b),b=1,nbasis)
C      do b = 1, Nbasis
C      Write(6,*) 
C         Write(6, "(5(1x,F15.10))") (out_dsh_xz(a,b), a=1, nbasis)
C      enddo
#endif
CSSS#endif      
      enddo ! ncenters
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The integrals for the total atoms" 

      Write(6,*) "The <mu| dsh_xx |nu > integral matrix"
      call output(out_dsh_xx, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
      Write(6,*) "The <mu| dsh_yx |nu > integral matrix"
      call output(out_dsh_yx, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
      Write(6,*) "The <mu| dsh_zx |nu > integral matrix"
      call output(out_dsh_zx, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)

      Write(6,*) "The <mu| dsh_xy |nu > integral matrix"
      call output(out_dsh_xy, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
      Write(6,*) "The <mu| dsh_yy |nu > integral matrix"
      call output(out_dsh_yy, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
      Write(6,*) "The <mu| dsh_zy |nu > integral matrix"
      call output(out_dsh_zy, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)

      Write(6,*) "The <mu| dsh_xz |nu > integral matrix"
      call output(out_dsh_xz, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
      Write(6,*) "The <mu| dsh_yz |nu > integral matrix"
      call output(out_dsh_yz, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
      Write(6,*) "The <mu| dsh_zz |nu > integral matrix"
      call output(out_dsh_zz, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)

#endif


      return
      end

