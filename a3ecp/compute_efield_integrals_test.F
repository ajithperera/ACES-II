      subroutine compute_efield_integrals_test
     &                              (a1,a2,b1,b2,scr,
     &                              maxblk,iscr,coords,coeffs,
     &                              alphas, ccbeg, 
     &                              ccend, end_nfps, nshells, nbasis,
     &                              zmax, intmax, npfps, ncfps, 
     &                              ivangmom,
     &                              ixalpha, indx_cc, Atm_4shell,
     &                              ixpcoef,ispherical, out_efx,
     &                              out_efy,out_efz, nsend, 
     &                              Nalpha, Npcoef,ncenters,charge,
     &                              max_centers)
c---------------------------------------------------------------------------
c   "Work" routine for the integral worker task.
c---------------------------------------------------------------------------
      implicit none
c
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
      double precision xn(max_centers),yn(max_centers),
     *                 zn(max_centers)

cPV we need this for NAI Lo, L-1, L+1
      integer ncenters
      integer iatom
      integer DERCX,DERCY,DERCZ
      integer DER1X, DER1Y,DER1Z,IXDERC,DER2X,DER2Y,DER2Z
      integer Ea1,Ea2,Eb1,Eb2
c end PV
      double precision coords(3,*), coeffs(*), alphas(*),
     *                 charge(*)
      double precision out_efx(a1:a2,b1:b2),
     *                 out_efy(a1:a2,b1:b2),out_efz(a1:a2,b1:b2)
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
C
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
c----------------------------------------------------------------------
      call lookup_shell(end_nfps, nshells, a1, m1)
      call lookup_shell(end_nfps, nshells, a2, m2)
      call lookup_shell(end_nfps, nshells, b1, n1)
      call lookup_shell(end_nfps, nshells, b2, n2)

      do m=1,nshells
         iatom=Atm_4shell(m)
         xn(iatom) = coords(1,iatom)
         yn(iatom) = coords(2,iatom)
         zn(iatom) = coords(3,iatom)
#ifdef _DEBUG_LVL0
         Write(6,*)
         Write(6,"(a)") "The center coordinates"
         write(6,"(1x,I3,3F12.6)") iatom, xn(iatom),yn(iatom),
     &                             zn(iatom)
#endif
      end do
C
C Loop over the centers and compute the electric field integral x, y and z
C component for each center.
C
      do iatm=1, Ncenters

         do b = 1, Nbasis
         do a = 1, Nbasis
            out_efx(a,b) = 0.d0
            out_efy(a,b) = 0.d0
            out_efz(a,b) = 0.d0
         end do
         end do
          
         xn = Coords(1,iatm)
         yn = Coords(2,iatm)
         zn = Coords(3,iatm)

         ixderc = Iatm
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

CSSS      Write(6,*) "@-entry to OED_batch"
CSSS      Write(6,"((a,2(1x,i4)))") "m, n :", m, n
           
           DERCX = 0
           DERCY = 0
           DERCZ = 0
           DER1X = 1
           DER1Y = 0
           DER1Z = 0
           DER2X = 0
           DER2Y = 0
           DER2Z = 0
           Write(6,*) "The first call"

           call oed__gener_nai_derv_batch_test(intmax,zmax,
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
C
              call add_integrals2(out_efx,a1,a2,b1,b2,
     *                            scr(nfirst),Ea1,Ea2,Eb1,Eb2,1.0d0)
            endif ! nints 
      call output(out_efx, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)

           DERCX = 0
           DERCY = 0
           DERCZ = 0
           DER1X = 0
           DER1Y = 0
           DER1Z = 0
           DER2X = 1
           DER2Y = 0
           DER2Z = 0

           Write(6,*) "The second call"
           call oed__gener_nai_derv_batch_test(intmax,zmax,
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
C
              call add_integrals2(out_efx,a1,a2,b1,b2,
     *                            scr(nfirst),Ea1,Ea2,Eb1,Eb2,1.0d0)
            endif ! nints


      call output(out_efx, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
         end do !n
         end do !m
C
#ifdef _DEBUG_LVL0
C      Write(6,*) "The <mu| ef_X |nu > integral matrix"
C      call output(out_efx, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
C      do b = 1, Nbasis
C      Write(6,*)
C         Write(6, "(5(1x,F15.10))") (out_efx(a,b), a=1, nbasis)
C      enddo
#endif
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
            DERCY = 0
            DERCZ = 0
            DER1X = 0
            DER1Y = 1
            DER1Z = 0
            DER2X = 0
            DER2Y = 0
            DER2Z = 0

            call oed__gener_nai_derv_batch_test(intmax,zmax,
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
C
              call add_integrals2(out_efy,a1,a2,b1,b2,
     *                            scr(nfirst),Ea1,Ea2,Eb1,Eb2,1.0d0)
            endif ! nints 

            DERCX = 0
            DERCY = 0
            DERCZ = 0
            DER1X = 0
            DER1Y = 0
            DER1Z = 0
            DER2X = 0
            DER2Y = 1
            DER2Z = 0

            call oed__gener_nai_derv_batch_test(intmax,zmax,
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
C
              call add_integrals2(out_efy,a1,a2,b1,b2,
     *                            scr(nfirst),Ea1,Ea2,Eb1,Eb2,1.0d0)
            endif ! nints



         end do !n
         end do !m
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The <mu| ef_Y |nu > integral matrix"
      call output(out_efy, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
C      do b = 1, Nbasis
C      Write(6,*)
C         Write(6, "(5(1x,F15.10))") (out_efy(a,b), a=1, nbasis)
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
            DERCZ = 0
            DER1X = 0
            DER1Y = 0
            DER1Z = 1
            DER2X = 0
            DER2Y = 0
            DER2Z = 0

            call oed__gener_nai_derv_batch_test(intmax,zmax,
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
C
              call add_integrals2(out_efz,a1,a2,b1,b2,
     *                            scr(nfirst),Ea1,Ea2,Eb1,Eb2,1.0d0)
            endif ! nints 
            DERCX = 0
            DERCY = 0
            DERCZ = 0
            DER1X = 0
            DER1Y = 0
            DER1Z = 0
            DER2X = 0
            DER2Y = 0
            DER2Z = 1 

            call oed__gener_nai_derv_batch_test(intmax,zmax,
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
C
              call add_integrals2(out_efz,a1,a2,b1,b2,
     *                            scr(nfirst),Ea1,Ea2,Eb1,Eb2,1.0d0)
            endif ! nints


         end do !n
         end do !m
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The <mu| ef_Z |nu > integral matrix"
      call output(out_efz, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
C      do b = 1, Nbasis
C      Write(6,*) 
C         Write(6, "(5(1x,F15.10))") (out_efz(a,b), a=1, nbasis)
C      enddo
#endif
     
      enddo ! ncenters
C
      return
      end

