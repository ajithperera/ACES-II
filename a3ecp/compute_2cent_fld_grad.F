      subroutine compute_2cent_fld_grad
     &                              (a1,a2,b1,b2,scr,
     &                              maxblk,iscr,coords,coeffs,
     &                              alphas, ccbeg, 
     &                              ccend, end_nfps, nshells, nbasis,
     &                              zmax, intmax, npfps, ncfps, 
     &                              ivangmom,
     &                              ixalpha, indx_cc, Atm_4shell,
     &                              ixpcoef,ispherical, 
     &                              out_fldgr_xx,out_fldgr_xy,
     &                              out_fldgr_xz,
     &                              out_fldgr_yx,out_fldgr_yy,
     &                              out_fldgr_yz,
     &                              out_fldgr_zx,out_fldgr_zy,
     &                              out_fldgr_zz,
     &                              nsend, 
     &                              Nalpha, Npcoef,ncenters,charge,
     &                              max_centers,ERD_index,
     &                              ERD_scale,ACES_4FLDGR)
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
      integer a,b,c,d,ixyz,jxyz 
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
      integer  ERD_index(*) 

      logical skip
      logical mn_symmetry
      logical rs_symmetry
      logical mn_rs_symmetry
      logical*8 l8true, l8spherical
      logical spherical, erfints

      double precision x1,y1,z1
      double precision x2,y2,z2
      double precision xg,yg,zg
      double precision ncharge, erd_scale(*) 
      
      double precision xn(max_centers), yn(max_centers), 
     *                 zn(max_centers)

      integer ncenters, component 
      integer iatom
      integer DERCX,DERCY,DERCZ
      integer DER1X, DER1Y,DER1Z,IXDERC,DER2X,DER2Y,DER2Z
      integer Ea1,Ea2,Eb1,Eb2
      integer XX, XY, XZ, YY, YZ, ZZ

      double precision coords(3,*), coeffs(*), alphas(*),
     *                 charge(*)

      double precision out_fldgr_xx(a1:a2,b1:b2),
     *                 out_fldgr_yx(a1:a2,b1:b2),
     *                 out_fldgr_zx(a1:a2,b1:b2),
     *                 out_fldgr_xy(a1:a2,b1:b2),
     *                 out_fldgr_yy(a1:a2,b1:b2),
     *                 out_fldgr_zy(a1:a2,b1:b2),
     *                 out_fldgr_xz(a1:a2,b1:b2),
     *                 out_fldgr_yz(a1:a2,b1:b2),
     *                 out_fldgr_zz(a1:a2,b1:b2),
     *                 aces_4fldgr(a1:a2,b1:b2)
   
      double precision scr(*), sum, scale_factor 
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

         Call Dzero(out_fldgr_xx, Nbasis * Nbasis)
         Call Dzero(out_fldgr_xy, Nbasis * Nbasis)
         Call Dzero(out_fldgr_xz, Nbasis * Nbasis)
         Call Dzero(out_fldgr_yy, Nbasis * Nbasis)
         Call Dzero(out_fldgr_yz, Nbasis * Nbasis)
         Call Dzero(out_fldgr_zz, Nbasis * Nbasis)
C
         DER1X = 0
         DER1Y = 0
         DER1Z = 0
         DER2X = 0
         DER2Y = 0
         DER2Z = 0

      Do Ixyz = 1, 3
         Do Jxyz = Ixyz+1, 3

         If (Ixyz .eq. 1 .and. jxyz .eq. 2) then
            DERCX = 1
            DERCY = 1
            DERCZ = 0
            XY    = 1
            XZ    = 0
            YZ    = 0
         Elseif (Ixyz .eq. 1 .and. jxyz .eq. 3) then
            DERCX = 1
            DERCY = 0
            DERCZ = 1
            XY    = 0
            XZ    = 1
            YZ    = 0
         Elseif (Ixyz .eq. 2 .and. jxyz .eq. 3) then
            DERCX = 0
            DERCY = 1
            DERCZ = 1
            XY    = 0
            XZ    = 0
            YZ    = 1
         Endif 

         ixderc = iatm
         ncharge = charge(iatm)

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
            
           If (nints .gt. 0) then
C
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
                 Write(6,"(a,3i4)") "XY,XZ,YZ:", XY, XZ, YZ
C
                 If (XY .eq. 1) call add_integrals2(out_fldgr_xy,a1,a2,
     *                          b1,b2,scr(nfirst),Ea1,Ea2,Eb1,Eb2,1.0D0)

                 If (XZ .eq. 1) call add_integrals2(out_fldgr_xz,a1,a2,
     *                          b1,b2,scr(nfirst),Ea1,Ea2,Eb1,Eb2,1.0D0)

                 If (YZ .eq. 1) call add_integrals2(out_fldgr_yz,a1,a2,
     *                          b1,b2,scr(nfirst),Ea1,Ea2,Eb1,Eb2,1.0D0)

           Endif 
C 
         end do 
         end do 
         end do 
         end do 
CSSS#ifdef _NOSKIP
C
C The digonal components 
C
         do ixyz = 1, 3 
   
         If (Ixyz .Eq. 1) then
            DERCX = 2
            DERCY = 0
            DERCZ = 0
            XX    = 1
            YY    = 0
            ZZ    = 0 
         Else if (Ixyz .eq. 2) then
            DERCX = 0
            DERCY = 2
            DERCZ = 0
            XX    = 0
            YY    = 1
            ZZ    = 0
         Else if (Ixyz .eq. 3) then
            DERCX = 0
            DERCY = 0
            DERCZ = 2
            XX    = 0
            YY    = 0
            ZZ    = 1
         Endif 
           
         ixderc = iatm
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

           If (nints .gt. 0) then
C
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

           Write(6,"(a,3i4)") "XX,YY,ZZ:", XX,YY,ZZ 

           If (XX .eq. 1) call add_integrals2(out_fldgr_xx,a1,a2,
     *                    b1,b2,scr(nfirst),Ea1,Ea2,Eb1,Eb2,2.0D0)

           If (YY .eq. 1) call add_integrals2(out_fldgr_yy,a1,a2,
     *                    b1,b2,scr(nfirst),Ea1,Ea2,Eb1,Eb2,2.0D0)

           If (ZZ .eq. 1) call add_integrals2(out_fldgr_zz,a1,a2,
     *                    b1,b2,scr(nfirst),Ea1,Ea2,Eb1,Eb2,2.0D0)
     
           Endif
       enddo
       enddo 

       If (Ixyz .Eq. 1) then
            DERCX = 0
            DERCY = 2
            DERCZ = 0
            XX    = 1
            YY    = 0
            ZZ    = 0
       Else if (Ixyz .eq. 2) then
            DERCX = 2
            DERCY = 0
            DERCZ = 0
            XX    = 0
            YY    = 1
            ZZ    = 0
      Else if (Ixyz .eq. 3) then
            DERCX = 2
            DERCY = 0
            DERCZ = 0
            XX    = 0
            YY    = 0
            ZZ    = 1
      Endif

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

           If (nints .gt. 0) then
C
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

           Write(6,"(a,3i4)") "XX,YY,ZZ:", XX,YY,ZZ 

           If (XX .eq. 1) call add_integrals2(out_fldgr_xx,a1,a2,
     *                    b1,b2,scr(nfirst),Ea1,Ea2,Eb1,Eb2,-1.0D0)

           If (YY .eq. 1) call add_integrals2(out_fldgr_yy,a1,a2,
     *                    b1,b2,scr(nfirst),Ea1,Ea2,Eb1,Eb2,-1.0D0)

           If (ZZ .eq. 1) call add_integrals2(out_fldgr_zz,a1,a2,
     *                    b1,b2,scr(nfirst),Ea1,Ea2,Eb1,Eb2,-1.0D0)
           Endif

       enddo
       enddo

       If (Ixyz .Eq. 1) then
            DERCX = 0
            DERCY = 0
            DERCZ = 2
            XX    = 1
            YY    = 0
            ZZ    = 0
       Else if (Ixyz .eq. 2) then
            DERCX = 0
            DERCY = 0
            DERCZ = 2
            XX    = 0
            YY    = 1
            ZZ    = 0
      Else if (Ixyz .eq. 3) then
            DERCX = 0
            DERCY = 2
            DERCZ = 0
            XX    = 0
            YY    = 0
            ZZ    = 1
      Endif

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

           If (nints .gt. 0) then
C
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

           Write(6,"(a,3i4)") "XX,YY,ZZ:", XX,YY,ZZ 

           If (XX .eq. 1) call add_integrals2(out_fldgr_xx,a1,a2,
     *                    b1,b2,scr(nfirst),Ea1,Ea2,Eb1,Eb2,-1.0D0)

           If (YY .eq. 1) call add_integrals2(out_fldgr_yy,a1,a2,
     *                    b1,b2,scr(nfirst),Ea1,Ea2,Eb1,Eb2,-1.0D0)

           If (ZZ .eq. 1) call add_integrals2(out_fldgr_zz,a1,a2,
     *                    b1,b2,scr(nfirst),Ea1,Ea2,Eb1,Eb2,-1.0D0)
           Endif

       enddo
       enddo

       enddo

       Scale_factor = 1.0D0/3.0D0
       Call Dscal(Nbasis * Nbasis, scale_factor, out_fldgr_xx, 1) 
       Call Dscal(Nbasis * Nbasis, scale_factor, out_fldgr_yy, 1) 
       Call Dscal(Nbasis * Nbasis, scale_factor, out_fldgr_zz, 1) 

CSSS#endif 
C
C Since there is only one location we can not check all components at 
C once by comparing to ACES II. Do it one at a time

#ifdef _DEBUG_LVLM
      If (XX .eq. 1)  Call rescale_ints(out_fldgr_xx, aces_4fldgr,Nbasis,
     &                                  ERD_index,ERD_scale)
#endif 
C
#ifdef _DEBUG_LVL0

      Write(6,*) "The <mu| fldgr_xy |nu > integral matrix"
      call output(out_fldgr_xy, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)

      Write(6,*) "The <mu| fldgr_xz |nu > integral matrix"
      call output(out_fldgr_xz, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)

      Write(6,*) "The <mu| fldgr_yz |nu > integral matrix"
      call output(out_fldgr_yz, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)

#endif

#ifdef _DEBUG_LVL0
CSSS#ifdef _NOSKIP
      Write(6,*) "The <mu| fldgr_xx |nu > integral matrix"
      call output(out_fldgr_xx, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)

      Write(6,*) "The <mu| fldgr_yy |nu > integral matrix"
      call output(out_fldgr_yy, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)

      Write(6,*) "The <mu| fldgr_zz |nu > integral matrix"
      call output(out_fldgr_zz, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
CSSS#endif 
#endif
      enddo ! ncenters

      return
      end

