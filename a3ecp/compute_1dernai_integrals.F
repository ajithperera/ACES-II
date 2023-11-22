C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
      subroutine compute_1dernai_integrals(a1,a2,b1,b2,scr,
     &                              maxblk,iscr,coords,coeffs,
     &                              alphas, ccbeg,
     &                              ccend, end_nfps, nshells, nbasis,
     &                              zmax, intmax, npfps, ncfps,
     &                              ivangmom,
     &                              ixalpha, indx_cc, Atm_4shell,
     &                              ixpcoef,ispherical,
     &                              out_nai_x,out_nai_y,out_nai_z,
     &                              out_nai,
     &                              nsend,
     &                              Nalpha, Npcoef,ncenters,charge,
     &                              max_centers,ERD_index,ERD_scale)

c---------------------------------------------------------------------------
c  Computes the set of PVP integrals
c---------------------------------------------------------------------------
      implicit none
      integer a1, a2, b1, b2, max_centers
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
      integer matom, natom, ratm, satm

      integer imin, zmin, iblk, zblk, iatm
      integer ivangmom(*), ixalpha(*), ixpcoef(*), npfps(*)
      integer ncfps(*), indx_cc(*), Atm_4shell(*)
      integer ERD_index(*) 

      logical skip
      logical mn_symmetry
      logical rs_symmetry
      logical mn_rs_symmetry
      logical*8 l8true, l8spherical
      logical spherical, erfints

      double precision x1,y1,z1
      double precision x2,y2,z2
      double precision xg,yg,zg
      double precision Nai,Fact
      double precision xn(max_centers), yn(max_centers),
     *                 zn(max_centers)
      double precision Ncharge, ERD_scale(*)
      integer ncenters, component
      integer iatom
      integer DERCX,DERCY,DERCZ
      integer DER1X, DER1Y,DER1Z,IXDERC,DER2X,DER2Y,DER2Z
      integer Ea1,Ea2,Eb1,Eb2
      integer icomponent,jcomponent

      double precision coords(3,*), coeffs(*), alphas(*),
     *                 charge(*)
      double precision out_nai_x(a1:a2,b1:b2)
      double precision out_nai_y(a1:a2,b1:b2)
      double precision out_nai_z(a1:a2,b1:b2)
      double precision out_nai(a1:a2,b1:b2)

      double precision scr(*), sum
      integer iscr(*)

      integer ccbeg(*), ccend(*), end_nfps(*)

      integer max_dim_coeff
      parameter (max_dim_coeff = 5000)
      integer ccbeg_pack(max_dim_coeff), ccend_pack(max_dim_coeff)
      double precision alpha_pack(max_dim_coeff),
     *                 pcoeff_pack(max_dim_coeff)
      save me,alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack

      spherical = (ispherical .eq. 1)
      l8spherical = spherical
      l8true = .false.

c --------------------------------------------------------------------------- 
c  Set up send data. 
c --------------------------------------------------------------------------- 
      adim = Nbasis 
      bdim = Nbasis 

      nsend = adim*bdim
      if (nsend .lt. 0) then
         print *,'ERROR IN COMP_FOCK_DER ',me,' nsend = ',nsend
         print *,'adim,bdim = ',adim,bdim
         call Errex 
      endif 
      Fact = 1.0D0

#ifdef  _DEBUG_LVLM
        Write(*,"(a)") "Reorder Index array"
        Write(*,"(10I5)") (ERD_Index(i), i=1, Nbasis)
        Write(*,"(a)") "The scale values"
        Write(*,"(1x,10F15.7)") (ERD_scale(i), i=1, Nbasis)
#endif

c-------------------------------------------------------------------------
c   Find the shells we need to loop over for this block.
c--------------------------------------------------------------------------

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
      enddo

#ifdef _DEBUG_LVL0
         Write(6,*)
         Write(6,"(a)") "The center coordinates"
         do iatom = 1, Ncenters
         write(6,"(1x,I3,3F12.6)") iatom, xn(iatom),yn(iatom),
     &                             zn(iatom)
         enddo
#endif
c-------------------------------------------------------------------------
c   Calculate the integrals over the necessary shell blocks.
c-------------------------------------------------------------------------
         do ixderc = 1, ncenters
            Iatom = ixderc

           do icomponent = 1, 3   

            do b = 1, Nbasis
            do a = 1, Nbasis
               Out_nai_x(a,b) = 0.0d0
               Out_nai_y(a,b) = 0.0d0
               Out_nai_z(a,b) = 0.0d0
            end do
            end do

            do m = m1, m2
               matom = Atm_4shell(m)
               x1 = coords(1,matom)
               y1 = coords(2,matom)
               z1 = coords(3,matom)
 
            do n = n1, n2
               natom = Atm_4shell(n)
               x2 = coords(1,natom)
               y2 = coords(2,natom)
               z2 = coords(3,natom)

            call pack_coeffs_oed(alphas, ixalpha, coeffs, ixpcoef,
     *                           ncfps, npfps, m, n,
     *                           alpha_pack, nalpha,
     *                           pcoeff_pack, npcoeff,
     *                           ccbeg, ccend, indx_cc,
     *                           ccbeg_pack, ccend_pack,max_dim_coeff)

            ncsum = ncfps(m) + ncfps(n)

c-------------------------------------------------------------------------
c   Clear the integral block.

            der1x = 0
            der1y = 0
            der1z = 0
            der2x = 0
            der2y = 0
            der2z = 0
            dercx = 0
            dercy = 0
            dercz = 0

            if (icomponent .eq. 1) dercx = 1
            if (icomponent .eq. 2) dercy = 1
            if (icomponent .eq. 3) dercz = 1

            if (ixderc .eq. matom) then
               if (icomponent .eq. 1) der1x = 1
               if (icomponent .eq. 2) der1y = 1
               if (icomponent .eq. 3) der1z = 1
            endif

            if (ixderc .eq. natom) then
               if (icomponent .eq. 1) der2x = 1
               if (icomponent .eq. 2) der2y = 1
               if (icomponent .eq. 3) der2z = 1
            endif

#ifdef _DEBUG_LVLM
      Write(6,*) "@-Entry to compute_1der_nai_integrals"
      Write(6,"((a,4(1x,i4)))") "m, n:", m, n
      Write(6,"(a,2(1x,I4))") "The Nalpha and Npcoeff:",
     &                         nalpha, npcoeff
      Write(6,"((a,2(1x,i4)))") "Ivanagmom:", Ivangmom(m),
     &      Ivangmom(n)
      Write(6,"((a,2(1x,i4)))") "Ncfps:", Ncfps(m),
     &      Ncfps(n)
      Write(6,"((a,2(1x,i4)))") "Npfps:", Npfps(m),
     &      Npfps(n)
      Write(6,*)
      Write(6,"(4(1x,F15.7))") (Alpha_pack(i),i=1,Nalpha)
      Write(6,*)
      Write(6,"(4(1x,F15.7))")(Pcoeff_pack(i),i=1,npcoeff)
      Write(6,*)
      Write(6,"(4(1x,I4))")(CCbeg_pack(i),i=1,Nalpha)
      Write(6,*)
      Write(6,"(4(1x,I4))")(CCend_pack(i),i=1,Nalpha)
      Write(6,*)
      Write(6,"(a,3(1x,I10))") "The intmax,Zmax,ncsum:",intmax, zmax,
     &                       ncsum
      Write(6,*) "The coordiantes"
      Write(6,*)
      Write(6,"(3(1x,F10.5))") x1, y1, z1
      Write(6,"(3(1x,F10.5))") x2, y2, z2
#endif 
            call oed__gener_nai_derv_batch
     *                (intmax, zmax, nalpha, npcoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 ncenters, xn,yn,zn, charge,ixderc,
     *                 der1x, der1y, der1z,
     *                 der2x, der2y, der2z,
     *                 dercx, dercy, dercz,
     *                 alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .false., iscr, nints, nfirst,
     *                 scr)

#ifdef _NOSKIP

           call oed__gener_nai_batch(intmax,zmax,
     *                 nalpha,npcoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 ncenters,xn,yn,zn,charge,
     *                 alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .false., iscr, nints, nfirst,
     *                 scr)
#endif 

#ifdef _DEBUG_LVLM
      Write(6,*) "@-compute_1der_nai_integrals"
      Write(6,"(a,i4)") " Nints: ", Nints 
      Write(6,"(4(1x,i4))") a1,a2,b1,b2
      Write(6,"(4(1x,i4))") Ea1,Ea2,Eb1,Eb2
#endif

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

              If (Icomponent .Eq. 1) Then
              call add_integrals2(out_nai_x,a1,a2,b1,b2,scr(nfirst),Ea1,
     *                            Ea2,Eb1,Eb2,1.0D0)
              Elseif (Icomponent .Eq. 2) Then
              call add_integrals2(out_nai_y,a1,a2,b1,b2,scr(nfirst),Ea1,
     *                            Ea2,Eb1,Eb2,1.0D0)
              Elseif (Icomponent .Eq. 3) Then
              call add_integrals2(out_nai_z,a1,a2,b1,b2,scr(nfirst),Ea1,
     *                            Ea2,Eb1,Eb2,1.0D0)
              Endif 

            endif ! nints 

         enddo
         enddo 

         enddo   ! n shells
         enddo   ! m shells
      

#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6,"(a)") " The d/drc(x)<mu|1/r-Rc|nu> integral"
      call output(out_nai_x, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
      Write(6,"(a)") " The d/drc(y)<mu|1/r-Rc|nu> integral"
      call output(out_nai_y, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
      Write(6,"(a)") " The d/drc(z)<mu|1/r-Rc|nu> integral"
      call output(out_nai_z, 1, nbasis, 1, nbasis, nbasis, nbasis, 1)
#endif 
#ifdef _DEBUG_LVLM
      Nai = 0.0d0
      Call rescale_2ints(out_nai_x,out_nai,Nbasis,Erd_index,ERD_scale)
      do a = 1, Nbasis 
      do b = 1, Nbasis 
         Nai = Nai  + Out_nai(a,b) * Out_nai(a,b)
      enddo
      enddo 
      Write(6,"(a,F20.10)")"Check sum of scaled Nai_x integrals = " ,Nai
      Nai = 0.0d0
      Call rescale_2ints(out_nai_y,out_nai,Nbasis,Erd_index,ERD_scale)
      do a = 1, Nbasis
      do b = 1, Nbasis
         Nai = Nai  + Out_nai(a,b) * Out_nai(a,b)
      enddo
      enddo
      Write(6,"(a,F20.10)")"Check sum of scaled Nai_x integrals = " ,Nai
      Nai = 0.0d0
      Call rescale_2ints(out_nai_z,out_nai,Nbasis,Erd_index,ERD_scale)
      do a = 1, Nbasis
      do b = 1, Nbasis
         Nai = Nai  + Out_nai(a,b) * Out_nai(a,b)
      enddo
      enddo
      Write(6,"(a,F20.10)")"Check sum of scaled Nai_z integrals = " ,Nai
#endif 

      return
      end
