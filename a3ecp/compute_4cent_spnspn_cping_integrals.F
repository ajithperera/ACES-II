      subroutine compute_4cent_spnspn_cpling_integrals
     &                              (a1,a2,b1,b2,c1,c2,d1,d2,scr,
     &                              maxblk,iscr,coords,coeffs,
     &                              alphas, ccbeg, 
     &                              ccend, end_nfps, nshells, nbasis,
     &                              zmax, intmax, npfps, ncfps, 
     &                              ivangmom,
     &                              ixalpha, indx_cc, Atm_4shell,
     &                              ixpcoef,ispherical, 
     &                              out_spnspn_xx,
     &                              out_spnspn_xy,
     &                              out_spnspn_xz,
     &                              out_spnspn_yx,
     &                              out_spnspn_yy,
     &                              out_spnspn_yz,
     &                              out_spnspn_zx,
     &                              out_spnspn_zy,
     &                              out_spnspn_zz, 
     &                              nsend, 
     &                              Nalpha, Npcoef,ncenters,charge,
     &                              max_centers,ERD_index,
     &                              ERD_scale,ACES_4DSH)
c---------------------------------------------------------------------------
c   "Work" routine for the integral worker task.
c---------------------------------------------------------------------------
      implicit none
c
      integer a1, a2, b1, b2, c1, c2, d1, d2
      integer aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2
      integer adim, bdim, cdim, ddim
      integer m1, m2, n1, n2, r1, r2, s1, s2
      integer i, j, n, m, r, s
      integer a,b,c,d
      integer nshells, nbasis, zmax, intmax, ispherical

      integer num_to_do, nsend, max_centers
      integer nints, maxblk
      integer nalpha_pack, npcoeff_pack,npcoeff
      integer ncsum, next, nfirst
      integer me, ierr, nalpha, npcoef, nprim_quad, ncoef_quad
      integer matm, natm, ratm, satm
      integer ixyz,jxyz,nbas4

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
      double precision x3,y3,z3
      double precision x4,y4,z4
      double precision x, y, z
      double precision xg,yg,zg
      double precision goc_x,goc_y,goc_z,s11,s22
      integer xx,xy,xz,yx,yy,yz,zx,zy,zz
      double precision xx_s,yx_s,yy_s,zx_s,zy_s,zz_s
      double precision xy_s,xz_s,yz_s
      double precision xx_l,yx_l,yy_l,zx_l,zy_l,zz_l
      double precision xy_l,xz_l,yz_l

      double precision xn(max_centers), yn(max_centers), 
     *                 zn(max_centers)
      double precision Ncharge, ERD_scale(*)

cPV we need this for NAI Lo, L-1, L+1
      integer ncenters, component 
      integer iatom
      integer Ea1,Ea2,Eb1,Eb2, icount
c end PV
      double precision coords(3,*), coeffs(*), alphas(*),
     *                 charge(*)
      double precision out_spnspn_xx(a1:a2,b1:b2,c1:c2,d1:d2),
     *                 out_spnspn_xy(a1:a2,b1:b2,c1:c2,d1:d2),
     *                 out_spnspn_xz(a1:a2,b1:b2,c1:c2,d1:d2),
     *                 out_spnspn_yx(a1:a2,b1:b2,c1:c2,d1:d2),
     *                 out_spnspn_yy(a1:a2,b1:b2,c1:c2,d1:d2),
     *                 out_spnspn_yz(a1:a2,b1:b2,c1:c2,d1:d2),
     *                 out_spnspn_zx(a1:a2,b1:b2,c1:c2,d1:d2),
     *                 out_spnspn_zy(a1:a2,b1:b2,c1:c2,d1:d2),
     *                 out_spnspn_zz(a1:a2,b1:b2,c1:c2,d1:d2),
     *                 ACES_4DSH(a1:a2,b1:b2,c1:c2,d1:d2)
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
      cdim = Nbasis
      ddim = Nbasis
      l8true = .true.
      spherical = (ispherical .eq. 1)
      l8spherical = spherical

      nsend = adim*bdim*cdim*ddim
      if (nsend .lt. 0) then
         print *,'ERROR IN INTEGRAL WORKER ',me,' nsend = ',nsend
         print *,'adim,bdim,cdim,ddim = ',adim,bdim,cdim,ddim
         call Errex
      endif
c----------------------------------------------------------------------------
c   Clear the output array.
c----------------------------------------------------------------------------

      Write(6,*)
      write(6,*) "Printing from compute_4cent_spnspn_cpling_integrals"
      Write(6,*) "The atom label for each shell and coordinates"
      Write(6,"(4(1x,I4))") (Atm_4shell(i), i=1, Nshells)
      Write(6,*)
         Do Iatm =1, Ncenters
         Write(6, "(6(1x,F12.8))") (Coords(i,Iatm), i=1,3)
      Enddo
      Write(6,*)

c-----------------------------------------------------------------------
c   Find the shell blocks for which we shall loop through.
c-----------------------------------------------------------------------
      call lookup_shell(end_nfps, nshells, a1, m1)
      call lookup_shell(end_nfps, nshells, a2, m2)
      call lookup_shell(end_nfps, nshells, b1, n1)
      call lookup_shell(end_nfps, nshells, b2, n2)
      call lookup_shell(end_nfps, nshells, c1, r1)
      call lookup_shell(end_nfps, nshells, c2, r2)
      call lookup_shell(end_nfps, nshells, d1, s1)
      call lookup_shell(end_nfps, nshells, d2, s2)

C         Write(6,*)
C         Write(6,*) "Loop over shell quads." 
C         Write(6,"(4(a,1x,I4,1x,I4,1x))") "m1 and m2: ", m1, m2,
C     &                                    "n1 and N2: ", n1, n2, 
C     &                                    "r1 and r2: ", r1, r2, 
C     &                                    "s1 and s2: ", s1, s2
C 
C
C Xg, Yg, and Zg are the gauge origin (currently set to the origin of
C Cartesian coordinate system.
C
      xg = 0.0D0
      yg = 0.0D0
      zg = 0.0D0
      do i =1, intmax
         iScr(i) = 0
      enddo
     
      do a = 1, Nbasis
         do b = 1, Nbasis
            do c = 1, Nbasis 
               do d = 1, Nbasis 
                  out_spnspn_xx(a,b,c,d) = 0.0d0
                  out_spnspn_xy(a,b,c,d) = 0.0d0
                  out_spnspn_xz(a,b,c,d) = 0.0d0
                  out_spnspn_yx(a,b,c,d) = 0.0d0
                  out_spnspn_yy(a,b,c,d) = 0.0d0
                  out_spnspn_yz(a,b,c,d) = 0.0d0
                  out_spnspn_zx(a,b,c,d) = 0.0d0
                  out_spnspn_zy(a,b,c,d) = 0.0d0
                  out_spnspn_zz(a,b,c,d) = 0.0d0
               enddo
            enddo 
         end do
      end do


C#ifdef _NOSKIP

      Write(6,*)  "Entering the compute_4cent_dshield_integrals"
      write(*,*)

      Do Ixyz = 1, 3
         Do Jxyz = 1, 3

            If (Ixyz .eq. 1 .and. jxyz .eq. 1) then
               XX=1
               XY=0
               XZ=0
               YX=0
               YY=0
               YZ=0
               ZX=0
               ZY=0
               ZZ=0
            ElseIf (Ixyz .eq. 1 .and. jxyz .eq. 2) then
               XX=0
               XY=1
               XZ=0
               YX=0
               YY=0
               YZ=0
               ZX=0
               ZY=0
               ZZ=0
            Elseif (Ixyz .eq. 1 .and. jxyz .eq. 3) then
               XX=0
               XY=0
               XZ=1
               YX=0
               YY=0
               YZ=0
               ZX=0
               ZY=0
               ZZ=0
            elseif (Ixyz .eq. 2 .and. jxyz .eq. 1) then
               XX=0
               XY=0
               XZ=0
               YX=1
               YY=0
               YZ=0
               ZX=0
               ZY=0
               ZZ=0
            elseif (Ixyz .eq. 2 .and. jxyz .eq. 2) then
               XX=0
               XY=0
               XZ=0
               YX=0
               YY=1
               YZ=0
               ZX=0
               ZY=0
               ZZ=0
            elseif (Ixyz .eq. 2 .and. jxyz .eq. 3) then
               XX=0
               XY=0
               XZ=0
               YX=0
               YY=0
               YZ=1
               ZX=0
               ZY=0
               ZZ=0
            elseif (Ixyz .eq. 3 .and. jxyz .eq. 1) then
               XX=0
               XY=0
               XZ=0
               YX=0
               YY=0
               YZ=0
               ZX=1
               ZY=0
               ZZ=0
            elseif (Ixyz .eq. 3 .and. jxyz .eq. 2) then
               XX=0
               XY=0
               XZ=0
               YX=0
               YY=0
               YZ=0
               ZX=0
               ZY=1
               ZZ=0
            elseif (Ixyz .eq. 3 .and. jxyz .eq. 3) then
               XX=0
               XY=0
               XZ=0
               YX=0
               YY=0
               YZ=0
               ZX=0
               ZY=0
               ZZ=1
            endif
 
c-------------------------------------------------------------------------
c   Build an array of coordinates per atom.
c-------------------------------------------------------------------------

      do m = m1, m2
         if (m .eq. 1) then
             aa1 = 1
         else
             aa1 = end_nfps(m-1) + 1
         endif
         aa2 = end_nfps(m)

         matm = Atm_4shell(m)
         x1 = coords(1,matm)
         y1 = coords(2,matm)
         z1 = coords(3,matm)

         do n = n1, n2
            if (n .eq. 1) then
               bb1 = 1
            else
               bb1 = end_nfps(n-1) + 1
            endif
            bb2 = end_nfps(n)

            natm = Atm_4shell(n)
            x2 = coords(1,natm)
            y2 = coords(2,natm)
            z2 = coords(3,natm)

            do r = r1, r2
              if (r .eq. 1) then
                  cc1 = 1
              else
                  cc1 = end_nfps(r-1) + 1
              endif
              cc2 = end_nfps(r)

              ratm = Atm_4shell(r)
              x3 = coords(1,ratm)
              y3 = coords(2,ratm)
              z3 = coords(3,ratm)

              do s = s1, s2
                satm = Atm_4shell(s)
                x4 = coords(1,satm)
                y4 = coords(2,satm)
                z4 = coords(3,satm)

                 call pack_coeffs(alphas, ixalpha, coeffs, ixpcoef,
     *                          ncfps, npfps, m, n,
     *                          r, s, alpha_pack, nalpha_pack,
     *                          pcoeff_pack, npcoeff_pack,
     *                          ccbeg, ccend, indx_cc,
     *                          ccbeg_pack, ccend_pack)

c---------------------------------------------------------------------------
c   Calling sequence for ERD version 2.
c---------------------------------------------------------------------------
       ncsum = ncfps(m) + ncfps(n) + ncfps(r) + ncfps(s)

C
C      Write(6,*) "@-entry to ERD_batch"
C      Write(6,"((a,4(1x,i4)))") "m, n, r, s:", m, n, r, s
C      Write(6,"(a,2(1x,I4))") "The Nalpha_pack and Npcoeff_pack:",
C     &                         nalpha_pack, npcoeff_pack
C      Write(6,"((a,4(1x,i4)))") "Ivanagmom:", Ivangmom(m),
C     &      Ivangmom(n), Ivangmom(r), Ivangmom(s)
C      Write(6,"((a,4(1x,i4)))") "Ncfps:", Ncfps(m),
C     &      Ncfps(n), Ncfps(r), Ncfps(s)
C      Write(6,"((a,4(1x,i4)))") "Npfps:", Npfps(m),
C     &      Npfps(n), Npfps(r), Npfps(s)
C      Write(6,*)
C      Write(6,"(4(1x,F15.7))") (Alpha_pack(i),i=1,Nalpha_pack)
C      Write(6,*)
C      Write(6,"(4(1x,F15.7))")(Pcoeff_pack(i),i=1,npcoeff_pack)
C      Write(6,*)
C      Write(6,"(4(1x,I4))")(CCbeg_pack(i),i=1,Nalpha_pack)
C      Write(6,*)
C      Write(6,"(4(1x,I4))")(CCend_pack(i),i=1,Nalpha_pack)
C      Write(6,*)
C      Write(6,"(a,3(1x,I10))") "The intmax,Zmax,ncsum:",intmax, zmax,
C     &                       ncsum
C      Write(6,*) "The coordiantes"
C      Write(6,*)
C      Write(6,"(3(1x,F10.5))") x1, y1, z1
C      Write(6,"(3(1x,F10.5))") x2, y2, z2
C      Write(6,"(3(1x,F10.5))") x3, y3, z3
C      Write(6,"(3(1x,F10.5))") x4, y4, z4
C
CSSS      if (m .le.2 .or. n .le. 2 .and. r .le. 2 .and.s.le.2) then
CSSS          Write(6,*)  "Entering ERD__GENER_spnspn_BATCH"

           call ERD__GENER_SPNSPN_CPLING_BATCH(intmax, zmax,
     *                nalpha_pack, npcoeff_pack, ncsum,
     *                ncfps(m),ncfps(n), ncfps(r), ncfps(s),
     *                npfps(m),npfps(n), npfps(r), npfps(s),
     *                ivangmom(m), ivangmom(n),
     *                ivangmom(r), ivangmom(s), x1,y1,z1,
     *                x2,y2,z2,x3,y3,z3,x4,y4,z4,
     *                xx,xy,xz,yx,yy,yz,zx,zy,zz,
     *                alpha_pack,
     *                pcoeff_pack, ccbeg_pack, ccend_pack,
     *                spherical, .true., iscr,
     *                nints,
     *                nfirst, scr)

           if (nints .gt. 0) then
              if (s .eq. 1) then
                    dd1 = 1
              else
                   dd1 = end_nfps(s-1) + 1
              endif
                  dd2 = end_nfps(s)

              if (XX .eq. 1) call add_4cent_integrals(out_spnspn_XX, 
     *                           a1,a2,b1,b2,c1,
     *                           c2,d1,d2,scr(nfirst),
     *                           aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,1.0d0,
     *                           0)  
              if (XY .eq. 1) call add_4cent_integrals(out_spnspn_XY, 
     *                           a1,a2,b1,b2,c1,
     *                           c2,d1,d2,scr(nfirst),
     *                           aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,1.0d0,
     *                           0)  
              if (XZ .eq. 1) call add_4cent_integrals(out_spnspn_XZ, 
     *                           a1,a2,b1,b2,c1,
     *                           c2,d1,d2,scr(nfirst),
     *                           aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,1.0d0,
     *                           0)  

              if (YX .eq. 1) call add_4cent_integrals(out_spnspn_YX, 
     *                           a1,a2,b1,b2,c1,
     *                           c2,d1,d2,scr(nfirst),
     *                           aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,1.0d0,
     *                           0)  

              if (YY .eq. 1) call add_4cent_integrals(out_spnspn_XX, 
     *                           a1,a2,b1,b2,c1,
     *                           c2,d1,d2,scr(nfirst),
     *                           aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,1.0d0,
     *                           0)  

              if (YZ .eq. 1) call add_4cent_integrals(out_spnspn_YZ, 
     *                           a1,a2,b1,b2,c1,
     *                           c2,d1,d2,scr(nfirst),
     *                           aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,1.0d0,
     *                           0)  

              if (ZX .eq. 1) call add_4cent_integrals(out_spnspn_ZX, 
     *                           a1,a2,b1,b2,c1,
     *                           c2,d1,d2,scr(nfirst),
     *                           aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,1.0d0,
     *                           0)  

              if (ZY .eq. 1) call add_4cent_integrals(out_spnspn_ZY, 
     *                           a1,a2,b1,b2,c1,
     *                           c2,d1,d2,scr(nfirst),
     *                           aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,1.0d0,
     *                           0)  

              if (ZZ .eq. 1) call add_4cent_integrals(out_spnspn_XX,
     *                           a1,a2,b1,b2,c1,
     *                           c2,d1,d2,scr(nfirst),
     *                           aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,1.0d0,
     *                           0)   
          endif 
CSSS          endif  !nints>0

      enddo   ! s
      enddo   ! r
      enddo   ! n
      enddo   ! m

      enddo ! ixyx
      enddo ! jxyz

      Write(6,*) "The <munu |(r1-r2)*r2/r12| lam sig> integral matrix"
      s22  = 0.0d0

      xx_s = 0.0d0
      xy_s = 0.0d0
      xz_s = 0.0d0
      yx_s = 0.0d0
      yy_s = 0.0d0
      yz_s = 0.0d0
      zx_s = 0.0d0
      zy_s = 0.0d0
      zz_s = 0.0d0

      xx_l = 0.0d0
      xy_l = 0.0d0
      xz_l = 0.0d0
      yx_l = 0.0d0
      yy_l = 0.0d0
      yz_l = 0.0d0
      zx_l = 0.0d0
      zy_l = 0.0d0
      zz_l = 0.0d0
 
      icount = 0

CSSS      do a = 1, Nbasis
CSSS      do b = 1, Nbasis
CSSS      do c = 1, Nbasis
CSSS      do d = 1, Nbasis
CSSS 
CSSS      ZZ  = 1
CSSS
CSSS      if (ZZ .EQ. 1) then
CSSS          Call rescale_ints(out_spnspn_ZZ, ACES_4DSH, Nbasis,
CSSS     &                      ERD_index,
CSSS     &                      ERD_scale)
CSSS
CSSS         if (d .le. c .and. b .le. a) then
CSSS          Write(6, "(4(1x,i3),1x,F15.10))")
Ca,b,c,d,out_spnspn_XX(a,b,c,d)
CSSS          s22 = s22 + ACES_4DSH(a,b,c,d)*ACES_4DSH(a,b,c,d)
CSSS          icount = icount + 1
CSSS         endif
CSSS      endif 

      do a = a1,a2
      do b = b1,b2
      do c = c1,c2
      do d = d1,d2

CSSS      if (d .le.c .and. b .le. a) then

CSSS      Write(*,*) "XX integrals"
CSSS      Write(6, "(5(1x,F15.10))") (Out_spnspn_XX(a,b,c,d), a=1, nbasis)
      xx_s = xx_s + Out_spnspn_XX(a,b,c,d) * Out_spnspn_XX(a,b,c,d) 
      xx_l = xx_l +  Out_spnspn_XX(a,b,c,d) 
C      write(*,*)

CSSS      Write(*,*) "XY integrals"
CSSS      Write(6, "(5(1x,F15.10))") (Out_spnspn_XY(a,b,c,d), a=1, nbasis)
      xy_s = xy_s + Out_spnspn_XY(a,b,c,d) * Out_spnspn_XY(a,b,c,d) 
      xy_l = xy_l + Out_spnspn_XY(a,b,c,d) 
C      write(*,*)

CSSS      Write(*,*) "XZ integrals"
CSSS      Write(6, "(5(1x,F15.10))") (Out_spnspn_XZ(a,b,c,d), a=1, nbasis)
      xz_s = xz_s + Out_spnspn_XZ(a,b,c,d) * Out_spnspn_XZ(a,b,c,d) 
      xz_l = xz_l + Out_spnspn_XZ(a,b,c,d) 
C      write(*,*)

CSSS      Write(*,*) "YX integrals"
CSSS      Write(6, "(5(1x,F15.10))") (Out_spnspn_YX(a,b,c,d), a=1, nbasis)
      yx_s = yx_s + Out_spnspn_YX(a,b,c,d) * Out_spnspn_YX(a,b,c,d) 
      yx_l = yx_l + Out_spnspn_YX(a,b,c,d) 
C      write(*,*)

CSSS      Write(*,*) "YY integrals"
CSSS      Write(6, "(5(1x,F15.10))") (Out_spnspn_YY(a,b,c,d), a=1, nbasis)
      yy_s = yy_s + Out_spnspn_YY(a,b,c,d) * Out_spnspn_YY(a,b,c,d) 
      yy_l = yy_l + Out_spnspn_YY(a,b,c,d) 
C      write(*,*)

CSSS      Write(*,*) "YZ integrals"
CSSS      Write(6, "(5(1x,F15.10))") (Out_spnspn_YZ(a,b,c,d), a=1, nbasis)
      yz_s = yz_s + Out_spnspn_YZ(a,b,c,d) * Out_spnspn_YZ(a,b,c,d) 
      yz_l = yz_l + Out_spnspn_YZ(a,b,c,d) 
C      write(*,*)

CSSS      Write(*,*) "ZX integrals"
CSSS      Write(6, "(5(1x,F15.10))") (Out_spnspn_ZX(a,b,c,d), a=1, nbasis)
      zx_s = zx_s + Out_spnspn_ZX(a,b,c,d) * Out_spnspn_ZX(a,b,c,d) 
      zx_l = zx_l + Out_spnspn_ZX(a,b,c,d) 
C      write(*,*)

CSSS      Write(*,*) "ZY integrals"
CSSS      Write(6, "(5(1x,F15.10))") (Out_spnspn_ZY(a,b,c,d), a=1, nbasis)
      zy_s = zy_s + Out_spnspn_ZY(a,b,c,d) * Out_spnspn_ZY(a,b,c,d) 
      zy_l = zy_l + Out_spnspn_ZY(a,b,c,d) 
C      write(*,*)

CSSS      Write(*,*) "ZZ integrals"
CSSS      Write(6, "(5(1x,F15.10))") (Out_spnspn_ZZ(a,b,c,d), a=1, nbasis)
      zz_s = zz_s + Out_spnspn_ZZ(a,b,c,d) * Out_spnspn_ZZ(a,b,c,d) 
      zz_l = zz_l + Out_spnspn_ZZ(a,b,c,d) 
CSSS      endif 

      enddo
      enddo
      enddo
      enddo

      Write(*,*)
      Write(6,"(a,a,F20.10,1x,F20.10)")"Check sum of the XX-spnspn",
     &                      "integrals = " ,xx_s,xx_l
      Write(6,"(a,a,F20.10,1x,F20.10)")"Check sum of the XY-spnspn",
     &                      "integrals = " ,xy_s,xy_l
      Write(6,"(a,a,F20.10,1x,F20.10)")"Check sum of the XZ-spnspn",
     &                      "integrals = " ,xz_s,xz_l
      Write(6,"(a,a,F20.10,1x,F20.10)")"Check sum of the YX-spnspn",
     &                      "integrals = " ,YX_s,YX_l
      Write(6,"(a,a,F20.10,1x,F20.10)")"Check sum of the YY-spnspn",
     &                      "integrals = " ,YY_s,YY_l
      Write(6,"(a,a,F20.10,1x,F20.10)")"Check sum of the YZ-spnspn",
     &                      "integrals = " ,YZ_s,YZ_l
      Write(6,"(a,a,F20.10,1x,F20.10)")"Check sum of the ZX-spnspn",
     &                      "integrals = " ,ZX_s,ZX_l
      Write(6,"(a,a,F20.10,1x,F20.10)")"Check sum of the ZY-spnspn",
     &                      "integrals = " ,ZY_s,ZY_l
      Write(6,"(a,a,F20.10,1x,F20.10)")"Check sum of the ZZ-spnspn",
     &                      "integrals = " ,ZZ_s,ZZ_l
      Write(6,*) 
      Write(6,"(a,a,F20.10,1x,F20.10)")"Check sum of scaled spnspn",
     &                      "integrals = " ,s22

C#endif 

      return
      end
