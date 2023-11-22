      subroutine compute_4cent_integrals_SOI
     &                              (a1,a2,b1,b2,
     *                               c1,c2,d1,d2,
     *                              scr,
     &                              maxblk,iscr,
     *                              coords,coeffs,
     &                              alphas, ccbeg, 
     &                              ccend, end_nfps, nshells, nbasis,
     &                              zmax, intmax, npfps, ncfps, 
     &                              ivangmom,
     &                              ixalpha, indx_cc, Atm_4shell,
     &                              ixpcoef,ispherical, soi_Jz, soi_Jx,
     *                              soi_Jy,nsend, 
     &                              Nalpha, Npcoef,ERD_index,
     &                              ERD_scale,ACES_SOI)
c---------------------------------------------------------------------------
c   "Work" routine for the integral worker task.
c---------------------------------------------------------------------------

      implicit none

      INTEGER     DER1X,DER2X,DER3X,DER4X
      INTEGER     DER1Y,DER2Y,DER3Y,DER4Y
      INTEGER     DER1Z,DER2Z,DER3Z,DER4Z

      integer a1, a2, b1, b2, c1, c2, d1, d2 
      integer aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2
      integer adim, bdim, cdim, ddim,icount
      integer m1, m2, n1, n2, r1, r2, s1, s2
      integer i, j, n, m, r, s
      integer a,b,c,d
      integer nshells, nbasis, zmax, intmax, ispherical

      integer num_to_do, nsend
      integer nints, maxblk
      integer nalpha_pack, npcoeff_pack
      integer ncsum, next, nfirst
      integer me, ierr, nalpha, npcoef, nprim_quad, ncoef_quad
      integer matm, natm, ratm, satm

      integer imin, zmin, iblk, zblk, iatm
      integer ivangmom(*), ixalpha(*), ixpcoef(*), npfps(*)
      integer ncfps(*), indx_cc(*), Atm_4shell(*)
      integer ERD_index(*)

      logical skip
      logical mn_symmetry
      logical rs_symmetry
      logical mn_rs_symmetry
      logical*8 l8true, l8spherical
      logical spherical

      double precision x1,y1,z1
      double precision x2,y2,z2
      double precision x3,y3,z3
      double precision x4,y4,z4
      double precision x,y,z
      double precision s11, s22

      double precision coords(3,*), coeffs(*), alphas(*)
      double precision soi_Jz(a1:a2,b1:b2,c1:c2,d1:d2),
     *                 soi_Jx(a1:a2,b1:b2,c1:c2,d1:d2),
     *                 soi_Jy(a1:a2,b1:b2,c1:c2,d1:d2),
     *                 ACES_SOI(a1:a2,b1:b2,c1:c2,d1:d2)

      double precision scr(*),ERD_scale(*)
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

C
      do d = 1, Nbasis
      do c = 1, Nbasis
      do b = 1, Nbasis
      do a = 1, Nbasis
         soi_Jz(a,b,c,d) = 0.d0
         soi_Jx(a,b,c,d) = 0.d0
         soi_Jy(a,b,c,d) = 0.d0
      enddo
      enddo
      enddo
      enddo
C
CSSS      Write(6,*)
CSSS      Write(6,*) "The atom label for each shell"
CSSS      Write(6,"(4(1x,I4))") (Atm_4shell(i), i=1, Nshells)
CSSS      Write(6,*)

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

CSSS         Write(6,*)
CSSS         Write(6,*) "Loop over shell quads." 
CSSS         Write(6,"(4(a,1x,I4,1x,I4,1x))") "m1 and m2: ", m1, m2,
CSSS     &                                    "n1 and N2: ", n1, n2, 
CSSS     &                                    "r1 and r2: ", r1, r2, 
CSSS     &                                    "s1 and s2: ", s1, s2

C   These are always undifferentiated

         DER3X = 0
         DER3Y = 0
         DER3Z = 0
         DER4X = 0
         DER4Y = 0
         DER4Z = 0

c constructing < mu nu | Jz | del gamma > 
c more precisely
c Jz = <mu^x nu | del^y gamma> - <mu^y nu | del^x gamma >
c Part A
c <X ..| Y ..>

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
c      Write(6,*) "@-entry to ERD_batch"
c      Write(6,"((a,4(1x,i4)))") "m, n, r, s:", m, n, r, s
c      Write(6,"(a,2(1x,I4))") "The Nalpha_pack and Npcoeff_pack:",
c     &                         nalpha_pack, npcoeff_pack
c      Write(6,"((a,4(1x,i4)))") "Ivanagmom:", Ivangmom(m), 
c     &      Ivangmom(n), Ivangmom(r), Ivangmom(s)
c      Write(6,"((a,4(1x,i4)))") "Ncfps:", Ncfps(m), 
c     &      Ncfps(n), Ncfps(r), Ncfps(s)
c      Write(6,"((a,4(1x,i4)))") "Npfps:", Npfps(m), 
c     &      Npfps(n), Npfps(r), Npfps(s)
c      Write(6,*)
c      Write(6,"(4(1x,F15.7))") (Alpha_pack(i),i=1,Nalpha_pack)
c      Write(6,*) 
c      Write(6,"(4(1x,F15.7))")(Pcoeff_pack(i),i=1,npcoeff_pack)
c      Write(6,*) 
c      Write(6,"(4(1x,I4))")(CCbeg_pack(i),i=1,Nalpha_pack)
c      Write(6,*) 
c      Write(6,"(4(1x,I4))")(CCend_pack(i),i=1,Nalpha_pack)
c      Write(6,*)
c      Write(6,"(a,3(1x,I10))") "The intmax,Zmax,ncsum:",intmax, zmax,
c     &                       ncsum
c      Write(6,*) "The coordiantes"
c      Write(6,*)
c      Write(6,"(3(1x,F10.5))") x1, y1, z1
c      Write(6,"(3(1x,F10.5))") x2, y2, z2
c      Write(6,"(3(1x,F10.5))") x3, y3, z3
c      Write(6,"(3(1x,F10.5))") x4, y4, z4

        DER1X = 1
        DER1Y = 0
        DER1Z = 0
        DER2X = 0
        DER2Y = 1
        DER2Z = 0

       call ERD__SOI_BATCH_XY(intmax, zmax,  
     *                nalpha_pack, npcoeff_pack, ncsum, 
     *                ncfps(m),ncfps(n), ncfps(r), ncfps(s),
     *                npfps(m),npfps(n), npfps(r), npfps(s),
     *                ivangmom(m), ivangmom(n), 
     *                ivangmom(r), ivangmom(s),
     *                x1,y1,z1,
     *                x2,y2,z2,x3,y3,z3,x4,y4,z4, 
     &                DER1X,DER1Y,DER1Z,
     +                DER2X,DER2Y,DER2Z,
     +                DER3X,DER3Y,DER3Z,
     +                DER4X,DER4Y,DER4Z,
     *                alpha_pack,
     *                pcoeff_pack, ccbeg_pack, ccend_pack, 
     *                spherical, .true., iscr,
     *                nints, 
     *                nfirst, scr)    

c---------------------------------------------------------------------------
c   Move the integrals into the output block.  
c---------------------------------------------------------------------------

            if (nints .gt. 0) then

c       Write(6,*)"The <munu | lam sig> integral matrix for shell block"
c       Write(6,"(5(1x,F15.10))") (scr(nfirst-1+i),i=1,nints)

               
               if (s .eq. 1) then
                  dd1 = 1
               else
                  dd1 = end_nfps(s-1) + 1
               endif
               dd2 = end_nfps(s)

               call add_integrals(soi_Jz, a1,a2,b1,b2,c1,c2,d1,d2,
     *                             scr(nfirst), 
     *                            aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,1.0d0)
            endif

c  <Y ...| X..>

        DER1X = 0
        DER1Y = 1
        DER1Z = 0
        DER2X = 1
        DER2Y = 0
        DER2Z = 0

         call ERD__SOI_BATCH_YX(intmax, zmax, 
     *                nalpha_pack, npcoeff_pack, ncsum,
     *                ncfps(m),ncfps(n), ncfps(r), ncfps(s),
     *                npfps(m),npfps(n), npfps(r), npfps(s),
     *                ivangmom(m), ivangmom(n),
     *                ivangmom(r), ivangmom(s), x1,y1,z1,
     *                x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                DER1X,DER1Y,DER1Z,
     +                DER2X,DER2Y,DER2Z,
     +                DER3X,DER3Y,DER3Z,
     +                DER4X,DER4Y,DER4Z,
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
      call add_integrals(soi_Jz, a1,a2,b1,b2,c1,c2,d1,d2,
     *                            scr(nfirst),
     *                            aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,-1.d0)
      endif

         enddo   ! s
         enddo   ! r

         enddo   ! n
         enddo   ! m
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The <munu |SOI_JZ| lamsig> integral matrix"
      s11=0.0d0
      s22=0.0d0
      icount = 0
      Call rescale_ints(SOI_Jz, ACES_SOI, Nbasis, ERD_index, 
     &                  ERD_scale)
      do a = 1, Nbasis
CSSS      Write(6,*)
      do b = 1, Nbasis
CSSS      Write(6,*)
      do c = 1, Nbasis
      do d = 1, Nbasis 
       if (Dabs(SOI_Jz(a,b,c,d)) .gt. 10D-12) then
       if (a .le. b .and. c .le. d) then
       Write(6, "(4(1x,i3),1x,F15.10))") a,b,c,d,SOI_Jz(a,b,c,d)
       s11 = s11 + SOI_jz(a,b,c,d)*SOI_jz(a,b,c,d)
       s22 = s22 + ACES_SOI(a,b,c,d)*ACES_SOI(a,b,c,d)
       icount = icount + 1
       endif
       endif
CSSS      Write(6, "(5(1x,F15.10))") (SOI_JZ(a,b,c,d), a=1, nbasis)
      enddo
      enddo
      enddo
      enddo
      Write(*,*)
      Write(6,"(a,i8,a,F20.10)")"L_z check sum of",icount," integrals= "
     &                           ,s11
     
      Write(6,"(a,i8,a,F20.10)")"L_z check sum of",icount,
     &                             " scaled integrals= "
     &                             ,s22
#endif
      

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c Jx part
c (Y...|Z...) - (Z..|Y..)

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

        DER1X = 0
        DER1Y = 1
        DER1Z = 0
        DER2X = 0
        DER2Y = 0
        DER2Z = 1

       call ERD__SOI_BATCH_XY(intmax, zmax,
     *                nalpha_pack, npcoeff_pack, ncsum,
     *                ncfps(m),ncfps(n), ncfps(r), ncfps(s),
     *                npfps(m),npfps(n), npfps(r), npfps(s),
     *                ivangmom(m), ivangmom(n),
     *                ivangmom(r), ivangmom(s),
     *                x1,y1,z1,
     *                x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                DER1X,DER1Y,DER1Z,
     +                DER2X,DER2Y,DER2Z,
     +                DER3X,DER3Y,DER3Z,
     +                DER4X,DER4Y,DER4Z,
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
       call add_integrals(soi_Jx, a1,a2,b1,b2,c1,c2,d1,d2,
     *                             scr(nfirst),
     *                            aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,1.0d0)
       endif
c  <Z ...| Y..>

        DER1X = 0
        DER1Y = 0
        DER1Z = 1
        DER2X = 0
        DER2Y = 1
        DER2Z = 0

         call ERD__SOI_BATCH_YX(intmax, zmax,
     *                nalpha_pack, npcoeff_pack, ncsum,
     *                ncfps(m),ncfps(n), ncfps(r), ncfps(s),
     *                npfps(m),npfps(n), npfps(r), npfps(s),
     *                ivangmom(m), ivangmom(n),
     *                ivangmom(r), ivangmom(s), x1,y1,z1,
     *                x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                DER1X,DER1Y,DER1Z,
     +                DER2X,DER2Y,DER2Z,
     +                DER3X,DER3Y,DER3Z,
     +                DER4X,DER4Y,DER4Z,
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
      call add_integrals(soi_Jx, a1,a2,b1,b2,c1,c2,d1,d2,
     *                            scr(nfirst),
     *                            aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,-1.d0)
      endif
         enddo   ! s
         enddo   ! r

         enddo   ! n
         enddo   ! m
#ifdef _DEBUG_LVL0
      Write(6,*) "The <munu |SOI_JX| lamsig> integral matrix"
      Call rescale_ints(SOI_JX, ACES_SOI, Nbasis, ERD_index, 
     &                  ERD_scale)
      s11 = 0.0d0
      s22 = 0.0d0
      icount = 0
      do a = 1, Nbasis
CSSS      Write(6,*) 
      do b = 1, Nbasis
CSSS      Write(6,*) 
      do c = 1, Nbasis
      do d = 1, Nbasis
       if (Dabs(SOI_Jx(a,b,c,d)) .gt. 10D-12) then
       if (a .le. b .and. c .le. d) then
       Write(6, "(4(1x,i3),1x,F15.10))") a,b,c,d,SOI_JX(a,b,c,d)
       s11 = s11 + SOI_JX(a,b,c,d)*SOI_JX(a,b,c,d)
       s22 = s22 + ACES_SOI(a,b,c,d)*ACES_SOI(a,b,c,d)
       icount = icount + 1
      endif
      endif 
      enddo 
      enddo
      enddo
      enddo
      Write(6,"(a,i8,a,F20.10)")"L_x check sum of",icount," integrals= "
     &                            ,s11
      Write(6,"(a,i8,a,F20.10)")"L_x check sum of",icount,
     &                             " scaled integrals= "
     &                             ,s22
#endif
c+++++++++++++Jx .......Jx .....done++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c Jy part
c (Z...|X...) - (X..|Z..)
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
       ncsum = ncfps(m) + ncfps(n) + ncfps(r) + ncfps(s)

        DER1X = 0
        DER1Y = 0
        DER1Z = 1
        DER2X = 1
        DER2Y = 0
        DER2Z = 0

       call ERD__SOI_BATCH_XY(intmax, zmax,
     *                nalpha_pack, npcoeff_pack, ncsum,
     *                ncfps(m),ncfps(n), ncfps(r), ncfps(s),
     *                npfps(m),npfps(n), npfps(r), npfps(s),
     *                ivangmom(m), ivangmom(n),
     *                ivangmom(r), ivangmom(s),
     *                x1,y1,z1,
     *                x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                DER1X,DER1Y,DER1Z,
     +                DER2X,DER2Y,DER2Z,
     +                DER3X,DER3Y,DER3Z,
     +                DER4X,DER4Y,DER4Z,
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
       call add_integrals(soi_Jy, a1,a2,b1,b2,c1,c2,d1,d2,
     *                             scr(nfirst),
     *                            aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,1.0d0)
       endif
c  <Z ...| Y..>
C      write(*,*) '<Z....|Y...> part'
        DER1X = 1
        DER1Y = 0
        DER1Z = 0
        DER2X = 0
        DER2Y = 0
        DER2Z = 1

         call ERD__SOI_BATCH_YX(intmax, zmax,
     *                nalpha_pack, npcoeff_pack, ncsum,
     *                ncfps(m),ncfps(n), ncfps(r), ncfps(s),
     *                npfps(m),npfps(n), npfps(r), npfps(s),
     *                ivangmom(m), ivangmom(n),
     *                ivangmom(r), ivangmom(s), x1,y1,z1,
     *                x2,y2,z2,x3,y3,z3,x4,y4,z4,
     &                DER1X,DER1Y,DER1Z,
     +                DER2X,DER2Y,DER2Z,
     +                DER3X,DER3Y,DER3Z,
     +                DER4X,DER4Y,DER4Z,
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
      call add_integrals(soi_Jy, a1,a2,b1,b2,c1,c2,d1,d2,
     *                            scr(nfirst),
     *                            aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,-1.d0)
      endif

      enddo   ! s
      enddo   ! r

      enddo   ! n
      enddo   ! m
        
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The <munu |SOI_JY| lamsig> integral matrix"
      Call rescale_ints(SOI_JY, ACES_SOI, Nbasis, ERD_index, 
     &                  ERD_scale)
      s11 = 0.0d0
      s22 = 0.0d0
      icount = 0
      do a = 1, Nbasis
CSSS      Write(6,*) 
      do b = 1, Nbasis
CSSS      Write(6,*) 
      do c = 1, Nbasis
      do d = 1, Nbasis
       if (Dabs(SOI_Jy(a,b,c,d)) .gt. 10D-12) then
       if (a .le. b .and. c .le. d) then
          Write(6, "(4(1x,i3),1x,F15.10))") a,b,c,d,SOI_JY(a,b,c,d)      
         s11 = s11 + SOI_JY(a,b,c,d)*SOI_JY(a,b,c,d)
         s22 = s22 + ACES_SOI(a,b,c,d)*ACES_SOI(a,b,c,d)
         icount = icount + 1
      endif
      endif
      enddo 
      enddo
      enddo
      enddo
      Write(6,"(a,i8,a,F20.10)")"L_y check sum of",icount," integrals= "
     &                        ,s11
      Write(6,"(a,i8,a,F20.10)")"L_x check sum of",icount,
     &                             " scaled integrals= ",
     &                             s22
#endif

      return
      end

