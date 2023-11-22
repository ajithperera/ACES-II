      subroutine compute_4cent_integrals
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
     &                              ixpcoef,ispherical, 
     *                              Repulsion, 
     &                              Nalpha, Npcoef)
c---------------------------------------------------------------------------
c   "Work" routine for the integral worker task.
c---------------------------------------------------------------------------

      implicit none

      integer a1, a2, b1, b2, c1, c2, d1, d2 
      integer aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2
      integer adim, bdim, cdim, ddim  
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
      double precision eri_s,eri_l

      double precision coords(3,*), coeffs(*), alphas(*),
     *                 Repulsion(a1:a2,b1:b2,c1:c2,d1:d2)

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

      do m = 1, Nbasis 
      do n = 1, Nbasis 
      do r = 1, Nbasis 
      do s = 1, Nbasis 
         Repulsion(m,n,r,s) = 0.0d0
      enddo 
      enddo 
      enddo 
      enddo   

CSSS         Write(6,*)
CSSS         Write(6,*) "Loop over shell quads." 
CSSS         Write(6,"(4(a,1x,I4,1x,I4,1x))") "m1 and m2: ", m1, m2,
CSSS     &                                    "n1 and N2: ", n1, n2, 
CSSS     &                                    "r1 and r2: ", r1, r2, 
CSSS     &                                    "s1 and s2: ", s1, s2

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

 
      Write(6,*) "@-entry to ERD_batch"
      Write(6,"((a,4(1x,i4)))") "m, n, r, s:", m, n, r, s
      Write(6,"(a,2(1x,I4))") "The Nalpha_pack and Npcoeff_pack:",
     &                         nalpha_pack, npcoeff_pack
      Write(6,"((a,4(1x,i4)))") "Ivanagmom:", Ivangmom(m), 
     &      Ivangmom(n), Ivangmom(r), Ivangmom(s)
      Write(6,"((a,4(1x,i4)))") "Ncfps:", Ncfps(m), 
     &      Ncfps(n), Ncfps(r), Ncfps(s)
      Write(6,"((a,4(1x,i4)))") "Npfps:", Npfps(m), 
     &      Npfps(n), Npfps(r), Npfps(s)
      Write(6,*)
      Write(6,"(4(1x,F15.7))") (Alpha_pack(i),i=1,Nalpha_pack)
      Write(6,*) 
      Write(6,"(4(1x,F15.7))")(Pcoeff_pack(i),i=1,npcoeff_pack)
      Write(6,*) 
      Write(6,"(4(1x,I4))")(CCbeg_pack(i),i=1,Nalpha_pack)
      Write(6,*) 
      Write(6,"(4(1x,I4))")(CCend_pack(i),i=1,Nalpha_pack)
      Write(6,*)
      Write(6,"(a,3(1x,I10))") "The intmax,Zmax,ncsum:",intmax, zmax,
     &                       ncsum
      Write(6,*) "The coordiantes"
      Write(6,*)
      Write(6,"(3(1x,F10.5))") x1, y1, z1
      Write(6,"(3(1x,F10.5))") x2, y2, z2
      Write(6,"(3(1x,F10.5))") x3, y3, z3
      Write(6,"(3(1x,F10.5))") x4, y4, z4

       call ERD__GENER_ERI_BATCH(intmax, zmax,  
     *                nalpha_pack, npcoeff_pack, ncsum, 
     *                ncfps(m),ncfps(n), ncfps(r), ncfps(s),
     *                npfps(m),npfps(n), npfps(r), npfps(s),
     *                ivangmom(m), ivangmom(n), 
     *                ivangmom(r), ivangmom(s),
     *                x1,y1,z1,
     *                x2,y2,z2,x3,y3,z3,x4,y4,z4, 
     *                alpha_pack,
     *                pcoeff_pack, ccbeg_pack, ccend_pack, 
     *                spherical, .true., iscr,
     *                nints, 
     *                nfirst, scr)    

c---------------------------------------------------------------------------
c   Move the integrals into the output block.  
c---------------------------------------------------------------------------

            if (nints .gt. 0) then
               
               if (s .eq. 1) then
                  dd1 = 1
               else
                  dd1 = end_nfps(s-1) + 1
               endif
               dd2 = end_nfps(s)
               call add_integrals(Repulsion, a1,a2,b1,b2,c1,c2,d1,d2,
     *                             scr(nfirst), 
     *                            aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,1.0d0)
            endif

      enddo   ! s
      enddo   ! r

      enddo   ! n
      enddo   ! m
        
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The <munu |1/r12| lam sig> integral matrix"
      eri_s = 0.0
      eri_l = 0.0
      do a = a1,a2
CSSS      Write(6,*) 
      do b = b1,b2
CSSS      Write(6,*) 
      do c = c1,c2
CSSS      do d = d1,d2
      Write(6, "(5(1x,F15.10))") (Repulsion(d,c,b,a), d=1, nbasis)
CSSS         eri_s = eri_s + Repulsion(a,b,c,d) * Repulsion(a,b,c,d)
CSSS         eri_l = eri_l + Repulsion(a,b,c,d) + Repulsion(a,b,c,d)
      enddo
      enddo
      enddo
CSSS      enddo
CSSS      Write(*,*)
CSSS      Write(6,"(a,a,F20.10,1x,F20.10)")"Check sum of the repulsion",
CSSS     &                      " integrals = " ,eri_s,eri_l
#endif

      return
      end

