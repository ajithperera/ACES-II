      subroutine compute_delta_integrals
     &                              (a1,a2,b1,b2,scr,
     &                              maxblk,iscr,coords,coeffs,
     &                              alphas, ccbeg, 
     &                              ccend, end_nfps, nshells, nbasis,
     &                              zmax, intmax, npfps, ncfps, 
     &                              ivangmom,
     &                              ixalpha, indx_cc, Atm_4shell,
     &                              ixpcoef,ispherical, delta_int,
     *                              nsend, 
     &                              Nalpha, Npcoef,nuclei, max_centers)
c---------------------------------------------------------------------------
c   "Work" routine for the integral worker task.
c    Writtent by Prakash Nov. 2012. 
c    Program calculate < mu(r-Xa) | delta(r-Rn) | nu(r-Xb) >
c                      = mu(Rn-Xa)*nu(Rn-Xb)
c   for the centre N.
c---------------------------------------------------------------------------
      implicit none
     
      integer max_centers 

      integer a1, a2, b1, b2 
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

      integer imin, zmin, iblk, zblk
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

cPV we need this for NAI Lo, L-1, L+1
       integer nuclei 
       double precision xnc(max_centers), ync(max_centers),
     *                  znc(max_centers), xn, yn, zn
      integer cent, iatom 
      integer lmnvala(3,36),lmnvalb(3,36),i1
      integer Ea1,Ea2,Eb1,Eb2
c end PV
      double precision coords(3,*), coeffs(*), alphas(*)
      double precision delta_int(a1:a2,b1:b2,1:Nuclei)
      double precision scr(*)   
      double precision erfcoef
      integer iscr(*)
      integer ix_old,jx_old
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

      do N=1,Nuclei 
      do b = 1, Nbasis
      do a = 1, Nbasis
         delta_int(a,b,N) = 0.d0
      end do
      end do
      end do 

c-----------------------------------------------------------------------
c   Find the shell blocks for which we shall loop through.
c-----------------------------------------------------------------------

      call lookup_shell(end_nfps, nshells, a1, m1)
      call lookup_shell(end_nfps, nshells, a2, m2)
      call lookup_shell(end_nfps, nshells, b1, n1)
      call lookup_shell(end_nfps, nshells, b2, n2)

cPV  obtain the nulear coordinate

      do m=1,nshells
         iatom=Atm_4shell(m)
         xnc(iatom) = coords(1,iatom)
         ync(iatom) = coords(2,iatom)
         znc(iatom) = coords(3,iatom)
#ifdef _DEBUG_LVL0
         Write(6,*)
         Write(6,"(a)") "The center coordinates"
         write(6,"(1x,I3,3F12.6)") iatom, xnc(iatom),ync(iatom),
     &                             znc(iatom)
#endif
      end do
cPV end prakash
   
       do cent = 1, Nuclei
   
         Xn = Xnc(cent)
         yn = Ync(cent)
         zn = Znc(cent)


      do m = m1, m2
CSSS      do m = 1, 2

        matm = Atm_4shell(m)

         x1 = coords(1,matm)
         y1 = coords(2,matm)
         z1 = coords(3,matm)

         do n = n1, n2
CSSS         do n = 1, 2

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
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "@-entry to ERD_batch"
      Write(6,"((a,4(1x,i4)))") "m, n :", m, n
      Write(6,"(a,2(1x,I4))") "The Nalpha_pack and Npcoeff_pack:",
     &                          nalpha, npcoeff
      Write(6,"((a,2(1x,i4)))") "Ivanagmom:", Ivangmom(m),
     &      Ivangmom(n)
      Write(6,"((a,2(1x,i4)))") "Ncfps:", Ncfps(m),
     &      Ncfps(n)
      Write(6,"((a,4(1x,i4)))") "Npfps:", Npfps(m),
     &      Npfps(n)
      Write(6,"(4(1x,F15.7))") (Alpha_pack(i),i=1,Nalpha)
      Write(6,*)
      Write(6,"(4(1x,F15.7))")(Pcoeff_pack(i),i=1,npcoeff)
      Write(6,*)
      Write(6,"(6(1x,I4))")(CCbeg_pack(i),i=1,Nalpha)
      Write(6,*)
      Write(6,"(6(1x,I4))")(CCend_pack(i),i=1,Nalpha)
      Write(6,*)
      Write(6,"(a,3(1x,I10))") "The intmax,Zmax,ncsum:",intmax, zmax,
     &                       ncsum
      Write(6,*) "The shell and nuclear coordiantes"
      Write(6,*)
      Write(6,"(3(1x,F12.6))") x1, y1, z1
      Write(6,"(3(1x,F12.6))") x2, y2, z2
      Write(6,"(3(1x,F12.6))") xn, yn, zn
#endif

      call oed__gener_delta_batch(intmax, zmax, nalpha, npcoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 Xn,Yn,Zn,alpha_pack, 
     *                 pcoeff_pack, ccbeg_pack, ccend_pack,spherical,
     *                 iscr,nints, nfirst, 
     *                 scr)

       if(nints .gt. 0) then

#ifdef _DEBUG_LVL0
       Write(6,*)
       Write(6,*)"The <mu|delta|nu> integral matrix for shell block"
       Write(6,"(5(1x,F15.10))") (scr(nfirst-1+i),i=1,nints)
       Write(6,*) 
#endif

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

         call move_integrals2_delta(delta_int, 1,Nbasis,1,Nbasis,
     *                              scr(nfirst),Ea1,Ea2,Eb1,Eb2,
     *                              cent)
      end if

       end do !n
       end do !m
C
#ifdef _DEBUG_LVL0
       Write(6,*) "The <mu|delta|nu > integral matrix for Center :", 
     &              Cent
       call output(delta_int(1,1,cent), 1, Nbasis, 1, Nbasis, Nbasis, 
     &             Nbasis, 1)
#endif

       end do !centers 
C
      return
      end

