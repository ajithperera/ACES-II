
      subroutine compute_oed_integrals
     &                              (a1,a2,b1,b2,scr,
     &                              maxblk,iscr,coords,coeffs,
     &                              coeffs_uc, alphas, ccbeg, 
     &                              ccend, end_nfps, nshells, nbasis,
     &                              zmax, intmax, npfps, ncfps, 
     &                              ivangmom,
     &                              ixalpha, indx_cc, Atm_4shell,
     &                              ixpcoef,ispherical,contract,
     &                              dkh_order,out_kin,
     &                              out_ovl,out_nai,out_pvp,
     &                              coreham,acesord, 
     &                              erd_index,erd_scale,
     &                              nsend, 
     &                              Nalpha,Nprims,Ncnfns,Npcoef,
     &                              npcoef_uc,
     &                              ncenters,
     &                              charge,
     &                              max_centers,zcore,maxcor)
c---------------------------------------------------------------------------
c   "Work" routine for the integral worker task.
c---------------------------------------------------------------------------
      implicit none
c
      integer a1, a2, b1, b2,max_centers 
      integer aa1,aa2,bb1,bb2, iorder 
      integer adim, bdim, cdim, ddim  
      integer m1, m2, n1, n2
      integer i, j, n, m, r, s,l,ml
      integer a,b,c,d
      integer nshells, nbasis, zmax, intmax, ispherical
      integer max_atoms, nprims, maxcor, n_confns
      integer dkh_order 

      integer num_to_do, nsend
      integer nints, maxblk
      integer nalpha_pack, npcoeff_pack,npcoeff
      integer ncsum, next, nfirst
      integer me, ierr, nalpha, npcoef, nprim_quad, ncoef_quad
      integer matm, natm, ratm, satm, npcoef_uc, ncnfns 
      integer erd_index(*)

      integer imin, zmin, iblk, zblk, iatm
      integer ivangmom(*), ixalpha(*), ixpcoef(*), npfps(*)
      integer ncfps(*), indx_cc(*), Atm_4shell(*)
      double precision xn(max_centers), yn(max_centers),
     *                 zn(max_centers)
      double precision erd_scale(*) 
      logical spherical, contract 

      double precision x1,y1,z1
      double precision x2,y2,z2

      integer ncenters
      integer iatom
      integer Ea1,Ea2,Eb1,Eb2

      double precision coords(3,*), coeffs(*), alphas(*),
     *                 charge(*), coeffs_uc(*)
      double precision out_kin(a1:a2,b1:b2)
      double precision out_ovl(a1:a2,b1:b2)
      double precision out_nai(a1:a2,b1:b2)
      double precision out_pvp(a1:a2,b1:b2)
      double precision coreham(a1:a2,b1:b2)
      double precision acesord(a1:a2,b1:b2)
      double precision scr(*)   
      double precision zcore(maxcor)
      integer iscr(*)

      integer ccbeg(*), ccend(*), end_nfps(*)

      integer max_dim_coeff
      parameter (max_dim_coeff = 5000)
      integer ccbeg_pack(max_dim_coeff), ccend_pack(max_dim_coeff)
      double precision alpha_pack(max_dim_coeff), 
     *                 pcoeff_pack(max_dim_coeff)
      save me,alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack

      adim = Nprims 
      bdim = Nprims 

      spherical = (ispherical .eq. 1)
          
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
      Do Iatm =1, Ncenters
         Write(6, "(6(1x,F12.8))") (Coords(i,Iatm), i=1,3)
      Enddo
      Write(6,"(a,l)") "Spherical?", Spherical
      write(6,*) "The contraction coefs." 
      Write(6, "(6(1x,F12.8))") (coeffs(i),i=1,npcoef)
      write(6,*) "The primitive functions." 
      Write(6,"(4(1x,I4))") (npfps(i), i=1, Nshells)
      write(6,*) "The contracted functions." 
      Write(6,"(4(1x,I4))") (ncfps(i), i=1, Nshells)
#endif

c-----------------------------------------------------------------------
c   Find the shell blocks for which we shall loop through.
c----------------------------------------------------------------------
      call lookup_shell(end_nfps, nshells, a1, m1)
      call lookup_shell(end_nfps, nshells, a2, m2)
      call lookup_shell(end_nfps, nshells, b1, n1)
      call lookup_shell(end_nfps, nshells, b2, n2)

      do b = 1, Nprims 
      do a = 1, Nprims 
         out_kin(a,b) = 0.d0
         out_ovl(a,b) = 0.d0
         out_nai(a,b) = 0.d0
         out_pvp(a,b) = 0.d0
      end do
      end do
        
      do m = 1, nshells
         iatom = Atm_4shell(m)
         xn(iatom) = coords(1,iatom)
         yn(iatom) = coords(2,iatom)
         zn(iatom) = coords(3,iatom)
      enddo
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

            call pack_coeffs_oed(alphas, ixalpha, coeffs_uc, ixpcoef, 
     *                           npfps, npfps, m, n, 
     *                           alpha_pack, nalpha, 
     *                           pcoeff_pack, npcoeff, 
     *                           ccbeg, ccend, indx_cc,
     *                           ccbeg_pack, ccend_pack,max_dim_coeff)

           ncsum = ncfps(m) + ncfps(n) 

CSSS      Write(6,*) "@-entry to OED_batch"
CSSS      Write(6,"((a,2(1x,i4)))") "m, n :", m, n

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
#endif 
C Kinetic energy

           call oed__gener_kin_batch(intmax,zmax,
     *                 nalpha,npcoeff,
     *                 ncsum, npfps(m), npfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
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
              call add_integrals2(out_kin, a1,a2,b1,b2,
     *                            scr(nfirst),Ea1,Ea2,Eb1,Eb2,1.0d0)
           endif 
C
C Overlap 
C 
           call oed__gener_ovl_batch(intmax,zmax,
     *                 nalpha,npcoeff,
     *                 ncsum, npfps(m), npfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
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
               call add_integrals2(out_ovl, a1,a2,b1,b2,
     *                             scr(nfirst),Ea1,Ea2,Eb1,Eb2,1.0d0)
            endif 
C
C Nuclear attraction
C
           call oed__gener_nai_batch(intmax,zmax,
     *                 nalpha,npcoeff,
     *                 ncsum, npfps(m), npfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 ncenters,xn,yn,zn,charge,
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
               call add_integrals2(out_nai, a1,a2,b1,b2,
     *                             scr(nfirst),Ea1,Ea2,Eb1,Eb2,1.0d0)
          endif 
          
C PVP batch
C
           call oed__gener_pvp_batch
     *                (intmax, zmax, nalpha, npcoeff,
     *                 ncsum, npfps(m), npfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 ncenters, xn,yn,zn, charge,
     *                 1,0,0,
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

               call move_integrals2(out_pvp, a1,a2,b1,b2,
     *                           scr(nfirst),Ea1,Ea2,Eb1,Eb2)
           Endif 

           call oed__gener_pvp_batch
     *                (intmax, zmax, nalpha, npcoeff,
     *                 ncsum, npfps(m), npfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 ncenters, xn,yn,zn, charge,
     *                 0,1,0,
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


               call add_integrals2(out_pvp, a1,a2,b1,b2,
     *                           scr(nfirst),Ea1,Ea2,Eb1,Eb2,1.0d0)
           endif 

           call oed__gener_pvp_batch
     *                (intmax, zmax, nalpha, npcoeff,
     *                 ncsum, npfps(m), npfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 ncenters, xn,yn,zn, charge,
     *                 0,0,1,
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

               call add_integrals2(out_pvp, a1,a2,b1,b2,
     *                           scr(nfirst),Ea1,Ea2,Eb1,Eb2,1.0d0)
            endif 

         end do !n
         end do !m

C      enddo
C
#ifdef _DEBUG_LVL0
      Call rescale_2ints(out_pvp,acesord,Nbasis,Erd_index,ERD_scale)
      Call dcopy(Nbasis*Nbasis,acesord,1,out_pvp,1) 
      Call rescale_2ints(out_nai,acesord,Nbasis,Erd_index,ERD_scale)
      Call dcopy(Nbasis*Nbasis,acesord,1,out_nai,1) 
      Call rescale_2ints(out_ovl,acesord,Nbasis,Erd_index,ERD_scale)
      Call dcopy(Nbasis*Nbasis,acesord,1,out_ovl,1) 
      Call rescale_2ints(out_kin,acesord,Nbasis,Erd_index,ERD_scale)
      Call dcopy(Nbasis*Nbasis,acesord,1,out_kin,1) 
      Write(6,*) "The <mu| P^2 |nu > integral matrix"
      call output(out_kin, 1, nprims, 1, nprims, nprims, nprims, 1)
      call checksum("KIN",out_kin,Nprims*Nprims)
      Write(6,*) "The <mu| S |nu >  matrix"
      call output(out_ovl, 1, nprims, 1, nprims, nprims, nprims, 1)
      call checksum("S",out_ovl,Nprims*Nprims)
      Write(6,*) "The <mu| PVP |nu >  matrix"
      call output(out_pvp, 1, nprims, 1, nprims, nprims, nprims, 1)
      call checksum("pvp",out_pvp,Nprims*Nprims)
      Write(6,*) "The <mu| NAI |nu >  matrix"
      call output(out_nai, 1, nprims, 1, nprims, nprims, nprims, 1)
      call checksum("Nai",out_nai,Nprims*Nprims)
C
#endif
      Iorder = dkh_order 
      Write(6,*)
      Write(6,"(a,1x,i1,1x,a)") " The DKH Hamilonian of order", Iorder,
     +                      "is computed" 
      Call dkh__init(Iorder,Nprims,out_kin,out_nai,out_pvp,out_ovl,
     +               zcore,maxcor,coreham)
#ifdef _DEBUG_LVLM
      Call rescale_2ints(coreham,acesord,NPrims,Erd_index,ERD_scale)
      call checksum("A2_ORD_UNC_DK",acesord,Nprims*Nprims)
#endif 
      If (Contract) then
      
         Call Dzero(out_pvp,nprims*nprims)
         Call dkh__contract(zcore,maxcor,coreham,out_pvp,coeffs,nshells,
     +                   nprims,ncnfns,npfps,npcoef,npcoef_uc,ncfps,
     +                   ivangmom,spherical)
         do i=1,Nshells
            l = ivangmom(i)
            if (spherical) then
               ml=2*l+1
               Ncfps(i) = Ncfps(i)*ml
            else
               ml = (l+1)*(l+2)/2
              Ncfps(i) = Ncfps(i)*ml
            endif
        enddo

        Call Gen_oed_to_vmol(Ncfps, Ivangmom, Nshells, ispherical,
     &                       Erd_index,ERD_Scale)
      Else
         Call Dzero(out_pvp,nprims*nprims)
         Call Dcopy(nprims*nprims,coreham,1,out_pvp,1)
      Endif 

#ifdef _DEBUG
C Contration routine debugging. Lets contract a known integral like overlap 

      Write(6,*) "Testing the contraction routines"
      Call Dzero(coreham,nprims*nprims)
      Call dkh__contract(zcore,maxcor,out_ovl,coreham,coeffs,nshells,
     +                   nprims,ncnfns,npfps,npcoef,npcoef_uc,ncfps,
     +                   ivangmom,spherical)
      n_confns = 0
      do i=1,Nshells
      l = ivangmom(i) 
      if (spherical) then
         ml=2*l+1
         Ncfps(i) = Ncfps(i)*ml
      else
         ml = (l+1)*(l+2)/2 
         Ncfps(i) = Ncfps(i)*ml
      endif 
      enddo 

      Call Gen_oed_to_vmol(Ncfps, Ivangmom, Nshells, ispherical,
     &                 Erd_index,ERD_Scale)
      Call rescale_2ints(coreham,out_ovl,Ncnfns,Erd_index,ERD_scale)
      Write(6,*) "The <mu| S |nu >  matrix"
      call output(out_ovl, 1, Ncnfns, 1, Ncnfns, Ncnfns, Ncnfns, 1)
      write(6,"(a,i4)") "The number of contracted funtions",ncnfns
      Call checksum("contracted Overlap ints",out_ovl,Ncnfns*Ncnfns)
#endif

      return
      end
