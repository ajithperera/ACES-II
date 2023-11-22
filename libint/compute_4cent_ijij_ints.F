      Subroutine Compute_4cent_ijij_ints(a1,a2,b1,b2,Scr,Iscr,Coords,
     &                                   Coeffs,pcoeff_pack,Alpha,
     &                                   Alpha_pack,Ccbeg,Ccbeg_pack,
     &                                   Ccend,Ccend_pack,End_nfps,
     &                                   Nshells,Nbasis,Max_cfps,
     &                                   Max_pfps,Zmax,Intmax,Npfps,
     &                                   Ncfps,Ivangmom,Ixalpha,Indx_cc,
     &                                   Atm_4shell,Ixpcoef,
     &                                   Spherical,Nalpha,Npcoef,
     &                                   Erd_index,Erd_scale,Ncenters,
     &                                   Nprims,Repulsion)

c---------------------------------------------------------------------------
c   "Work" routine for the integral worker task.
c---------------------------------------------------------------------------

      implicit none

      integer a1,a2,b1,b2
      integer m1,m2,n1,n2
      integer i,j,n,m,r,s
      integer Max_pfps,Max_cfps
      integer a,b,c,d
      integer aa1,aa2,bb1,bb2
      integer nshells,nbasis,zmax,intmax,Nprims

      integer nints,maxblk,Ncenters
      integer nalpha_pack,npcoeff_pack
      integer ncsum,next,nfirst
      integer nalpha,npcoef
      integer matm,natm,ratm,satm

      integer imin, zmin, iblk, zblk, iatm
      integer ivangmom(Nshells),ixalpha(Nshells),ixpcoef(Nshells)
      integer npfps(Nshells),ncfps(Nshells),indx_cc(Nshells)
      integer atm_4shell(Nshells)

      logical spherical
      logical degenerate 

      double precision x1,y1,z1
      double precision x2,y2,z2
      double precision x3,y3,z3
      double precision x4,y4,z4
      double precision x,y,z
      double precision eri_s,eri_l

      double precision coords(3,Nshells)
      double precision coeffs(Npcoef),alpha(Nalpha)
      double precision pcoeff_pack(4*Max_cfps*Max_pfps)
      double precision alpha_pack(4*Max_pfps)
      double precision Repulsion(a1:a2,b1:b2,a1:a2,b1:b2)
      double precision Dzero

      double precision scr(Zmax),Erd_scale(Nbasis)
      integer iscr(Intmax),Erd_index(Nbasis)
      integer ccbeg(Nprims),ccend(Nprims),end_nfps(Nshells)
      integer ccbeg_pack(Max_cfps),ccend_pack(Max_cfps)

      Data Dzero /0.0D0/
c-----------------------------------------------------------------------
c   Clear the output array.
c   Find the shell blocks for which we shall loop through.
c-----------------------------------------------------------------------

      call lookup_shell(end_nfps, nshells, a1, m1)
      call lookup_shell(end_nfps, nshells, a2, m2)
      call lookup_shell(end_nfps, nshells, b1, n1)
      call lookup_shell(end_nfps, nshells, b2, n2)

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

            ratm = Atm_4shell(n)
            x3 = coords(1,ratm)
            y3 = coords(2,ratm)
            z3 = coords(3,ratm)

            call pack_coeffs(alpha,ixalpha,coeffs,ixpcoef, 
     &                       ncfps,npfps,n,m,n,m,alpha_pack, 
     &                       nalpha_pack,pcoeff_pack,
     &                       npcoeff_pack,ccbeg,ccend,
     &                       indx_cc,ccbeg_pack,ccend_pack,
     &                       Nprims,Nshells,Nalpha,Npcoef,
     &                       Max_cfps,Max_pfps)

c---------------------------------------------------------------------------
c   Calling sequence for ERD version 2.
c---------------------------------------------------------------------------
            ncsum = ncfps(n) + ncfps(m) + ncfps(n) + ncfps(m)

#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6,*) "@-entry to ERD_batch"
      Write(6,"((a,4(1x,i4)))") "m, n, r, s:", n,m,n,m
      Write(6,"(a,2(1x,I4))") "The Nalpha_pack and Npcoeff_pack:",
     &                         nalpha_pack, npcoeff_pack
      Write(6,"((a,4(1x,i4)))") "Ivanagmom:", Ivangmom(m), 
     &      Ivangmom(r), Ivangmom(m), Ivangmom(r)
      Write(6,"((a,4(1x,i4)))") "Ncfps    :", Ncfps(n), 
     &      Ncfps(m), Ncfps(n), Ncfps(m)
      Write(6,"((a,4(1x,i4)))") "Npfps    :", Npfps(n), 
     &      Npfps(m), Npfps(n), Npfps(m)
      Write(6,*)
      Write(6,"(4(1x,F15.7))") (Alpha_pack(i),i=1,Nalpha_pack)
      Write(6,*) 
      Write(6,"(4(1x,F15.7))")(Pcoeff_pack(i),i=1,npcoeff_pack)
      Write(6,*) 
      Write(6,"(4(1x,I4))")(CCbeg_pack(i),    i=1,Nalpha_pack)
      Write(6,*) 
      Write(6,"(4(1x,I4))")(CCend_pack(i),    i=1,Nalpha_pack)
      Write(6,*)
      Write(6,"(a,3(1x,I10))") " The intmax,Zmax,ncsum:",intmax, zmax,
     &                          ncsum
      Write(6,*) "The shell coordinates"
      Write(6,*)
      Write(6,"(3(1x,F10.5))") x1, y1, z1
      Write(6,"(3(1x,F10.5))") x1, y1, z1
      Write(6,"(3(1x,F10.5))") x3, y3, z3
      Write(6,"(3(1x,F10.5))") x4, y3, z3
#endif 
            call ERD__GENER_ERI_BATCH(intmax,zmax,nalpha_pack,
     &                                npcoeff_pack,ncsum,
     &                                ncfps(n),ncfps(m),ncfps(n),
     &                                ncfps(m),npfps(n),npfps(m),
     &                                npfps(n),npfps(m),
     &                                ivangmom(n),ivangmom(m), 
     &                                ivangmom(n),ivangmom(m),
     &                                x1,y1,z1,x3,y3,z3,x1,y1,z1,
     &                                x3,y3,z3,alpha_pack,
     &                                pcoeff_pack,ccbeg_pack,
     &                                ccend_pack,spherical,
     &                               .true.,iscr,nints,nfirst,
     &                                scr)
c---------------------------------------------------------------------------
c   Move the integrals into the output block.  
c---------------------------------------------------------------------------

            if (nints .gt. 0) then
               
               call move_4cent_integrals(Repulsion,a1,a2,b1,b2,
     &                                   a1,a2,b1,b2,scr(nfirst),
     &                                   aa1,aa2,bb1,bb2,aa1,aa2,
     &                                   bb1,bb2,Erd_index,
     &                                   Erd_scale,Nbasis)
            endif

         enddo   
      enddo   

#ifdef _DEBUG_LVL0
      Write(6,*) "The <mu mu|1/r12| mu nu> integral matrix"
      eri_s = 0.0
      eri_l = 0.0
      do a = a1,a2
      do b = b1,b2
      Write(6, "(5(1x,F15.10))") ((Repulsion(a,b,c,d),c=a1,a2),
     &                             d=b1,b2)
C      eri_s = eri_s + Repulsion(a,b,c,d) * Repulsion(a,b,c,d)
C      eri_l = eri_l + Repulsion(a,b,c,d) + Repulsion(a,b,c,d)
      enddo
      enddo
      Write(*,*)
C      Write(6,"(a,a,F20.10,1x,F20.10)")"Check sums of the diagonal",
C     &                      " repulsion integrals = " ,eri_s,eri_l
#endif
      return
      end

