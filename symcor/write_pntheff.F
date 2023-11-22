      subroutine write_pntheff(scr)
c     
c     In this subroutine the PNTHEFF file is written for each electronic structure method.
c     This file can be processed further by the program vibron (using option refine_coupling)
c     
      implicit none

#include <aces.h>
#include <vibron.com>
#include <machsp.com>
#include <mrcc_heff.com>
#include "flags.h"
#include "jodaflags.com"
c     
      integer i_sided, npoint, npoint2, nmode2, ijunk, ione,
     $     irp, irpi, irpj, iunit2, ndim, ndim2,
     $     iref, ipnt, ipnt2, i, nspecial, nrest, ioff, ioff2,
     $     icount, icase
      M_REAL stpsiz, scr(*)
      character*4 doit
      character*6 ndimheff0, nameheff0
c     
 888  format(A50,A20)
c     
c     get some of the basic info from Jobarc / Flags
c     
      ione = 1
      i_sided = iflags2(h_IFLAGS2_grid_vibron)
      CALL GETREC(20,'JOBARC','NUMPOINT',IONE,npoint)
c     
      CALL GETREC(-1,'JOBARC','DANGERUS',IONE,IJUNK)
      if(iflags(79).eq.0.and.ijunk.eq.0)doit='FULL'
      if(iflags(79).eq.1.or.ijunk.eq.1)doit='COMP'
      CALL GETREC(20,'JOBARC',DOIT//'NIRX',IONE,nirrepf)
      call getrec(20, 'JOBARC', 'NMODES_F', nirrepf, numirpf)
c
      if (IFLAGS(57) .lt. 0) then
         write(6,*)  ' We used displacements along ',
     $        'dimensionless reduced coordinates.'
         stpsiz = DFLOAT(IFLAGS(57)) / 1000.0d0
         write(6,511) stpsiz
 511     format(T3,' Elementary stepsize ', F8.3, ' Dimensionless ')
      else
         STPSIZ=DFLOAT(IFLAGS(57))*10.0D-5
         write(6,*) ' Use Normal mode displacements '
         WRITE(6,500)STPSIZ
 500     FORMAT(T3,'Step size ',F8.5,' amu**(1/2) * bohr.')
      endif
c     
      write(6,*) ' @write_pntheff : normal mode information'
      write(6,800) nirrepf
 800  format('  Number of symmetries : ', I8)
      write(6,*) ' Number of normal modes per symmetry block '
      write(6,*) '   IRREP       # of Modes '
      do irp = 1, nirrepf
         write(6,801) irp, numirpf(irp)
      enddo
 801  format(2I8)
c     
      write(6,*) 
c     
      nmode = 0
      do irp = 1, nirrepf
         nmode = nmode + numirpf(irp)
      enddo
      write(6,802) nmode
 802  format( '  Total number of normal modes ', I8)
c     
      nmode2 = 0
      if (i_sided .ge. 3 .and. i_sided .ne. 5) then
         do irp = 1, nirrepf
            nmode2 = nmode2 + numirpf(irp)*(numirpf(irp)-1) / 2
         enddo
      endif
      if (i_sided .eq. 4 .or. i_sided .eq. 7) then
         do irpi = 1, nirrepf
            do irpj = 1, irpi - 1
               nmode2 = nmode2 + numirpf(irpi) * numirpf(irpj)
            enddo
         enddo
      endif
c     
c     calculate full set of points, including points that we know by symmetry.
c     
      npoint2 = nmode
      if (i_sided .ge. 2) then
         npoint2 = npoint2 + nmode
      endif
      if (i_sided .ge. 3 .and. i_sided .ne. 5) then
         npoint2 = npoint2 + 2*nmode2
      endif
      if (i_sided .ge. 5) then
         npoint2 = npoint2 + 2 * nmode
      endif
c
      if (i_sided .eq. 8) then
         npoint2 = npoint
         nmode2 = npoint - 2 * nmode
      elseif (i_sided .ge. 9) then
         npoint2 = npoint
         nmode2 = npoint - 4 * nmode
      endif
c
      write(6,804) nmode2
 804  format(' Total number of quadratic displacements ', I8)
c     
 806     format(' Unique number of points: ', I8)
         write(6,806) npoint
         write(6,807) npoint2
 807     format( ' Expanded number of points: ', I8)

      nspecial = nmode
      nrest = nmode - nspecial

      iunit2 = 185
      open (unit=iunit2,file='PNTHEFF',form='formatted')
      rewind (iunit2)
c     
c write parent energy array in usual format
c
c
         ndimheff0 = 'NEFF_0'
         nameheff0 = 'HEFF_0'
         ndim = 1
         write(6,805) ndimheff0, ndim
c     
         ndim2 = ndim * ndim
         iref = 1
         call getrec(20,'JOBARC','REFPAR_E',
     $        iintfp, scr(iref))
c     
         write(iunit2,888) ' Heff points for method ', nameheff0
         write(iunit2,888) ' Reference hamiltonian'
         do i = 1, ndim2
            write(iunit2,887) scr(iref+i-1)
         enddo
         write(iunit2,889) i_sided, ndim2, nmode2, npoint2, stpsiz
c     
         ipnt = 1
         ipnt2 = ipnt + ndim2*npoint
         call getrec(20,'JOBARC','PNTENERG',
     $        npoint*iintfp, scr(ipnt))
c     
            do i = 1, ndim2*npoint
               write(iunit2,887) scr(ipnt+i-1)
            enddo
c
c write fake transition moments for convenience in VIBRON
c
c     
c     process transition moments
c     
         iref = 1
c
         call dzero(scr(iref), 3*ndim)
c         call getrec(20,'JOBARC',reftmom(icase),
c     $        ndim*3*iintfp, scr(iref))
c     
         write(iunit2,888) ' TMOM points for method ', nameheff0
         write(iunit2,888) ' Reference transition moments'
         do i = 1, 3*ndim
            write(iunit2,887) scr(iref+i-1)
         enddo
c     
         write(iunit2,889) i_sided, ndim*3, nmode2, npoint2, stpsiz

         ipnt = 1
         ipnt2 = ipnt + ndim*3*npoint
         call dzero(scr(ipnt), 3*ndim*npoint)
c         call getrec(20,'JOBARC',pnttmom(icase),
c     $        ndim*3*npoint*iintfp, scr(ipnt))
c     
         do i = 1, ndim*3*npoint
            write(iunit2,887) scr(ipnt+i-1)
         enddo
c
         write(iunit2,888) ' End point info ', nameheff0
c     
         write(6,808) nameheff0
c
c     done with parent energy
c
      do 10 icase = 1, ncases
c     
         ndim = 0
         call getrec(-1, 'JOBARC', ndimheff(icase), ione, ndim)
         write(6,805) ndimheff(icase), ndim
 805     format(' Process vibronic_coupling ', A12, /,
     $        '   # of electronic states: ', I10)
         if (ndim .eq. 0) goto 10
c     
         ndim2 = ndim * ndim
         iref = 1
         call getrec(20,'JOBARC',refheff(icase),
     $        ndim2*iintfp, scr(iref))
c     
         write(iunit2,888) ' Heff points for method ', nameheff(icase)
         write(iunit2,888) ' Reference hamiltonian'
         do i = 1, ndim2
            write(iunit2,887) scr(iref+i-1)
         enddo
         write(iunit2,889) i_sided, ndim2, nmode2, npoint2, stpsiz
 889     format(4I12, E12.4)
 887     format(F25.12)
c     
         ipnt = 1
         ipnt2 = ipnt + ndim2*npoint
         call getrec(20,'JOBARC',pntheff(icase),
     $        ndim2*npoint*iintfp, scr(ipnt))
c     
         if (npoint .ne. npoint2) then
            ioff = ipnt
            ioff2 = ipnt2
            icount = 1
            write(6,*) ' copy from point , # ', icount, nmode
            call B_COPY(nmode*ndim2, scr(ioff), 1,
     $           scr(ioff2), 1)
            ioff2 = ioff2 + nmode*ndim2
            ioff = ioff + nmode*ndim2
            icount = icount + nmode
            if (i_sided .ge. 2) then
               write(6,*) ' copy from point , # ', icount, nspecial
               call B_COPY(nspecial*ndim2, scr(ioff), 1,
     $              scr(ioff2), 1)            
               ioff2 = ioff2 + nspecial * ndim2
               ioff = ioff + (nspecial-nmode) * ndim2
               icount = icount + (nspecial -nmode)
               write(6,*) ' copy from point , # ', icount, nrest
               call B_COPY(nrest*ndim2, scr(ioff), 1,
     $              scr(ioff2), 1)            
               ioff2 = ioff2 + nrest * ndim2
               ioff = ioff + nmode * ndim2
               icount = icount + nmode
            endif
            if (i_sided .ge. 3 .and. i_sided .ne. 5) then
               write(6,*) ' copy from point , # ', icount, 2*nmode2
               call B_COPY(2*nmode2*ndim2, scr(ioff), 1,
     $              scr(ioff2), 1)            
               ioff2 = ioff2 + 2*nmode2 * ndim2
               ioff = ioff + 2*nmode2 * ndim2
               icount = icount + 2*nmode2
            endif
            if (i_sided .ge. 5) then
               write(6,*) ' copy from point , # ', icount, nmode
               call B_COPY(nmode*ndim2, scr(ioff), 1,
     $              scr(ioff2), 1)            
               ioff2 = ioff2 + nmode * ndim2
               ioff = ioff + nmode*ndim2
               icount = icount + nmode
               write(6,*) ' copy from point , # ', icount, nspecial
               call B_COPY(nspecial*ndim2, scr(ioff), 1,
     $              scr(ioff2), 1)            
               ioff2 = ioff2 + nspecial * ndim2
               ioff = ioff + (nspecial-nmode) * ndim2
               icount = icount + (nspecial - nmode)
               write(6,*) ' copy from point , # ', icount, nrest
               call B_COPY(nrest*ndim2, scr(ioff), 1,
     $              scr(ioff2), 1)            
               ioff2 = ioff2 + nrest * ndim2
               ioff = ioff + nmode * ndim2
               icount = icount + nmode
            endif            
            write(6,*) ' final icount, npoint ', icount, npoint
            if (ioff .ne. ipnt + npoint*ndim2) then
               write(6,*) ' Something wrong in calc_couplings a'
               call aces_exit(1)
            endif
            if (ioff2 .ne. ipnt2 + npoint2*ndim2) then
               write(6,*) ' Something wrong in calc_couplings b'
               call aces_exit(1)
            endif
            do i = 1, ndim2*npoint2
               write(iunit2,887) scr(ipnt2+i-1)
            enddo
         else 
            do i = 1, ndim2*npoint
               write(iunit2,887) scr(ipnt+i-1)
            enddo
         endif
c     
c     process transition moments
c     
         iref = 1
c     
         call getrec(20,'JOBARC',reftmom(icase),
     $        ndim*3*iintfp, scr(iref))
c     
         write(iunit2,888) ' TMOM points for method ', nameheff(icase)
         write(iunit2,888) ' Reference transition moments'
         do i = 1, 3*ndim
            write(iunit2,887) scr(iref+i-1)
         enddo
c     
         write(iunit2,889) i_sided, ndim*3, nmode2, npoint2, stpsiz

         ipnt = 1
         ipnt2 = ipnt + ndim*3*npoint
         call getrec(20,'JOBARC',pnttmom(icase),
     $        ndim*3*npoint*iintfp, scr(ipnt))
c     
         if (npoint .ne. npoint2) then
            ioff = ipnt
            ioff2 = ipnt2
            call B_COPY(nmode*ndim*3, scr(ioff), 1,
     $           scr(ioff2), 1)
            ioff2 = ioff2 + nmode*ndim*3
            ioff = ioff + nmode*ndim*3
            if (i_sided .ge. 2) then
               call B_COPY(nspecial*ndim*3, scr(ioff), 1,
     $              scr(ioff2), 1)            
               ioff2 = ioff2 + nspecial * ndim*3
               ioff = ioff + (nspecial-nmode) * ndim*3
               call B_COPY(nrest*ndim*3, scr(ioff), 1,
     $              scr(ioff2), 1)            
               ioff2 = ioff2 + nrest * ndim*3
               ioff = ioff + nmode * ndim*3
            endif
            if (i_sided .ge. 3 .and. i_sided .ne. 5) then
               call B_COPY(2*nmode2*ndim*3, scr(ioff), 1,
     $              scr(ioff2), 1)            
               ioff2 = ioff2 + 2*nmode2 * ndim*3
               ioff = ioff + 2*nmode2 * ndim*3
            endif
            if (i_sided .ge. 5) then
               call B_COPY(nmode*ndim*3, scr(ioff), 1,
     $              scr(ioff2), 1)            
               ioff2 = ioff2 + nmode * ndim*3
               ioff = ioff + nmode*ndim*3
               call B_COPY(nspecial*ndim*3, scr(ioff), 1,
     $              scr(ioff2), 1)            
               ioff2 = ioff2 + nspecial * ndim*3
               ioff = ioff + (nspecial-nmode) * ndim*3
               call B_COPY(nrest*ndim*3, scr(ioff), 1,
     $              scr(ioff2), 1)            
               ioff2 = ioff2 + nrest * ndim*3
               ioff = ioff + nmode * ndim*3
            endif            
            if (ioff .ne. ipnt + npoint*ndim*3) then
               write(6,*) ' Something wrong in calc_couplings a'
               call aces_exit(1)
            endif
            if (ioff2 .ne. ipnt2 + npoint2*ndim*3) then
               write(6,*) ' Something wrong in calc_couplings b'
               call aces_exit(1)
            endif
            do i = 1, ndim*3*npoint2
               write(iunit2,887) scr(ipnt2+i-1)
            enddo
         else 
            do i = 1, ndim*3*npoint
               write(iunit2,887) scr(ipnt+i-1)
            enddo
         endif
c     
         write(iunit2,888) ' End point info ', nameheff(icase)
c     
         write(6,808) nameheff(icase)
 808     format(' Finished processing case :', A20)
c     
 10   continue
c     
      write(6,*) ' All done write_pntheff before close'
      close(unit=iunit2, status = 'KEEP')
      write(6,*) ' All done write_pntheff '
c     
      return
      end
