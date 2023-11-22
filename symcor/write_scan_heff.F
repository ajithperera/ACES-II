      subroutine write_scan_heff(numvib, vector, scr, nscr)
c     
c     In this subroutine the SCANHEFF file is written for each electronic structure method.
c     This file can be processed further by the program vibron (using option process_scan)
c     
      implicit none

#include <aces.h>
#include <vibron.com>
#include <machsp.com>
#include <mrcc_heff.com>
#include "flags.h"
#include "jodaflags.com"
c     
      integer nscr, numvib
      double precision vector(numvib), scr(nscr)
c     
      integer min_scan, max_scan, npoint, iunit, iunit2,
     $     i, icase, ndim, ndim2, iref, ipnt, ione,
     $     iunit3, idata, idiag, ioffq, j
      double precision scan_step, reference_e, xnorm
      logical scan_info_prsnt
c     
      ione = 1
      nmode = numvib
c     
 888  format(A50,A8)
 895  format(A50,F16.4)
c     
c     get some of the basic info from SCAN_INFO file
c     
      INQUIRE (FILE='SCAN_INFO', EXIST=SCAN_INFO_PRSNT)
      if (SCAN_INFO_PRSNT) then
         write(6,*) '  @write_scan_heff: Read in Scan info'
      else
         write(6,*) ' At present Scan calculations require '
         write(6,*) ' Scan_info file '
         write(6,*) ' Unexpected error '
         call aces_exit(1)
      endif
c     
      iunit = 41
c     
c     SCAN_INFO should contain the following:
c     min_scan, max_scan, scan_step 
c     nmode
c     reference energy in a.u. 
c     displacement vector along reduced normal 
c     coordinates of parent state. One entry per line.
c     
      open(iunit, file='SCAN_INFO')
      rewind(iunit)
      read(iunit, *) min_scan, max_scan, scan_step
      read(iunit,*) numvib
      if (numvib .ne. nmode) then
         write(6,*) ' # of modes on SCAN_INFO does not match NMODE'
         write(6,*) ' error in stepts4_scan_pes.F '
         call errex
      endif
      read(iunit,*) reference_e
      do i = 1, nmode
         read(iunit,*) vector(i)
      enddo
c     
c     vector(i) contains the scan direction. It is supposed to be normalized.
c     
      xnorm = 0.0d0
      do i = 1, nmode
         xnorm = xnorm + vector(i)**2
      enddo
C     
      if (abs(xnorm) .lt. 1.0d-2) then
         write(6,*) ' Input norm of scanning vector is too low'
         write(6,*) ' check SCAN_INFO file'
         call errex
      endif
      if (abs(xnorm - 1.0d0) .gt. 1.0d-10) then
         write(6,*) ' Input vector from SCAN_INFO will be normalized'
         xnorm = 1.0d0 / sqrt(xnorm)
         do i = 1, nmode
            vector(i) = vector(i) * xnorm
         enddo
      endif
c
      write(6,*) ' min_scan, max_scan, scan_step ',
     $       min_scan, max_scan, scan_step 
c     
#ifdef _DEBUG_LVL0
      write(6,*) ' displacements from SCAN_INFO'
      call output(vector, 1, nmode, 1, 1, nmode, 1, 1)
#endif
c     
      npoint = max_scan - min_scan + 1
      if (min_scan .lt. 0 .and. max_scan .gt. 0) then
         npoint = npoint - 1
      endif
c     
      iunit2 = 185
      open (unit=iunit2,file='SCANHEFF',form='formatted')
      rewind (iunit2)
      iunit3 = 186
      open (unit=iunit3,file='SCANPLOT',form='formatted')
      rewind (iunit3)
c     
      do 10 icase = 1, ncases
c     
         ndim = 0
         call getrec(-1, 'JOBARC', ndimheff(icase), ione, ndim)
         write(6,805) ndimheff(icase), ndim
 805     format(' Process Scanned Heff ', A12, /,
     $        '   # of electronic states: ', I10)
         if (ndim .eq. 0) goto 10
c     
         ndim2 = ndim * ndim
         iref = 1
         idata = iref + ndim2
         idiag = idata + ndim2 + 1
         do i = 0, ndim2 - 1
            scr(idiag+i) = 0.0d0
         enddo
         ioffq = idiag
         do i = 1, ndim
            scr(ioffq) = reference_e
            ioffq = ioffq+ndim+1
         enddo

         call getrec(20,'JOBARC',refheff(icase),
     $        ndim2*iintfp, scr(iref))
         do i = 1, ndim2
            if (abs(scr(iref+i-1)) .lt. 1.0d-10) scr(iref+i-1)=0.0d0
         enddo
c     
         write(iunit2,888) ' Heff points for method ', nameheff(icase)
         write(iunit2,888) ' Reference hamiltonian'
         do i = 1, ndim2
            write(iunit2,887) scr(iref+i-1)
         enddo
         write(iunit2,889) min_scan, max_scan, ndim2, npoint, scan_step
         write(iunit2,888) ' Normal mode direction '
         write(iunit2, 889) nmode
         call put_matrix(vector, 1, nmode, 8, iunit2, .false.)
 889     format(4I12, E12.4)
 887     format(F25.12)
c     
c     also write info to SCANPLOT file
c     
 886     format(500E14.6)
         write(iunit3,888) ' Heff points for method ', nameheff(icase)
         write(iunit3,888) ' Reference hamiltonian'

         scr(idata) = 0.0d0
         call dcopy(ndim2, scr(iref), 1, scr(idata+1), 1)
         call daxpy(ndim2, -1.0d0, scr(idiag), 1,
     $              scr(idata+1), 1)
         call dscal(ndim2, 27.2113956d0, scr(idata+1),1)
         write(iunit3,886) (scr(idata+j), j=0,ndim2)
c     
         ipnt = 1
         call getrec(20,'JOBARC',pntheff(icase),
     $        ndim2*npoint*iintfp, scr(ipnt))
         do i = 1, ndim2*npoint
            if (abs(scr(ipnt+i-1)) .lt. 1.0d-10) scr(ipnt+i-1)=0.0d0
         enddo
c     
         do i = 1, ndim2*npoint
            write(iunit2,887) scr(ipnt+i-1)
         enddo
c     
         idata = ipnt + ndim2*npoint
         idiag = idata + ndim2 + 1
         do i = 0, ndim2 - 1
            scr(idiag+i) = 0.0d0
         enddo
         ioffq = idiag
         do i = 1, ndim
            scr(ioffq) = reference_e
            ioffq = ioffq+ndim+1
         enddo
c     
         write(iunit3,895) ' Other points', scan_step
         ioffq = ipnt
         do i = min_scan, max_scan, 1
            if (i .ne. 0) then
               scr(idata) = i * scan_step
               call dcopy(ndim2, scr(ioffq), 1, scr(idata+1), 1)
               call daxpy(ndim2, -1.0d0, scr(idiag), 1,
     $              scr(idata+1), 1)
c               write(6,*) ' current step i before dscal', i, i*scan_step
c               call output(scr(idata), 1, 1, 1, ndim2+1, 1, ndim2+1, 1)
               call dscal(ndim2, 27.2113956d0, scr(idata+1),1)
c               write(6,*) ' current step i', i, i*scan_step
c               call output(scr(idata), 1, 1, 1, ndim2+1, 1, ndim2+1, 1)
               write(iunit3,886) (scr(idata+j),j=0,ndim2)
               ioffq = ioffq + ndim2
            else
               write(iunit3,888) 'zero point'
            endif
         enddo
c     
c     process transition moments
c     
         iref = 1
c     
         call getrec(20,'JOBARC',reftmom(icase),
     $        ndim*3*iintfp, scr(iref))
         do i = 1, ndim*3
            if (abs(scr(iref+i-1)) .lt. 1.0d-10) scr(iref+i-1)=0.0d0
         enddo
c     
         write(iunit2,888) ' TMOM points for method ', nameheff(icase)
         write(iunit2,888) ' Reference transition moments'
         do i = 1, 3*ndim
            write(iunit2,887) scr(iref+i-1)
         enddo
c     
         write(iunit2,889) min_scan, max_scan, ndim*3, npoint, scan_step
c     
         write(iunit3,888) ' TMOM points for method ', nameheff(icase)
         write(iunit3,888) ' Reference transition moments'
         idata = iref + 3*ndim
         scr(idata) = 0.0d0
         call dcopy(3*ndim, scr(iref), 1, scr(idata+1), 1)
         write(iunit3,886) (scr(idata+j), j=0,ndim*3)

         ipnt = 1
         call getrec(20,'JOBARC',pnttmom(icase),
     $        ndim*3*npoint*iintfp, scr(ipnt))
         do i = 1, ndim*3*npoint
            if (abs(scr(ipnt+i-1)) .lt. 1.0d-10) scr(ipnt+i-1)=0.0d0
         enddo
c     
         do i = 1, ndim*3*npoint
            write(iunit2,887) scr(ipnt+i-1)
         enddo
c     
         write(iunit2,888) ' End point info ', nameheff(icase)

c     
         idata = ipnt + ndim*3*npoint
c     
         write(iunit3,888) ' Other points'
         ioffq = ipnt
         do i = min_scan, max_scan
            if (i .ne. 0) then
               scr(idata) = i * scan_step
               call dcopy(3*ndim, scr(ioffq), 1, scr(idata+1), 1)
               write(iunit3,886) (scr(idata+j), j=0,3*ndim)
               ioffq = ioffq + 3*ndim
            else
               write(iunit3,888) ' '
            endif
         enddo
c     
         write(iunit3,888) ' End point info ', nameheff(icase)
c     
         write(6,808) nameheff(icase)
 808     format(' Finished processing case :', A20)
c     
 10   continue
c     
      close(unit=iunit2, status = 'KEEP')
      close(unit=iunit3, status = 'KEEP')
      write(6,*) ' All done write_scanheff '
c     
      return
      end
