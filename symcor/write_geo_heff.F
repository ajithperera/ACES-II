      subroutine write_geo_heff(scr, nscr)
c     
c     In this subroutine the sample_HEFF file is written for each electronic structure method.
c     This file can be processed further by the program vibron (using option sample_surface)
c     
      implicit none

#include <aces.h>
#include <vibron.com>
#include <machsp.com>
#include <mrcc_heff.com>
#include "flags.h"
#include "jodaflags.com"
c     
      integer nscr
      double precision  scr(nscr)
c     
      integer nradial, nequi, nsample, npoint, iunit, iunit2,
     $     i, icase, ndim, ndim2, iref, ipnt, ione,
     $     iunit3, idata, idiag, ioffq, j, ievec, ieval
      double precision scan_step, reference_e, xnorm
      logical scan_info_prsnt
c     
      ione = 1
c     
 888  format(A50,A8)
 895  format(A50,F16.4)
c     
c     get some of the basic info from GEO_SAMPLE file
c     
      iunit = 41
      write(6,*) 'Entered write_geo_heff '
c     
      open(iunit, file='GEO_SAMPLE')
      rewind(iunit)
      read(iunit,*)
      read(iunit, *) nradial, nequi, nsample
c
      close(iunit,status='KEEP')
      call getrec(1, 'JOBARC', 'REFER_E ', IINTFP,reference_e)
      write(6,*) ' reference_e from JOBARC', reference_e
c
      npoint = nsample
c     
      iunit2 = 185
      open (unit=iunit2,file='GEO_ADIAB',form='formatted')
      rewind (iunit2)
      iunit3 = 186
      open (unit=iunit3,file='GEO_HEFF',form='formatted')
      rewind (iunit3)
c     
      do 10 icase = 1, ncases
c     
         ndim = 0
         write(6,805) ndimheff(icase), ndim
         write(6,*) 'before'
         call getrec(-1, 'JOBARC', ndimheff(icase), ione, ndim)
         write(6,805) ndimheff(icase), ndim
 805     format(' Process Sampled Heff ', A12, /,
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
c     
         write(iunit2,888) ' Adiabatic energies for method ',
     $        nameheff(icase)
c     
 886     format(I6,500F12.6)
         write(iunit3,888) ' Heff points for method ', nameheff(icase)

         ipnt = 1
         call getrec(20,'JOBARC',pntheff(icase),
     $        ndim2*npoint*iintfp, scr(ipnt))
         do i = 1, ndim2*npoint
            if (abs(scr(ipnt+i-1)) .lt. 1.0d-10) scr(ipnt+i-1)=0.0d0
         enddo
c     
         idata = ipnt + ndim2*npoint
         idiag = idata + ndim2 
         ievec = idiag + ndim2
         ieval = ievec + ndim2
c
         do i = 0, ndim2 - 1
            scr(idiag+i) = 0.0d0
         enddo
         ioffq = idiag
         do i = 1, ndim
            scr(ioffq) = reference_e
            ioffq = ioffq+ndim+1
         enddo
c     
         ioffq = ipnt
         do i = 1, nsample
               call dcopy(ndim2, scr(ioffq), 1, scr(idata), 1)
               call daxpy(ndim2, -1.0d0, scr(idiag), 1,
     $              scr(idata), 1)
               call dscal(ndim2, 27.2113956d0, scr(idata),1)
               write(iunit3,886) i, (scr(idata+j-1),j=1,ndim2)
               call eig(scr(idata), scr(ievec), ndim, ndim, 1)
               call dcopy(ndim, scr(idata), ndim+1, scr(ieval), 1)
               write(iunit2,886) i, (scr(ieval+j-1),j=1,ndim)
               ioffq = ioffq + ndim2
         enddo
c     
c     process transition moments
c     
         iref = 1
c
         write(iunit3,*)
         write(iunit3,888) ' TMOM points for method ', nameheff(icase)
         write(iunit3,*)

         ipnt = 1
         call getrec(20,'JOBARC',pnttmom(icase),
     $        ndim*3*npoint*iintfp, scr(ipnt))
         do i = 1, ndim*3*npoint
            if (abs(scr(ipnt+i-1)) .lt. 1.0d-10) scr(ipnt+i-1)=0.0d0
         enddo
c     
         idata = ipnt + ndim*3*npoint
c     
         ioffq = ipnt
         do i = 1, nsample
               call dcopy(3*ndim, scr(ioffq), 1, scr(idata+1), 1)
               write(iunit3,886) i, (scr(idata+j), j=1,3*ndim)
               ioffq = ioffq + 3*ndim
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
      write(6,*) ' All done write_geo_heff '
c     
      return
      end
