#include <flags.h>
      subroutine upd_fds4_heff(NATOM, DOIT, IMORE, FD_POINTS, 
     &                         IPTTYPE, DSCR, NPOINT, NDSCR)
c     
c this routine determines the next point which must be run in a
c numerical frequency calculation, and the calculation type
c (single point energy or gradient), by inspecting information
c on the jobarc file.  
c     
      implicit double precision(a-h,o-z)
C
#include <mrcc_heff.com>
#include <machsp.com>
#include <jodaflags.com>
#include <control.com>
#include "parallel_aces.com"
 
      character*4 doit
C
      double precision  DSCR(NDSCR), FD_POINTS(3*NATOM*NPOINT)
      logical print
C
      integer ipttype(NPOINT)
c     
      data tol   /1.d-8/
      data one   /1.0/
      data zilch /0.0/
c     
      ione=1
      print = .true.
      nsize=3*natom
c   o read point type record and the list of displacements
      write(6,*) ' @upd_fds4_heff: npoint : ', npoint
      CALL GETREC(20,'JOBARC','FDCALCTP',NPOINT,IPTTYPE)
      if (print) then
      write(6,*) ' @upd_fds4_heff, ipttype'
      do i = 1, npoint
         write(6,*) i, ipttype(i)
      enddo
      endif

c Point types can be one of the following:
c  > 0 : this point will be done (1+iRank will be done by this process)
c  = 0 : this point is skipped
c  < 0 : this point is done (-1-iRank was done by this process)
c
c   o find first entry which must be calculated
      inext = 1
      do while ((ipttype(inext).ne.1+irank).and.(inext.le.npoint))
         inext = inext + 1
      end do

c   o find last entry which was calculated
      ilast = inext-1
      do while ((ipttype(ilast).ne.-1-irank).and.(ilast.gt.0))
         ilast = ilast - 1
      end do

      write(6,*) ' @upd_fds4_heff, inext  : ', inext
      if (print) then
      write(6,*) ' @upd_fds4_heff, ilast  : ', ilast
      write(6,*) ' @upd_fds4_heff, npoint : ', npoint
      endif
c
c   o prepare info for the next calculation
      IF (INEXT.NE.NPOINT+1) THEN
         CALL GETREC(20,'JOBARC','FDCOORDS',IINTFP*NSIZE*INEXT,
     &               FD_POINTS)
         ILOC=1+(INEXT-1)*NSIZE
         CALL PUTREC(20,'JOBARC','NEXTGEOM',IINTFP*NSIZE,
     &               FD_POINTS(ILOC))
         if (print) then
            write(6,*) ' NEXTGEOM in upd_fds4_heff'
            call output(FD_POINTS(iloc), 1, 1, 1, nsize, 1,
     $           nsize, 1)
         endif
         IPTTYPE(INEXT)=-IPTTYPE(INEXT)
         CALL PUTREC(20,'JOBARC','FDCALCTP',INEXT,IPTTYPE)
         IMORE=1
c      o tag the last displacement
         inext = inext+1
         do while ((ipttype(inext).ne.1+irank).and.(inext.le.npoint))
            inext = inext + 1
         end do
         if (inext.eq.npoint+1) then
            CALL PUTREC(20,'JOBARC','LASTGEOM',1,1)
         end if
      ELSE
         IMORE = 0 
      END IF

      if (Ilast .NE. 0) Then
c
c if ACES2 calculation:
c
         write(6,*) ' @upd_fds4_heff: Program flag ',
     $        IFLAGS2(h_IFLAGS2_program)
         IF (IFLAGS2(h_IFLAGS2_program).NE.3) THEN
c      o update energy vector
         call getrec(20, 'JOBARC', 'PNTENERG',
     $        npoint*iintfp, DSCR)
               iloc = 1+(ILAST-1)
         IF (IFLAGS(h_IFLAGS_excite).EQ.0) THEN
            CALL GETREC(20,'JOBARC','TOTENERG',IINTFP,DSCR(iloc))
         ELSE
            CALL GETREC(20,'JOBARC','TOTENER2',IINTFP,DSCR(iloc))
         END IF
         call putrec(20, 'JOBARC', 'PNTENERG',
     $        npoint*iintfp, DSCR)
         write(6,*) ' @upd_fds4_heff: update ACES2 PNTENERG '
         call output(dscr, 1, 1, 1, npoint, 1, npoint, 1)
         go to 999
         endif

c    
c     process information from last calculation
c     
c process the parent energy
c
         call getrec(20, 'JOBARC', 'PNTENERG',
     $        npoint*iintfp, DSCR)
               iloc = 1+(ILAST-1)
         call getrec(20, 'JOBARC', 'PARENERG',
     $        iintfp, DSCR(iloc))
         call putrec(20, 'JOBARC', 'PNTENERG',
     $        npoint*iintfp, DSCR)
c     
c loop over a variety of effective Hamiltonians and t-moments
c     
         do icase = 1, ncases
            ndim = 0
            call getrec(20, 'JOBARC', ndimheff(icase), ione, ndim)
#ifdef _DEBUG_LVL0

#endif
            if (ndim .ne. 0) then
c     
c process heff part
c
               write(6,*) ' Process icase ', icase, ndimheff(icase),
     $              ndim
#ifdef _DEBUG_LVL0
               write(6,*) ' Process icase ', icase, ndimheff(icase)
#endif
               ndim2 = ndim*ndim
               call getrec(20,'JOBARC',pntheff(icase),
     $              ndim2*npoint*iintfp, DSCR)
               iloc = 1+(ILAST-1)*ndim2
               call getrec(20,'JOBARC',nameheff(icase),
     $              ndim2*iintfp, DSCR(iloc))
#ifdef _DEBUG_LVL0
               write(6,880)
     $              nameheff(icase), ilast, iloc
 880           format(' Read Heff from Jobarc ',A12, 2I10)
               call output(DSCR(iloc), 1, ndim, 1, ndim, ndim,
     $              ndim, 1)
#endif
               call putrec(20,'JOBARC',pntheff(icase),
     $              ndim2*npoint*iintfp, DSCR) 
c     
c now process transition moments
c     
               call getrec(20,'JOBARC',pnttmom(icase),
     $              ndim*3*npoint*iintfp, DSCR)
               iloc = 1+(ILAST-1)*ndim*3
               call getrec(20,'JOBARC',nametmom(icase),
     $              ndim*3*iintfp, DSCR(iloc))
#ifdef _DEBUG_LVL0
               write(6,881) pnttmom(icase), ilast, iloc
 881           format(' Read Tmom from Jobarc ',A12, 2I10)
               call output(DSCR(iloc), 1, 3, 1, ndim, 3,
     $              ndim, 1)
#endif
               call putrec(20,'JOBARC',pnttmom(icase),
     $              ndim*3*npoint*iintfp, DSCR)
c     
            endif
         enddo
      endif
c     
c markzero indicates the position from where Jobarc can be zeroed out
c
 999  continue
      CALL PUTREC(20,'JOBARC','MARKZERO',ione, ione)
c
      return
      end
