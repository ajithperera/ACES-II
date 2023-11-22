c------------------------------------------------------------------------------
c     Modifications for "Release 2" : There are better mechanisms are in place
c     now to document the changes from version to another also version control
c     is more systematic (Anthony Yau is responsible for version control and
c     documentation).
c
c     NO  3/29/94  HF-DFT single points
c     JDW 3/31/94  5th order energy corrections, tdhf, EOM-EA single points
c     JDW 9/25/95  KKB's' dropped core derivatives
c     SG  3/4/96   P-EOM-MBPT(2) gradients
c     AP  12/07/99 KS-DFT
c     AP  07/2000  seward, alaska, and mckinley
c------------------------------------------------------------------------------
c
#include "flags.h"

      Program ACES2

#include "aces2.h"

      integer t0, t1, itmp,iIRC_search
      logical Ignore,IRC_present 
      character*80 szJOBARC,szIRC_search

c   o gettimeofday
      call c_gtod(t0,itmp)
      dTimeStart = t0 + 1.d-6*itmp

c   o start up the parallel evironment
      call parenv_start

c   o print the ACES banner
      call title
      call flushout
C
C Check whether we are doing a manual parallled finite differnce. Well,
C if there is a JOBARC, this could be a manual finite difference or 
C a restart. This logic related to IRC_scan is added on 04/2019.
C IRC searches also have a JOBARC, but can not skip runing joda. 
C 
      call gfname('JOBARC',szJOBARC,iJOBARC)
      call gfname('IRC_SCAN',szIRC_search,iIRC_search)
      inquire(file=szJOBARC(1:iJOBARC),exist=ignore)
      inquire(file=szIRC_SEARCH(1:iIRC_search),exist=IRC_present)
      if (Ignore) Manual_FD = .TRUE. 
      if (IRC_present) Manual_FD = .FALSE.

      If (Manual_FD  .AND. IRC_present) Then
         Write(6,"(a)") " Do not combine manual finite difference and",
     &                  " the IRC path following"
         Call Errex
      Endif 

      call runit('rm -f FILES')
      if (.NOT. Manual_FD) call runit(xjoda)

c   o initialize the ACES environment and turn off the job archive subsystem
      call crapsi(ijunk,iuhf,-1)
      call aces_ja_fin

      call init
      if (iflags2(113).eq.1.and.
     &    iflags(3)   .eq.1.and.
     &    iflags(68)  .eq.1     ) then
         call polyrate
         write(*,*)
         write(*,*) '@ACES2: The Polyrate routine has',
     &              ' completed successfully.'
         write(*,*)
         call c_exit(0)
      end if
c
c   o the functionalities in mrcc branch is developed by Marcel Nooijen
c
      if (mrcc) then
         write(6,*) ' @aces2.main: enter mrcc_branch '
         call mrcc_branch
       else if (iflags2(h_IFLAGS2_prep_vibron).eq.1 .or.
     $         iflags2(h_IFLAGS2_geo_sample).eq.1 .or.
     $         iflags2(h_IFLAGS2_gen_quartic).eq.1 .or.
     $          iflags2(h_IFLAGS2_scan_normal).eq.1 ) then
         write(6,*) ' @aces2.main: enter aces2_heff'
         call aces2_heff
      else
         call Job_control
      end if

c   o shut down the parallel evironment
      call parenv_stop

      call c_gtod(t1,itmp)
      print '(/)'
      print *, '@ACES2: ',
     &         'The ACES2 program has completed successfully in ',
     &         t1+1-t0,' seconds.'
      print '(/)'

      call c_exit(0)
      end

