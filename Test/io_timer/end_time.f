
      subroutine end_time
      implicit none

c TIMING VARIABLES
      integer stime_sec, stime_usec
      integer utime_sec, utime_usec
      integer rtime_sec, rtime_usec
      double precision         ss, st,         us, ut,         rs, rt
      integer          sh, sm,         uh, um,         rh, rm

c COMMON BLOCKS



c aces_time.com : begin

c These six times hold the timing data from aces_init (*_in) to
c aces_fin (*_out). utime and stime count the total number of user
c and system seconds since the start of the process. rtime counts
c the total number of real seconds since 1 January 1970. ame_timed
c is a logical flag set in aces_init telling aces_fin to print out
c a timing summary.


      external aces_bd_aces_time


      double precision   ame_utime_in,  ame_stime_in,  ame_rtime_in,
     &                   ame_utime_out, ame_stime_out, ame_rtime_out
      logical            ame_timed
      common /aces_time/ ame_utime_in,  ame_stime_in,  ame_rtime_in,
     &                   ame_utime_out, ame_stime_out, ame_rtime_out,
     &                   ame_timed
      save   /aces_time/

c aces_time.com : end



c ----------------------------------------------------------------------


c   o print the total execution time since aces_init
      call c_rutimes(utime_sec,utime_usec,stime_sec,stime_usec)
      call c_gtod   (rtime_sec,rtime_usec)
      ame_stime_out = (1.d-6 * stime_usec) + stime_sec
      ame_utime_out = (1.d-6 * utime_usec) + utime_sec
      ame_rtime_out = (1.d-6 * rtime_usec) + rtime_sec

      if (ame_timed) then
         st = ame_stime_out - ame_stime_in
         ut = ame_utime_out - ame_utime_in
         rt = ame_rtime_out - ame_rtime_in
         call ds2hms(st,sh,sm,ss)
         call ds2hms(ut,uh,um,us)
         call ds2hms(rt,rh,rm,rs)
CSSS         if (nprocs.eq.1) then
            print '(/)'
            print '(3a)', "ACES ---SYSTEM--- (total sec)",
     &                        " ----USER---- (total sec)",
     &                         " -WALLCLOCK- (total sec)"
 10         format (a,2(1x,i3,":",i2,":",f5.2," (",f9.1,")"),
     &                 (1x,i3,":",i2,":",f4.1," (",f9.1,")") )
            print 10, "TIME",sh,sm,ss,st,uh,um,us,ut,rh,rm,rs,rt
CSSS         else
CSSS            if (irank.eq.0) then
CSSS 11            format (a,(1x,i3,":",i2,":",f4.1," (",f9.1,")") )
CSSS               print 11,
CSSS     &       "Total wallclock time for the root process was",rh,rm,rs,rt
CSSS            end if
CSSS         end if
      end if

      return
c     end subroutine aces_fin
      end

