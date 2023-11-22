


















      subroutine start_time
      implicit none

c TIMING VARIABLES
      integer stime_sec, stime_usec
      integer utime_sec, utime_usec
      integer rtime_sec, rtime_usec

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

c   o initialize timing information
      call c_rutimes(utime_sec,utime_usec,stime_sec,stime_usec)
      call c_gtod   (rtime_sec,rtime_usec)
      ame_stime_in = (1.d-6 * stime_usec) + stime_sec
      ame_utime_in = (1.d-6 * utime_usec) + utime_sec
      ame_rtime_in = (1.d-6 * rtime_usec) + rtime_sec
      ame_timed = .true.

      return
      end

c ----------------------------------------------------------------------

