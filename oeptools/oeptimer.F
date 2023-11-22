      subroutine oeptimer (
     &     name, job)
c
c This routine either turns timer on or turns it off and remembers the
c timing. The timing are stored in the hash, with the key given by NAME.
c
c Argument JOB can take the following values:
c   "IN"  : initialize the OEPTIMER common block
c   "ON " : turn the timer on for NAME
c   "OF"  : turn the timer off for NAME and calculate the timing
c   "RE"  : report timings for all NAMEs registred
c 
c The scope of the routine has nothing to do with OEP, the reason for
c the name is that I do not want to conflict with VSCF (or any module)
c namespace.
c     
c It's somewhat poor implementation, just want to get the thing going.
c
c Igor Schweigert, Sep 2003 ($Id: oeptimer.FPP,v 1.1.1.1 2008/05/29 19:35:39 taube Exp $)
c 
      implicit none
c
c     Arguments
c 
      character*(*) name
      character*(*) job
c
c     Global variables (stored in a common block)
c 
      integer size, last
      parameter (size=50)
      integer ntimes (size)
      double precision rtimes (size), utimes (size), stimes (size)
      character*25 names (size)
      common /oeptime/ last, ntimes, rtimes, utimes, stimes, names
c
c     Local variables
c 
      integer
     &     i, rec, r_sec, r_msec, u_sec, u_msec, s_sec, s_msec
c 
      double precision
     &     rtime, utime, stime
c
c     Messages to the user
c     
 5900 format (
     &     t3,'@OEPTIMER-W: Cannot create another record ',a,/)
 5920 format (
     &     t3,'@OEPTIMER-W: Timer has never been turned on for ',a,/)
 5940 format (
     &     t3,'@OEPTIMER-I:',/,
     &     3x,58('-'),/,
     &     6x,26x,'Wallclock',6x,'User',4x,'System',1x,' Ave.'/,
     &     3x,58('-'))
 5945 format (
     &     6x,a25,1x,f9.3,1x,f9.3,1x,f9.3,1x,i4)
 5947 format (
     &     3x,58('-'),/)
 5960 format (
     &     t3,'@OEPTIME-W: Argument JOB has illegal value ', a)
      
c
c     If job is "IN" then we initilize the OEPTIMER's common block
c     
      if (job(1:2) .eq. "IN") then
         last = 0
         do i=1, size
            ntimes (i) = 0
            rtimes (i) = 0.d0
            utimes (i) = 0.d0
            stimes (i) = 0.d0
            names (i) = ''
         enddo
c
c     If job is "ON" then we:
c       1/ get the current wallclock, user, and system times;
c       2/ check if the record for NAME already exists;
c       3/ overwrite it if it does;
c       4/ check if there is a space for a new record and create it
c     
      elseif (job(1:2) .eq. "ON") then
c
         call c_gtod (r_sec, r_msec)
         call c_rutimes (u_sec, u_msec, s_sec, s_msec)
         rtime = dble (r_sec) + 1.d-6 * dble (r_msec)
         utime = dble (u_sec) + 1.d-6 * dble (u_msec)
         stime = dble (s_sec) + 1.d-6 * dble (s_msec)
c     
         rec = 0
         do i=1, last
            if (name .eq. names (i)) rec = i
         enddo
c     
         if (rec .gt. 0) then
            ntimes (rec) = ntimes (rec) + 1
            rtimes (rec) = rtimes (rec) - rtime
            utimes (rec) = utimes (rec) - utime
            stimes (rec) = stimes (rec) - stime
         else
            if (last .ge. size) then
               write (6, 5900) name
            else
               last = last + 1
               names  (last) = name
               ntimes (last) = 1
               rtimes (last) = -rtime
               utimes (last) = -utime
               stimes (last) = -stime
            endif
         endif
c
c     If JOB is "OFF" we:
c       1/ get the current wallclock time
c       2/ check if the record for NAME exists
c       3/ if it does then substract its value from the current time
c          and store the difference, it it does not then complain and
c          die
c     
      elseif (job(1:2) .eq. "OF") then
c
         call c_gtod (r_sec, r_msec)
         call c_rutimes (u_sec, u_msec, s_sec, s_msec)
         rtime = dble (r_sec) + 1.d-6 * dble (r_msec)
         utime = dble (u_sec) + 1.d-6 * dble (u_msec)
         stime = dble (s_sec) + 1.d-6 * dble (s_msec)
c     
         rec = 0
         do i=1, last
            if (name .eq. names (i)) rec = i
         enddo
c     
         if (rec .gt. 0) then
            rtimes (rec) = rtime + rtimes (rec)
            utimes (rec) = utime + utimes (rec)
            stimes (rec) = stime + stimes (rec)
         else
            write (6, 5920) name
         endif
c     
c     If JOB is "RE" we print out the NAME, TIME values of all the
c     records.
c     
      elseif (job (1:2) .eq. "RE") then
c     
         if (last .gt. 0) then
            write (6, 5940)
            do i=1, last
               write (6, 5945)
     &              names (i),
     &              rtimes (i)/dble (ntimes (i)),
     &              utimes (i)/dble (ntimes (i)),
     &              stimes (i)/dble (ntimes (i)),
     &              ntimes (i)
            enddo
            write (6, 5947)
         endif
c     
c     If JOB is not "ON", "OF", nor "RE" then print the warning and do
c     nothing.
c     
      else
         write (6, 5960) job
      endif
c
      return 
      end
         
         
         
         
      
      

