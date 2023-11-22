      subroutine zclock (string,i)
C
C     The Zclock routine prints the elapsed time.
C     If I=0, the wall clock is recorded for future use.
C     If I=1, the current timing is printed.
C
C#include <aces.h>

      double precision cpu0, cpu, sc, wall, sw, percent,
     &    timein,timenow,timetot,timenew,string
      integer hc, mc, hw, mw
      integer i, y0,mo0,d0,h0,mi0, y,mo,d,h,mi,s,s0
      logical set_up,print
      integer days(12),day0,day1
c
      common /timeinfo/ timein,timenow,timetot,timenew
      common/flags/iflags(100)
c
      data days /0,31,59,90,120,151,181,212,243,273,304,334/
      data set_up/.false./
c
      save y0,mo0,d0,h0,mi0,s0, cpu0      
c
c      print=iflags(1).gt.11
      print=.false.
      if (i.ne.0 .and. .not. set_up) then
         print*,' Warning: Zclock called with I=1 before set up'
         return
      end if
      if (i.eq.0) then
         call wallclock (y0,mo0,d0,h0,mi0,s0)
         call timer(1)
         cpu0 = timenow
c        you could do other initializations here, if needed.
         set_up=.true.
      else if (i.eq.1) then
         call wallclock (y,mo,d,h,mi,s)
c        Runs longer than one year are unprobable and not supported.
         day0=days(mo0)+d0
         day1=days(mo )+d
         if (y0.ne.y) day1=day1+365
         if (mo0.gt.2 .and. (mod(y0,4).eq.0)) day0=day0+1
         if (mo .gt.2 .and. (mod(y ,4).eq.0)) day1=day1+1
         wall=(day1-day0)*86400+(h-h0)*3600+(mi-mi0)*60+(s-s0)
         call timer(1)
         cpu = timenow - cpu0
         call hms (cpu ,hc,mc,sc)
         call hms (wall,hw,mw,sw)
         if (wall.gt.0.0) then
            percent=cpu/wall*100.0
         else
            percent=100.0
         endif
         if(print)then
         print 9000,string,hc,mc,sc,hw,mw,int(sw),percent
         endif
         call flush(6) 
 9000    format(A8, I4,' h ',I2,' min ',F5.2,' sec;  ',
     .   'real time:', I4,' h ',I2,' min ',I2,' sec  (',F5.1,'%)')
c        CPU TIME: 999 H 23 MIN 65.32 SEC   REAL TIME: 999 H 23 MIN 76
c        SEC  (100.0%)
      end if
      return
      end
c
c     Auxiliary routine to split a time into h:m:s
c
      subroutine hms (time,h,m,s)
c#include <aces.h>      
c      
      double precision
     &    time,s
      integer h,m
      
      h=time/3600.0
      m=(time-3600*h)/60.0
      s=time-3600*h-60*m
      return
      end
c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End: c
