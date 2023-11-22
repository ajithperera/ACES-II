      subroutine zclock (string,i)
C
C     The Zclock routine prints the elapsed time.
C     If I=0, the wall clock is recorded for future use.
C     If I=1, the current timing is printed.
C
C#include <aces.h>

      double precision cpu0, cpu, sc, wall, sw, percent,
     &    timein,timenow,timetot,timenew
      character*8 string
      integer hc, mc, hw, mw
      integer i, y0,mo0,d0,h0,mi0, y,mo,d,h,mi,s,s0
      logical set_up
      integer days(12),day0,day1
c
      common /timeinfo/ timein,timenow,timetot,timenew
c
      data days /0,31,59,90,120,151,181,212,243,273,304,334/
      data set_up/.false./
c
      save y0,mo0,d0,h0,mi0,s0, cpu0      
c
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
         print 9000,string,hc,mc,sc,hw,mw,int(sw),percent
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
         print 9000,string,hc,mc,sc,hw,mw,int(sw),percent
         call flush(6) 
 9000    format(A8, I4,' h ',I2,' min ',F5.2,' sec;  ',
     .   'real time:', I4,' h ',I2,' min ',i2,' sec  (',F5.1,'%)')
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
      subroutine wallclock(year,month,day,hour,minute,second)

C
C     This subroutine should return the current time information
C     in its parameters. Add another IFDEF group for your
C     architecture/OS and enjoy!
C
c#include <aces.h>
c
      integer year,month,day,hour,minute,second
      integer*4 iarray (3)

c#ifdef M_RS6000      
c     Please note the underscores!
c      integer*4 tarray(9)
c      integer*4 currtime,time_
c      currtime=time_()
c      call ltime_(currtime, tarray)
c      second=tarray(1)
c      minute=tarray(2)
c      hour=tarray(3)
c      day=tarray(4)
c      month=tarray(5)
c      year=tarray(6)+1900
c      day of week: tarray(7)
c      day of year: tarray(8)
c      daylight savings time: tarray(9)
c#endif

c#ifdef M_SPARC

c      call idate(iarray)
      day=iarray(1)
      month=iarray(2)
      year=iarray(3)
c      call itime(iarray)
      hour=iarray(1)
      minute=iarray(2)
      second=iarray(3)
c#endif

c#ifdef M_CRAY
c      character*8 cdate, ctime
c      integer conv_char_int
c      call date(cdate)
c      month = 10*conv_char_int(cdate(1:1))
c     &    + conv_char_int(cdate(2:2))
c      day = 10*conv_char_int(cdate(4:4))
c     &    + conv_char_int(cdate(5:5))
c      year = 10*conv_char_int(cdate(7:7))
c     &    + conv_char_int(cdate(8:8))
c      call clock(ctime)
c      hour = 10*conv_char_int(ctime(1:1))
c     &    + conv_char_int(ctime(2:2))
c      minute = 10*conv_char_int(ctime(4:4))
c     &    + conv_char_int(ctime(5:5))
c      second = 10*conv_char_int(ctime(7:7))
c     &    + conv_char_int(ctime(8:8))
c#endif

c#ifdef M_SGI
c      integer iarray (3)
c      call idate (month, day, year)
c      call itime (iarray)
c      hour=iarray(1)
c      minute=iarray(2)
c      second=iarray(3)
c#endif

c#ifdef M_IP22
c      integer iarray (3)
c      call idate (month, day, year)
c      call itime (iarray)
c      hour=iarray(1)
c      minute=iarray(2)
c      second=iarray(3)
c#endif

CC#ifdef M_HPUX
C      integer iarray (3)
C      call idate(month,day,year)
C      year=year+1900
C      call itime(iarray)
C      hour=iarray(1)
C      minute=iarray(2)
C      second=iarray(3)
CC..#ELSE
C      This subroutine needs FPP processing.
C#ENDIF
      return
      end
 
c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End: c
