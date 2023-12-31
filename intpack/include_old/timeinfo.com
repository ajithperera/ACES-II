
#ifndef _TIMEINFO_COM_
#define _TIMEINFO_COM_

c This common block keeps track of timing information.  It is only required
c in the crapsi and crapso routines for the time being as only the starting
c and initial times are determined.  In the future, it may be nice to find
c the time required to do different pieces of a calculation, and the data
c here is sufficiently flexible to do so.

c timein   the time of the first call to timer
c timenow  the time of the previous/current call to timer (timenow is set
c          to the time returned, so up to the actual call, timenow actually
c          contains the value from the previous call)
c timetot  total time elapsed since the first call to timer
c timenew  total time elapsed since the last call to timer

      M_REAL
     &    timein,timenow,timetot,timenew

      common /timeinfo/ timein,timenow,timetot,timenew
      save /timeinfo/

#endif /* _TIMEINFO_COM_ */

