      function ncindex(str,char,ind)
#include <aces.h>
      integer ncindex
      character *(*) str
      character*1 char
      integer ind,strlen,i,slen
c      callstack_curr='NCINDEX'
      ncindex=0
      if (ind.lt.0) return
      slen=strlen(str)
      i=ind+1
 10   if (i.gt.slen) then
        ncindex=0
        return
      end if
      if (str(i:i).eq.char) then
        ncindex=i
        return
      end if
      i=i+1
      goto 10
      end

c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End:
