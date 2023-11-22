      function fcindex(str,char)
#include <aces.h>
      integer fcindex
      character *(*) str
      character*1 char
      integer i,strlen,slen
c      callstack_curr='FCINDEX'
      i=1
      slen=strlen(str)
 10   if (i.gt.slen) then
        fcindex=0
        return
      end if
      if (str(i:i).eq.char) then
        fcindex=i
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
