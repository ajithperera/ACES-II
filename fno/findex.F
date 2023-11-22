      function findex(str,substr)
#include <aces.h>
      integer findex
      character*(*) str,substr
      integer strlen,sublen,strmin,ind,ncindex,slen
c      callstack_curr='FINDEX'
      slen=strlen(str)
      sublen=strlen(substr)
      if (sublen.gt.slen) then
        findex=0
        return
      end if
      strmin=slen-sublen+1
      ind=0
 10   findex=ncindex(str,substr(1:1),ind)
      if (findex.eq.0) return
      if (findex.gt.strmin) then
        findex=0
        return
      end if
      if (str(findex:findex+sublen-1).eq.substr(1:sublen)) return
      ind=findex
      goto 10
      end

c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End:
