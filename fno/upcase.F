      subroutine upcase(str)
#include <aces.h>
      character*(*) str
      integer l,i
c      callstack_curr='UPCASE'
      l=len(str)
      do 10 i=1,l
        if (str(i:i).ge.'a' .and. str(i:i).le.'z')
     &      str(i:i)=char(ichar(str(i:i))+ichar('A')-ichar('a'))
   10 continue
      return
      end

c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End:
