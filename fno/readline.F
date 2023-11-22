      character *(*) function readline(line,unit,eof)

c This reads the next non-comment non-blank line.  If none is found,
c eof is set to true.  Otherwise, the line is returned.  A comment is
c defined as a line starting with a "#".  linenumber is the line number
c in the file.  It should be initalized to zero before this is called
c for the first time.

#include <aces.h>
#include <callstack.com>

      integer line,unit
      logical eof

      integer strlen

      callstack_curr='READLINE'

      eof=.false.
   10 read(unit,'(a)',end=900) readline
      line=line+1
      if (strlen(readline).eq.0 .or.  readline(1:1).eq.'#') goto 10
      return

  900 eof=.true.
      return
      end

c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End:
