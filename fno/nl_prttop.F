      subroutine nl_prttop(nltitle)
c This prints the top of a table of namelist values.
#include <aces.h>
#include <callstack.com>
#include <io.par>
      character*(*) nltitle
      callstack_curr='NL_PRTTOP'
      write(stdout,990)
      write(stdout,900)
      write(stdout,910) nltitle
      write(stdout,920) 'KEYWORD','TYPE','VALUE','DEFAULT'
      write(stdout,930)
      write(stdout,920) 'print_nl','logical','true'
  900 format(76('='))
  910 format('NAMELIST: ',a)
  920 format(a20,2x,a10,2x,a20,2x,a20)
  930 format(20('-'),2x,10('-'),2x,20('-'),2x,20('-'),2x)
  990 format(' ')
      return
      end
 
c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End:
