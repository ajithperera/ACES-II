      subroutine nl_str(key,def,var)
c This reads a string value.  If it isn't in the namelist, it defaults
c to the value passed in.
c
c Variables which have default values may safely be passed in both
c as the default value and as the value read in.
#include <aces.h>
#include <callstack.com>
#include <io.par>
#include <namelist.com>
      character *(*) key,def,var

      character*(nllinelen) default,str1,str2
      logical nl_key,present
      integer strlen

      callstack_curr='NL_STR'
      default=def
      present=nl_key(key,var)
      if (.not.present) var=default

      if (prt_nl) then
        if (strlen(var).le.20) then
          str1=var
        else
          str1=var(1:18)//' $'
        endif
        if (strlen(default).le.20) then
          str2=default
        else
          str2=default(1:18)//' $'
        endif
        if (var(1:strlen(var)).eq.default(1:strlen(default))) then
          write(stdout,900) key,'string',str1
  900     format(a20,2x,a10,2x,a20)
        else
          write(stdout,910) key,'string',str1,str2
  910     format(a20,2x,a10,2x,a20,2x,a20)
        endif
      endif
      return
      end

c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End:
