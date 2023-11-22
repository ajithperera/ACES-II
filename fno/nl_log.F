      subroutine nl_log(key,def,var)
c This reads a logical value.  If the keyword has a value, the value must
c be one of the (case-insensitive) following:
c   yes, true,   1,  on,  y,  t
c   no,  false,  0,  off  n,  f
c If it does not have a value, it is set to .true.
c
c Variables which have default values may safely be passed in both
c as the default value and as the value read in.
#include <aces.h>
#include <callstack.com>
#include <io.par>
#include <namelist.com>
      character *(*) key
      logical var,default, def

      character*5 val,defval
      logical nl_key,present,streq
      integer strlen

      callstack_curr='NL_LOG'
      default = def
      val=' '
      present=nl_key(key,val)
      if (present) then
        if (strlen(val).gt.0) then
          if (streq(val,'yes',.true.) .or.
     &        streq(val,'true',.true.) .or.
     &        streq(val,'y',.true.) .or.
     &        streq(val,'t',.true.) .or.
     &        streq(val,'1',.true.) .or.
     &        streq(val,'on',.true.)) then
            var=.true.
          else if (streq(val,'no',.true.) .or.
     &          streq(val,'false',.true.) .or.
     &          streq(val,'n',.true.) .or.
     &          streq(val,'f',.true.) .or.
     &          streq(val,'0',.true.) .or.
     &          streq(val,'off',.true.)) then
            var=.false.
          else
            write(stdout,'(a)')
     &          '@NL_LOG-F: invalid value for ',key
            call errex
          end if
        else
          var=.true.
        endif
      else
        var=default
      endif

      if (prt_nl .and. .not.streq(key,'print_nl',.true.)) then
        val='false'
        if (var) val=' true'
        defval='false'
        if (default) defval=' true'
        if (val.eq.defval) then
          write(stdout,900) key,'logical',val
  900     format(a20,2x,a10,2x,a20)
        else
          write(stdout,910) key,'logical',val,defval
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
