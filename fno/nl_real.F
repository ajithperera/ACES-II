      subroutine nl_real(key,def,var)
c This reads an integer value.  If it isn't in the namelist, it defaults
c to the value passed in.
c
c Variables which have default values may safely be passed in both
c as the default value and as the value read in.
#include <aces.h>
#include <callstack.com>
#include <io.par>
#include <namelist.com>
      character *(*) key
      M_REAL
     &    def,var

      character *(nllinelen) val
      logical nl_key,present
      integer i
      M_REAL
     &    c_atof,default

      callstack_curr='NL_REAL'
      default=def
      err=0
      val=' '
      present=nl_key(key,val)
      if (.not.present) then
        var=default
      else
        i=len(val)
        do while (i.gt.0.and.val(i:i).eq.' ')
           i = i-1
        end do
        if (i.eq.0) then
          write(stdout,'(a)') '@NL_REAL-F: invalid real for ',key
          call errex
        endif
        val(i+1:i+1) = achar(0)
        var=c_atof(val)
      endif

      if (prt_nl) then
        if (var.eq.default) then
          write(stdout,900) key,'real',var
  900     format(a20,2x,a10,2x,f20.8)
        else
          write(stdout,910) key,'real',var,default
  910     format(a20,2x,a10,2x,f20.8,2x,f20.8)
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
