
c This routine initializes individual bytes of dcore. Ideally,
c these would be set to garbage, e.g. 0xff, so any attempt to reference
c uninitialized data would result in an exception.

c INPUT
c    character*1 dcore(dim) : the array to initialize
c    integer     dim        : the byte length of the array

#include "iachar.h"

      subroutine init_dcore(dcore,dim)
      implicit none

c ARGUMENT LIST
      integer dim
      character*1 dcore(dim)

c INTERNAL VARIABLES
      integer i

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

      if (dim.gt.0) then
         do i = 1, dim
c            dcore(i) = achar(_IACHAR_NULL)
            dcore(i) = achar(255)
         end do
      end if

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

