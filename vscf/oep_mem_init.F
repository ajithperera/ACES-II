      subroutine oep_mem_init(maxmem, irppmem, imemget)
c 
      implicit double precision (a-h, o-z)
c
#include <machsp.com>
#include <flags.com>
#include <flags2.com>
      maxmem  = iflags (36) - irppmem - iflags(44)*iflags(37)
      maxmem  = maxmem      + mod (maxmem, 2)
      imemget = maxmem      + irppmem + iflags(44)*iflags(37)
      imemget = imemget     + mod (imemget, 2)
     
      Return
      End

