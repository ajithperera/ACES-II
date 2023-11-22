#ifndef _XCC_COM_
#define _XCC_COM_
c xcc.com : begin
      double precision overlap
      integer dcoresize, freecore
      integer ndx_t1(8,2), ndx_t2(8,3)
      common /xcc_com/ overlap, dcoresize, freecore, ndx_t1, ndx_t2
#include "syminf.com"
#include "sympop.com"
c xcc.com : end
#endif /* _XCC_COM_ */
