#ifndef _SYMOFF_COM_
#define _SYMOFF_COM_
c symoff.com : begin
      integer  Ioff_oo(8,2),Ioff_vv(8,2)
      common /symoff/ Ioff_oo,Ioff_vv
c symoff.com : end
#endif /* _SYMOFF_COM_ */
