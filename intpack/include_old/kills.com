
#ifndef _KILLS_COM_
#define _KILLS_COM_

      integer killflg(3)
      common /kills/ killflg

#endif /* _KILLS_COM_ */

