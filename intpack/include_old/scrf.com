
#ifndef _SCRF_COM_
#define _SCRF_COM_

      M_REAL
     &    ffact,gfact,eborn

      common /scrf/ ffact,gfact,eborn

#endif /* _SCRF_COM_ */

