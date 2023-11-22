
#ifndef _GSCOMC_COM_
#define _GSCOMC_COM_

      logical rohfmo
      common /gscomc/ rohfmo

#endif /* _GSCOMC_COM_ */

