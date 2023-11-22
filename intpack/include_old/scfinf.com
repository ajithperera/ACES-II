
#ifndef _SCFINF_COM_
#define _SCFINF_COM_

c This stores some information about the current SCF iteration

      M_REAL
     &    dampavg,dampe0

      common /scfinfd/ dampavg,dampe0

#endif /* _SCFINF_COM_ */

