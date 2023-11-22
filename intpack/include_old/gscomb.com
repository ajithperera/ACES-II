
#ifndef _GSCOMB_COM_
#define _GSCOMB_COM_

      integer swap(4,8,2),lock(8,2),iprtgs(8,2),istop1,istop2,
     &    readmo,writmo,iuhfrhf,lugss

      common /gscomb/ swap,lock,iprtgs,istop1,istop2,
     &    readmo,writmo,iuhfrhf,lugss

#endif /* _GSCOMB_COM_ */

