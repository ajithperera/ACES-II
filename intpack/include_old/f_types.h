
#ifndef _F_TYPES_H_
#define _F_TYPES_H_

#ifdef __fortran

#ifndef F_INT
#define F_INT INTEGER
#endif

#ifndef F_ADR
#define F_ADR INTEGER
#endif

#else /* C */

#ifdef F_64BIT
typedef long f_int;
#else
typedef int  f_int;
#endif

#ifdef F_64BIT
typedef long    f_adr;
#else
#ifdef F_ADR /* assume F_INT is INT*4 and F_ADR is INT*8 */
#include <inttypes.h>
typedef int64_t f_adr;
#else
typedef int     f_adr;
#endif
#endif

#endif

#endif /* _F_TYPES_H_ */

