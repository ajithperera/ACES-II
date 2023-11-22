#ifndef _ACES_H_
#define _ACES_H_

#undef M_REAL
#undef M_SINGLE
#undef M_DOUBLE

#undef M_IMPLICIT
#undef M_TRACEBACK

#define F_INTEGER   0 /* ACES type handles */
#define F_REAL      1
#define F_COMPLEX   2
#define F_LOGICAL   3
#define F_CHARACTER 4 /* M_ defines are machine-/compiler- dependent types. */

#ifdef __fortran

#   ifdef __fortran77

c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler

#ifdef _UNICOS
#   define M_REAL real
#   define M_SINGLE
#   define M_IMPLICIT implicit none
#   define M_TRACEBACK trbk
#else
#   define M_REAL double precision
#   define M_DOUBLE
#   define M_IMPLICIT implicit none
#endif

#include <blas.h>
#include <matx.h>
c#include <aces.par>

#   else /* __fortran77 */

! INSERT FORTRAN9x CODE HERE

#   endif /* __fortran77 */

#else /* __fortran */

#define M_REAL double

#endif /* __fortran */

#endif /* _ACES_H_ */
