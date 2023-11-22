
#ifndef _ACES_H_
#define _ACES_H_

c**************************************************************************
c The following comments may be ignored when appearing in a .f file.
c They are only meaningful in a .F file.
c**************************************************************************

c  Macros beginning with M_ are machine dependant macros.
c  Macros beginning with F_ are general fortran macros.
c  Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"

#define M_REAL double precision
#define M_IMPLICITNONE

c Different fortran data types (used in parsing routines).
#define F_INTEGER   0
#define F_REAL      1
#define F_COMPLEX   2
#define F_LOGICAL   3
#define F_CHARACTER 4

#include <blas.h>
#include <matx.h>

c**************************************************************************
c End of special comment section
c**************************************************************************

#ifdef M_IMPLICITNONE
      implicit none
#else
      implicit logical (a-z)
#endif

#include <aces.par>

#endif /* _ACES_H_ */

