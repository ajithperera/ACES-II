#ifndef _MPWF_H_
#define _MPWF_H_

#if !defined(_HAVE_MPI) || defined(_64BIT_MPI)
#define F_MPI_INT integer
#define F_MPI_LOG logical
#else
#define F_MPI_INT integer*4
#define F_MPI_LOG logical*4
#endif

#endif /* _MPWF_H_ */
