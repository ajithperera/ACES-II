
#include "mpi.h"

#include "f77_name.h"
#include "f_types.h"

int
MPI_Barrier
(MPI_Comm comm)
{ return 0; }

void
F77_NAME(mpi_barrier,MPI_BARRIER)
(f_int * comm, f_int * ierr)
{ *ierr=MPI_SUCCESS; return; }

