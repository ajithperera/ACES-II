
#include "mpi.h"

#include "f77_name.h"
#include "f_types.h"

int
MPI_Finalize
()
{ return 0; }

void
F77_NAME(mpi_finalize,MPI_FINALIZE)
(f_int * ierror)
{ *ierror=MPI_SUCCESS; return; }

