
#include <stdlib.h>

#include "mpi.h"

#include "f77_name.h"
#include "f_types.h"

int
MPI_Abort
(MPI_Comm comm, int errorcode)
{ abort(); }

void
F77_NAME(mpi_abort,MPI_ABORT)
(f_int * comm, f_int * errorcode, f_int * ierror)
{ abort(); }

