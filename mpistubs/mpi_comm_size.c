
#include "mpi.h"

#include "f77_name.h"
#include "f_types.h"

int
MPI_Comm_size
(MPI_Comm comm, int *size)
{ *size=1; return 0; }

void
F77_NAME(mpi_comm_size,MPI_COMM_SIZE)
(f_int * comm, f_int * size, f_int * ierror)
{
    *size   = 1;
    *ierror = MPI_SUCCESS;
    return;
}

