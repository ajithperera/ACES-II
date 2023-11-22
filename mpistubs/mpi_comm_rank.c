
#include "mpi.h"

#include "f77_name.h"
#include "f_types.h"

int
MPI_Comm_rank
(MPI_Comm comm, int *rank)
{ *rank=0; return 0; }

void
F77_NAME(mpi_comm_rank,MPI_COMM_RANK)
(f_int * comm, f_int * rank, f_int * ierror)
{
    *rank   = 0;
    *ierror = MPI_SUCCESS;
    return;
}

