
#include "mpi.h"

#include "f77_name.h"
#include "f_types.h"

int
MPI_Comm_split
(MPI_Comm oldcomm, int color, int key, MPI_Comm *newcomm)
{ return 0; }

void
F77_NAME(mpi_comm_split,MPI_COMM_SPLIT)
(f_int * oldcomm, f_int * color, f_int * key, f_int * newcomm, f_int * ierror)
{ *ierror=MPI_SUCCESS; return; }

