
#include "mpi.h"

#include "f77_name.h"
#include "f_types.h"

int
MPI_Init
(int *argc, char ***argv)
{ return 0; }

void
F77_NAME(mpi_init,MPI_INIT)
(f_int * ierror)
{ *ierror=MPI_SUCCESS; return; }

