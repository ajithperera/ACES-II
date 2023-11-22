
#include <stdio.h> /* for sprintf */

#include "mpi.h"

int
MPI_Get_processor_name
(char *name, int *resultlen)
{ sprintf(name,"localhost"); *resultlen=9; return 0; }

