
#include <unistd.h>
#include "mpi.h"

int main(int argc, char *argv[])
{
    char szCPUName[MPI_MAX_PROCESSOR_NAME];
    int nProcs, iRank, iCPUName; pid_t pid = getpid(), ppid = getppid();
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD,&iRank);
    MPI_Get_processor_name(&szCPUName[0],&iCPUName);

    printf("%s:%i: I am client %i.\n",
           szCPUName,iRank,
           (int)pid,(int)ppid);
#ifdef _AIX
    if (iRank < 2)
         system("/ufl/qtp/rjb/ay/scr/mpi_cli_ser/power/server true");
    else
         system("/ufl/qtp/rjb/ay/scr/mpi_cli_ser/power/server false");
    system("/ufl/qtp/rjb/ay/scr/mpi_cli_ser/power/server all");
#else
    if (iRank) system("/usr/local/bin/mpirun -np 4 ./server 1");
#endif

    MPI_Finalize();
    return 0;
}

