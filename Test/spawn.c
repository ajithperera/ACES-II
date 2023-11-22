
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

    printf("%s:%i: I am %i, child of %i.\n",
           szCPUName,iRank,
           (int)pid,(int)ppid);
    system("./spawn");

    MPI_Finalize();
    return 0;
}

