
#include <unistd.h>
#include "mpi.h"

int main(int argc, char *argv[])
{
    char szCPUName[MPI_MAX_PROCESSOR_NAME];
    int nProcs, iRank, iCPUName; pid_t pid = getpid(), ppid = getppid();
    int iSend, iRecv;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD,&iRank);
    MPI_Get_processor_name(&szCPUName[0],&iCPUName);

    printf("%s:%i: I am %i, serving %i.\n",
           szCPUName,iRank,
           (int)pid,(int)ppid);
    printf("%i: last arg is '%s'\n",(int)pid,argv[argc-1]);

    iSend = iRank;
    MPI_Allreduce(&iSend,&iRecv,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    printf("%s:%i: %i has %i\n",szCPUName,iRank,(int)pid,iRecv);

    MPI_Finalize();
    return 0;
}

