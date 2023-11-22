
#include <string.h>
#include <stdio.h>

#include "mpi.h"

#include "f77_name.h"
#include "f_types.h"

int
MPI_Allgatherv
(void *sendbuf, int  sendcount,  MPI_Datatype sendtype,
 void *recvbuf, int *recvcounts, int *displs,
                                 MPI_Datatype recvtype, MPI_Comm comm)
{
    if (recvbuf != sendbuf)
    {
        size_t bytes;

        if (sendtype != recvtype)
        {
            printf("\nMPI_Allgatherv: moving different data types is not "
                     "supported\n");
            abort();
        }

        if (sendcount != *recvcounts)
        {
            printf("\nMPI_Allgatherv: sending and receiving different "
                     "amounts of data\n"
                     "                are not supported\n");
            abort();
        }

        if (*displs != 0)
        {
            printf("\nMPI_Allgatherv: receiving data beyond the receive "
                     "buffer address\n"
                     "                is not supported\n");
            abort();
        }

        switch (sendtype)
        {
            case MPI_DOUBLE_PRECISION:
            { bytes = sendcount*sizeof(double); break; }
            case MPI_INTEGER:
            { bytes = sendcount*sizeof(f_int); break; }
            default:
            { printf("\nMPI_Allgatherv: must add type=%i\n",(int)sendtype);
              abort(); break; }
        }

        memmove(recvbuf,sendbuf,bytes);

    }
    return 0;
}

void
F77_NAME(mpi_allgatherv,MPI_ALLGATHERV)
(void * sendbuf, f_int * sendcount,                  f_int * sendtype,
 void * recvbuf, f_int * recvcounts, f_int * displs, f_int * recvtype,
 f_int * comm, f_int * ierror)
{
    int m_sendcount  = *sendcount;
    int m_recvcounts = *recvcounts;
    int m_displs     = *displs;
    MPI_Datatype m_sendtype = *sendtype;
    MPI_Datatype m_recvtype = *recvtype;
    MPI_Comm m_comm = *comm;

    MPI_Allgatherv(sendbuf, m_sendcount,           m_sendtype,
                   recvbuf,&m_recvcounts,&m_displs,m_recvtype,m_comm);
    *ierror = MPI_SUCCESS;
    return;
}

