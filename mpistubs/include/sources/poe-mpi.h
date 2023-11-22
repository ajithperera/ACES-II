/* IBM_PROLOG_BEGIN_TAG                                                   */
/* This is an automatically generated prolog.                             */
/*                                                                        */
/*                                                                        */
/*                                                                        */
/* Licensed Materials - Property of IBM                                   */
/*                                                                        */
/* (C) COPYRIGHT International Business Machines Corp. 1994,2000          */
/* All Rights Reserved                                                    */
/*                                                                        */
/* US Government Users Restricted Rights - Use, duplication or            */
/* disclosure restricted by GSA ADP Schedule Contract with IBM Corp.      */
/*                                                                        */
/* IBM_PROLOG_END_TAG                                                     */
#ifndef _H_MPI
#define _H_MPI
/****************************************************************************
@(#) 1.31.1.18 src/ppe/poe/include/mpi.h, ppe.poe.mpi, ppe_rmoh, rmoht5du 00/06/12 10:08:52

 Name: mpi.h

 Description:  C/C++ header file describing the Message Passing Interface 
               implemented by the Parallel Environment                                                         
                                                                         
 Notes:  "#ifdef _THREAD_SAFE" is used when making MPI-2 definitions etc. 
         that will be available in the MPI threaded library but not the 
         signal library. 
         On AIX "_THREAD_SAFE" is set when the mpcc_r or mpCC_r script is 
         used to compile an application.  
         OS/390 supports only one libary version:
         V2R4 - V2R6: the signal libray (i.e. "_THREAD_SAFE must not be used)
         V2R7 and following releases: the threaded library 
             (i.e. "_THREAD_SAFE"  must always be used -- mpcc and mpCC set 
             "_THREAD_SAFE" by default.).

****************************************************************************/

#define MPI_VERSION 1
#define MPI_SUBVERSION 2

#define MPI_SUCCESS 0
#ifdef _THREAD_SAFE
enum {MPI_ERR_BUFFER=50,MPI_ERR_COUNT,MPI_ERR_TYPE,MPI_ERR_TAG,
      MPI_ERR_COMM,MPI_ERR_RANK,MPI_ERR_REQUEST,MPI_ERR_ROOT,MPI_ERR_GROUP,
      MPI_ERR_OP,MPI_ERR_TOPOLOGY,MPI_ERR_DIMS,MPI_ERR_ARG,MPI_ERR_UNKNOWN,
      MPI_ERR_TRUNCATE,MPI_ERR_OTHER,MPI_ERR_INTERN,MPI_ERR_IN_STATUS,
      MPI_ERR_PENDING,MPI_ERR_INFO_KEY,MPI_ERR_INFO_VALUE,MPI_ERR_INFO_NOKEY,
      MPI_ERR_INFO,MPI_ERR_FILE,MPI_ERR_NOT_SAME,MPI_ERR_AMODE,
      MPI_ERR_UNSUPPORTED_DATAREP,MPI_ERR_UNSUPPORTED_OPERATION,
      MPI_ERR_NO_SUCH_FILE,MPI_ERR_FILE_EXISTS,MPI_ERR_BAD_FILE,MPI_ERR_ACCESS,
      MPI_ERR_NO_SPACE,MPI_ERR_QUOTA,MPI_ERR_READ_ONLY,MPI_ERR_FILE_IN_USE,
      MPI_ERR_DUP_DATAREP,MPI_ERR_CONVERSION,MPI_ERR_IO,MPI_ERR_WIN,
      MPI_ERR_BASE,MPI_ERR_SIZE,MPI_ERR_DISP,MPI_ERR_LOCKTYPE,MPI_ERR_ASSERT,
      MPI_ERR_RMA_CONFLICT,MPI_ERR_RMA_SYNC,MPI_ERR_NO_MEM};
#define MPI_ERR_LASTCODE 500
#else
enum {MPI_ERR_BUFFER=50,MPI_ERR_COUNT,MPI_ERR_TYPE,MPI_ERR_TAG,
      MPI_ERR_COMM,MPI_ERR_RANK,MPI_ERR_REQUEST,MPI_ERR_ROOT,MPI_ERR_GROUP,
      MPI_ERR_OP,MPI_ERR_TOPOLOGY,MPI_ERR_DIMS,MPI_ERR_ARG,MPI_ERR_UNKNOWN,
      MPI_ERR_TRUNCATE,MPI_ERR_OTHER,MPI_ERR_INTERN,MPI_ERR_IN_STATUS,
      MPI_ERR_PENDING,MPI_ERR_NOT_SAME=74};
#define MPI_ERR_LASTCODE 250
#endif
#define MPI_PENDING MPI_ERR_PENDING

#define MPI_BOTTOM ((void*) 0)
#define MPI_PROC_NULL -3
#define MPI_ANY_SOURCE -1
#define MPI_ANY_TAG -1
#define MPI_UNDEFINED -1

#define MPI_SOURCE source
#define MPI_TAG tag
#define MPI_ERROR error

enum { MPI_ERRORS_ARE_FATAL, MPI_ERRORS_RETURN, MPE_ERRORS_WARN };

#define MPI_MAX_PROCESSOR_NAME 256
#define MPI_MAX_ERROR_STRING   128
#define MPI_BSEND_OVERHEAD      23
#ifdef _THREAD_SAFE
#define MPI_MAX_FILE_NAME      1023
#define MPI_MAX_DATAREP_STRING 255
#define MPI_MAX_INFO_KEY       127
#define MPI_MAX_INFO_VAL       1023
#endif

#ifdef _THREAD_SAFE

/* combiner values used for datatype decoding functions */
#define MPI_COMBINER_NAMED            0  /* a named predefined datatype    */
#define MPI_COMBINER_DUP              1  /* MPI_TYPE_DUP                   */
#define MPI_COMBINER_CONTIGUOUS       2  /* MPI_TYPE_CONTIGUOUS            */
#define MPI_COMBINER_VECTOR           3  /* MPI_TYPE_VECTOR                */
#define MPI_COMBINER_HVECTOR_INTEGER  4  /* MPI_TYPE_HVECTOR from Fortran  */
#define MPI_COMBINER_HVECTOR          5  /* MPI_TYPE_HVECTOR from C and    */
                                         /*   C++ and in some case Fortran */
                                         /*   or MPI_TYPE_CREATE_HVECTOR   */
#define MPI_COMBINER_INDEXED          6  /* MPI_TYPE_INDEXED               */
#define MPI_COMBINER_HINDEXED_INTEGER 7  /* MPI_TYPE_HINDEXED from Fortran */
#define MPI_COMBINER_HINDEXED         8  /* MPI_TYPE_HINDEXED from C and   */
                                         /*   C++ and in some case Fortran */
                                         /*   or MPI_TYPE_CREATE_HINDEXED  */
#define MPI_COMBINER_INDEXED_BLOCK    9  /* MPI_TYPE_CREATE_INDEXED_BLOCK  */
#define MPI_COMBINER_STRUCT_INTEGER   10 /* MPI_TYPE_STRUCT from Fortran   */
#define MPI_COMBINER_STRUCT           11 /* MPI_TYPE_STRUCT from C and     */
                                         /*   C++ and in some case Fortran */
                                         /*   or MPI_TYPE_CREATE_STRUCT    */
#define MPI_COMBINER_SUBARRAY         12 /* MPI_TYPE_CREATE_SUBARRAY       */
#define MPI_COMBINER_DARRAY           13 /* MPI_TYPE_CREATE_DARRAY         */
#define MPI_COMBINER_F90_REAL         14 /* MPI_TYPE_CREATE_F90_REAL       */
#define MPI_COMBINER_F90_COMPLEX      15 /* MPI_TYPE_CREATE_F90_COMPLEX    */
#define MPI_COMBINER_F90_INTEGER      16 /* MPI_TYPE_CREATE_F90_INTEGER    */
#define MPI_COMBINER_RESIZED          17 /* MPI_TYPE_CREATE_RESIZED        */

#endif /* _THREAD_SAFE */

/* MPI special purpose datatypes  */
#define MPI_LB                 0
#define MPI_UB                 1
#define MPI_BYTE               2
#define MPI_PACKED             3
/* MPI Predefined Datatypes for C language bindings */
#define MPI_CHAR               4
#define MPI_UNSIGNED_CHAR      5
#define MPI_SIGNED_CHAR        6
#define MPI_SHORT              7
#define MPI_INT                8
#define MPI_LONG               9
#define MPI_UNSIGNED_SHORT     10
#define MPI_UNSIGNED           11
#define MPI_UNSIGNED_LONG      12
#define MPI_FLOAT              13
#define MPI_DOUBLE             14
#define MPI_LONG_DOUBLE        15

#define MPI_LONG_LONG_INT      39
#define MPI_LONG_LONG          MPI_LONG_LONG_INT
#define MPI_UNSIGNED_LONG_LONG 40
#define MPI_WCHAR              41

/* MPI Predefined Datatypes for Fortran language bindings */
#define MPI_INTEGER1           16
#define MPI_INTEGER2           17
#define MPI_INTEGER4           18
#define MPI_INTEGER            MPI_INTEGER4
#define MPI_REAL4              19
#define MPI_REAL               MPI_REAL4
#define MPI_REAL8              20
#define MPI_DOUBLE_PRECISION   MPI_REAL8
#define MPI_REAL16             21
#define MPI_COMPLEX8           22
#define MPI_COMPLEX            MPI_COMPLEX8
#define MPI_COMPLEX16          23
#define MPI_DOUBLE_COMPLEX     MPI_COMPLEX16
#define MPI_COMPLEX32          24
#define MPI_LOGICAL1           25
#define MPI_LOGICAL2           26
#define MPI_LOGICAL4           27
#define MPI_LOGICAL            MPI_LOGICAL4
#define MPI_CHARACTER          28

#define MPI_INTEGER8           42
#define MPI_LOGICAL8           43

/* MPI Predefined datatypes for reduction functions */
           /* C  Reduction Types  */
#define MPI_FLOAT_INT          29  /* {MPI_FLOAT, MPI_INT} */
#define MPI_DOUBLE_INT         30  /* {MPI_DOUBLE, MPI_INT} */
#define MPI_LONG_INT           31  /* {MPI_LONG, MPI_INT} */
#define MPI_2INT               32  /* {MPI_INT, MPI_INT} */
#define MPI_SHORT_INT          33  /* {MPI_SHORT, MPI_INT */
#define MPI_LONG_DOUBLE_INT    34  /* {MPI_LONG_DOUBLE, MPI_INT} */
           /* Fortran Reduction Types */
#define MPI_2REAL              35  /* {MPI_REAL, MPI_REAL} */
#define MPI_2DOUBLE_PRECISION  36  /* {MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION} */
#define MPI_2INTEGER           37  /* {MPI_INTEGER, MPI_INTEGER} */
#define MPI_2COMPLEX           38  /* {MPI_COMPLEX, MPI_COMPLEX} */

enum { MPI_COMM_WORLD, MPI_COMM_SELF };

/* the order of these is important! */
enum { MPI_IDENT, MPI_CONGRUENT, MPI_SIMILAR, MPI_UNEQUAL };

enum { MPI_TAG_UB, MPI_IO, MPI_HOST, MPI_WTIME_IS_GLOBAL, MPI_WIN_BASE, MPI_WIN_SIZE,
       MPI_WIN_DISP_UNIT };

enum { MPI_MAX, MPI_MIN, MPI_SUM, MPI_PROD, MPI_MAXLOC, MPI_MINLOC,
       MPI_BAND, MPI_BOR, MPI_BXOR, MPI_LAND, MPI_LOR, MPI_LXOR, MPI_REPLACE, MAX_OP };

#define MPI_HANDLE_NULL     -1
#define MPI_GROUP_NULL      MPI_HANDLE_NULL
#define MPI_COMM_NULL       MPI_HANDLE_NULL
#define MPI_DATATYPE_NULL   MPI_HANDLE_NULL
#define MPI_REQUEST_NULL    MPI_HANDLE_NULL
#define MPI_OP_NULL         MPI_HANDLE_NULL
#define MPI_ERRHANDLER_NULL MPI_HANDLE_NULL
#define MPI_KEYVAL_INVALID  MPI_HANDLE_NULL /* this is missing from the spec */
#ifdef _THREAD_SAFE
#define MPI_INFO_NULL       MPI_HANDLE_NULL
#define MPI_FILE_NULL       MPI_HANDLE_NULL
#define MPI_WIN_NULL        MPI_HANDLE_NULL
#endif

#define MPI_GROUP_EMPTY 0

#ifdef _THREAD_SAFE
#define MPI_NON_ATOMIC              0
#define MPI_ATOMIC                  1

#define MPI_DISPLACEMENT_CURRENT    -1LL

#define MPI_DISTRIBUTE_NONE         0
#define MPI_DISTRIBUTE_BLOCK        1
#define MPI_DISTRIBUTE_CYCLIC       2
#define MPI_DISTRIBUTE_DFLT_DARG    0

#define MPI_ORDER_C                 1
#define MPI_ORDER_FORTRAN           2

#define MPI_SEEK_SET                0
#define MPI_SEEK_CUR                1
#define MPI_SEEK_END                2

#define MPI_MODE_RDONLY             0x000001
#define MPI_MODE_WRONLY             0x000002
#define MPI_MODE_RDWR               0x000004
#define MPI_MODE_CREATE             0x000008
#define MPI_MODE_APPEND             0x000010
#define MPI_MODE_EXCL               0x000020
#define MPI_MODE_DELETE_ON_CLOSE    0x000040
#define MPI_MODE_UNIQUE_OPEN        0x000080
#define MPI_MODE_SEQUENTIAL         0x000100

#define MPI_LOCK_EXCLUSIVE          0
#define MPI_LOCK_SHARED             1

#define MPI_MODE_NOCHECK            0x000200
#define MPI_MODE_NOSTORE            0x000400
#define MPI_MODE_NOPUT              0x000800
#define MPI_MODE_NOPRECEDE          0x001000
#define MPI_MODE_NOSUCCEED          0x002000
#endif

enum { MPI_GRAPH, MPI_CART };

typedef long MPI_Aint;
typedef int MPI_Handle;
typedef MPI_Handle MPI_Group;
typedef MPI_Handle MPI_Comm;
typedef MPI_Handle MPI_Datatype;
typedef MPI_Handle MPI_Request;
typedef MPI_Handle MPI_Op;
typedef MPI_Handle MPI_Errhandler;        /* this is missing from the spec */

typedef struct {
   int source;
   int tag;
   int error;
   MPI_Aint val1;
   int val2;
   int val3;
   int val4;
   int val5;
} MPI_Status;

#ifdef _THREAD_SAFE
typedef MPI_Handle MPI_Info;
typedef MPI_Handle MPI_File;
typedef MPI_Handle MPI_Win;
typedef long long int MPI_Offset;
#endif

#if defined(__cplusplus)
   extern "C" {
#endif

typedef int MPI_Copy_function(MPI_Comm,int,void *,void *,void *,int *);
typedef int MPI_Delete_function(MPI_Comm,int,void *,void *);

typedef void MPI_Handler_function(MPI_Comm *,int *,...);

typedef void MPI_User_function(void *invec,void *inoutvec,int *len,
                               MPI_Datatype *datatype);
#ifdef _THREAD_SAFE
typedef int MPI_Comm_copy_attr_function(MPI_Comm oldcomm, int comm_keyval, void *extra_state,
              void *attribute_val_in, void *attribute_val_out, int *flag);
typedef int MPI_Comm_delete_attr_function(MPI_Comm comm, int comm_keyval, void *attribute_val,
              void *extra_state);
typedef int MPI_Win_copy_attr_function(MPI_Win oldwin, int win_keyval, void *extra_state,
              void *attribute_val_in, void *attribute_val_out, int *flag);
typedef int MPI_Win_delete_attr_function(MPI_Win win, int win_keyval, void *attribute_val,
              void *extra_state);
typedef int MPI_Type_copy_attr_function(MPI_Datatype oldwin, int win_keyval, void *extra_state,
              void *attribute_val_in, void *attribute_val_out, int *flag);
typedef int MPI_Type_delete_attr_function(MPI_Datatype win, int win_keyval, void *attribute_val,
              void *extra_state);

typedef void MPI_File_errhandler_fn(MPI_File *,int *,...);
typedef void MPI_Win_errhandler_fn(MPI_Win *,int *,...);
typedef void MPI_Comm_errhandler_fn(MPI_Comm *,int *,...);

typedef int MPI_Datarep_extent_function(MPI_Datatype,MPI_Aint *,void *);
typedef int MPI_Datarep_conversion_function(void *,MPI_Datatype,int,void *,
                                            MPI_Offset,void *);
#define MPI_CONVERSION_FN_NULL  0
#endif

#if defined(__cplusplus)
   } 
#endif

#define MPI_NULL_COPY_FN         0
#define MPI_NULL_DELETE_FN       0
#define MPI_COMM_NULL_COPY_FN    0
#define MPI_COMM_NULL_DELETE_FN  0
#define MPI_TYPE_NULL_COPY_FN    0
#define MPI_TYPE_NULL_DELETE_FN  0
#define MPI_WIN_NULL_COPY_FN     0
#define MPI_WIN_NULL_DELETE_FN   0
MPI_Copy_function       _mpi_dup_fn;
#define MPI_DUP_FN      _mpi_dup_fn
#define MPI_COMM_DUP_FN _mpi_dup_fn
#define MPI_TYPE_DUP_FN _mpi_dup_fn
#define MPI_WIN_DUP_FN  _mpi_dup_fn

#if defined(__cplusplus)
  extern "C" {
#endif


   /* C Bindings for Point-to-Point Communication */

int MPI_Send(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
             MPI_Comm comm);
int MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source, int tag,
             MPI_Comm comm, MPI_Status *status);
int MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int *count);
int MPI_Bsend(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
             MPI_Comm comm);
int MPI_Ssend(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
             MPI_Comm comm);
int MPI_Rsend(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
             MPI_Comm comm);
int MPI_Buffer_attach( void* buffer, int size);
int MPI_Buffer_detach( void* buffer, int* size);
int MPI_Isend(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request);
int MPI_Ibsend(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request);
int MPI_Issend(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request);
int MPI_Irsend(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request);
int MPI_Irecv(void* buf, int count, MPI_Datatype datatype, int source, int tag,
              MPI_Comm comm, MPI_Request *request);
int MPI_Wait(MPI_Request *request, MPI_Status *status);
int MPI_Test(MPI_Request *request, int *flag, MPI_Status *status);
int MPI_Request_free(MPI_Request *request);
int MPI_Waitany(int count, MPI_Request *array_of_requests, int *index,
              MPI_Status *status);
int MPI_Testany(int count, MPI_Request *array_of_requests, int *index, int *flag,
              MPI_Status *status);
int MPI_Waitall(int count, MPI_Request *array_of_requests,
              MPI_Status *array_of_statuses);
int MPI_Testall(int count, MPI_Request *array_of_requests, int *flag,
              MPI_Status *array_of_statuses);
int MPI_Waitsome(int incount, MPI_Request *array_of_requests, int *outcount,
              int *array_of_indices, MPI_Status *array_of_statuses);
int MPI_Testsome(int incount, MPI_Request *array_of_requests, int *outcount,
              int *array_of_indices, MPI_Status *array_of_statuses);
int MPI_Iprobe(int source, int tag, MPI_Comm comm, int *flag, MPI_Status *status);
int MPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status *status);
int MPI_Cancel(MPI_Request *request);
int MPI_Test_cancelled(MPI_Status *status, int *flag);
int MPI_Send_init(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request);
int MPI_Bsend_init(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request);
int MPI_Ssend_init(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request);
int MPI_Rsend_init(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request);
int MPI_Recv_init(void* buf, int count, MPI_Datatype datatype, int source, int tag,
              MPI_Comm comm, MPI_Request *request);
int MPI_Start(MPI_Request *request);
int MPI_Startall(int count, MPI_Request *array_of_requests);
int MPI_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest,
              int sendtag, void *recvbuf, int recvcount, MPI_Datatype recvtype,
              int source, int recvtag, MPI_Comm comm, MPI_Status *status);
int MPI_Sendrecv_replace(void* buf, int count, MPI_Datatype datatype, int dest,
              int sendtag, int source, int recvtag, MPI_Comm comm, MPI_Status *status);
int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype);
int MPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype,
              MPI_Datatype *newtype);
int MPI_Type_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype,
              MPI_Datatype *newtype);
int MPI_Type_indexed(int count, int *array_of_blocklengths, int *array_of_displacements,
              MPI_Datatype oldtype, MPI_Datatype *newtype);
int MPI_Type_hindexed(int count, int *array_of_blocklengths, MPI_Aint *array_of_displacements,
              MPI_Datatype oldtype, MPI_Datatype *newtype);
int MPI_Type_struct(int count, int *array_of_blocklengths, MPI_Aint *array_of_displacements,
              MPI_Datatype *array_of_types, MPI_Datatype *newtype);

#ifdef _THREAD_SAFE
int MPI_Type_get_envelope(MPI_Datatype datatype, int *num_integers, int *num_addresses,
              int *num_datatypes, int *combiner );
int MPI_Type_get_contents(MPI_Datatype datatype, int max_integers, int max_addresses,
              int max_datatypes, int array_of_integers[], MPI_Aint array_of_addresses[],
              MPI_Datatype array_of_datatypes[] );
int MPI_Type_create_keyval(MPI_Type_copy_attr_function *type_copy_attr_fn,
              MPI_Type_delete_attr_function *type_delete_attr_fn, int *type_keyval,
              void *extra_state);
int MPI_Type_delete_attr(MPI_Datatype type, int type_keyval);
int MPI_Type_free_keyval(int *type_keyval);
int MPI_Type_get_attr(MPI_Datatype type, int type_keyval, void *attribute_val, int *flag);
int MPI_Type_set_attr(MPI_Datatype type, int type_keyval, void *attribute_val);
int MPI_Type_dup(MPI_Datatype type, MPI_Datatype *newtype);
#endif

int MPI_Address(void* location, MPI_Aint *address);
int MPI_Type_extent(MPI_Datatype datatype, MPI_Aint *extent);
int MPI_Type_size(MPI_Datatype datatype, int *size);
int MPI_Type_lb(MPI_Datatype datatype, MPI_Aint *displacement);
int MPI_Type_ub(MPI_Datatype datatype, MPI_Aint *displacement);
int MPI_Type_commit(MPI_Datatype *datatype);
int MPI_Type_free(MPI_Datatype *datatype);
int MPI_Get_elements(MPI_Status *status, MPI_Datatype datatype, int *count);
int MPI_Pack(void* inbuf, int incount, MPI_Datatype datatype, void *outbuf, int outsize,
              int *position,  MPI_Comm comm);
int MPI_Unpack(void* inbuf, int insize, int *position, void *outbuf, int outcount,
              MPI_Datatype datatype, MPI_Comm comm);
int MPI_Pack_size(int incount, MPI_Datatype datatype, MPI_Comm comm, int *size);

    /*  C Bindings for Collective Communication   */

int MPI_Barrier(MPI_Comm comm );
int MPI_Bcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm );
int MPI_Gather(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
              int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) ;
int MPI_Gatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
              int *recvcounts, int *displs, MPI_Datatype recvtype, int root,
              MPI_Comm comm);
int MPI_Scatter(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
              int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
int MPI_Scatterv(void* sendbuf, int *sendcounts, int *displs, MPI_Datatype sendtype,
              void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
int MPI_Allgather(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
              int recvcount, MPI_Datatype recvtype, MPI_Comm comm);
int MPI_Allgatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
              int *recvcounts, int *displs, MPI_Datatype recvtype, MPI_Comm comm);
int MPI_Alltoall(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
              int recvcount, MPI_Datatype recvtype, MPI_Comm comm);
int MPI_Alltoallv(void* sendbuf, int *sendcounts, int *sdispls, MPI_Datatype sendtype,
              void* recvbuf, int *recvcounts, int *rdispls, MPI_Datatype recvtype,
              MPI_Comm comm);
int MPI_Reduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype,
              MPI_Op op, int root, MPI_Comm comm);
int MPI_Op_create(MPI_User_function *function, int commute, MPI_Op *op);
int MPI_Op_free( MPI_Op *op);
int MPI_Allreduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype,
              MPI_Op op, MPI_Comm comm);
int MPI_Reduce_scatter(void* sendbuf, void* recvbuf, int *recvcounts,
              MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int MPI_Scan(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype,
              MPI_Op op, MPI_Comm comm );

    /* C  Bindings for Groups, Contexts, and Communicators  */

int MPI_Group_size(MPI_Group group, int *size);
int MPI_Group_rank(MPI_Group group, int *rank);
int MPI_Group_translate_ranks (MPI_Group group1, int n, int *ranks1,
              MPI_Group group2, int *ranks2);
int MPI_Group_compare(MPI_Group group1,MPI_Group group2, int *result);
int MPI_Comm_group(MPI_Comm comm, MPI_Group *group);
int MPI_Group_union(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup);
int MPI_Group_intersection(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup);
int MPI_Group_difference(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup);
int MPI_Group_incl(MPI_Group group, int n, int *ranks, MPI_Group *newgroup);
int MPI_Group_excl(MPI_Group group, int n, int *ranks, MPI_Group *newgroup);
int MPI_Group_range_incl(MPI_Group group, int n, int ranges[][3], MPI_Group *newgroup);
int MPI_Group_range_excl(MPI_Group group, int n, int ranges[][3], MPI_Group *newgroup);
int MPI_Group_free(MPI_Group *group);
int MPI_Comm_size(MPI_Comm comm, int *size);
int MPI_Comm_rank(MPI_Comm comm, int *rank);
int MPI_Comm_compare(MPI_Comm comm1, MPI_Comm comm2, int *result);
int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm);
int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm);
int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm);
int MPI_Comm_free(MPI_Comm *comm);
int MPI_Comm_test_inter(MPI_Comm comm, int *flag);
int MPI_Comm_remote_size(MPI_Comm comm, int *size);
int MPI_Comm_remote_group(MPI_Comm comm, MPI_Group *group);
int MPI_Intercomm_create(MPI_Comm local_comm, int local_leader, MPI_Comm peer_comm,
              int remote_leader, int tag, MPI_Comm *newintercomm);
int MPI_Intercomm_merge(MPI_Comm intercomm, int high, MPI_Comm *newintracomm);
int MPI_Keyval_create(MPI_Copy_function *copy_fn, MPI_Delete_function *delete_fn,
              int *keyval, void* extra_state);
int MPI_Keyval_free(int *keyval);
int MPI_Attr_put(MPI_Comm comm, int keyval, void *attribute_val);
int MPI_Attr_get(MPI_Comm comm, int keyval, void *attribute_val, int *flag);
int MPI_Attr_delete(MPI_Comm comm, int keyval);
#ifdef _THREAD_SAFE
int MPI_Comm_create_errhandler(MPI_Comm_errhandler_fn *function, MPI_Errhandler *errhandler);
int MPI_Comm_get_errhandler(MPI_Comm comm, MPI_Errhandler *errhandler);
int MPI_Comm_set_errhandler(MPI_Comm comm, MPI_Errhandler errhandler);
int MPI_Comm_create_keyval(MPI_Comm_copy_attr_function *comm_copy_attr_fn,
              MPI_Comm_delete_attr_function *comm_delete_attr_fn, int *comm_keyval,
              void *extra_state);
int MPI_Comm_delete_attr(MPI_Comm comm, int comm_keyval);
int MPI_Comm_free_keyval(int *comm_keyval);
int MPI_Comm_get_attr(MPI_Comm comm, int comm_keyval, void *attribute_val, int *flag);
int MPI_Comm_set_attr(MPI_Comm comm, int comm_keyval, void *attribute_val);
#endif


   /*   C Bindings for Process Topologies   */

int MPI_Cart_create(MPI_Comm comm_old, int ndims, int *dims, int *periods,
              int reorder, MPI_Comm *comm_cart);
int MPI_Dims_create(int nnodes, int ndims, int *dims);
int MPI_Graph_create(MPI_Comm comm_old, int nnodes, int *index, int *edges,
              int reorder, MPI_Comm *comm_graph);
int MPI_Topo_test(MPI_Comm comm, int *status);
int MPI_Graphdims_get(MPI_Comm comm, int *nnodes, int *nedges);
int MPI_Graph_get(MPI_Comm comm, int maxindex, int maxedges, int *index, int *edges);
int MPI_Cartdim_get(MPI_Comm comm, int *ndims);
int MPI_Cart_get(MPI_Comm comm, int maxdims, int *dims, int *periods, int *coords);
int MPI_Cart_rank(MPI_Comm comm, int *coords, int *rank);
int MPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int *coords);
int MPI_Graph_neighbors_count(MPI_Comm comm, int rank, int *nneighbors);
int MPI_Graph_neighbors(MPI_Comm comm, int rank, int maxneighbors, int *neighbors);
int MPI_Cart_shift(MPI_Comm comm, int direction, int disp, int *rank_source, int *rank_dest);
int MPI_Cart_sub(MPI_Comm comm, int *remain_dims, MPI_Comm *newcomm);
int MPI_Cart_map(MPI_Comm comm, int ndims, int *dims, int *periods, int *newrank);
int MPI_Graph_map(MPI_Comm comm, int nnodes, int *index, int *edges, int *newrank);


  /*    C bindings for Environmental Inquiry   */

int MPI_Get_version(int *version, int *subversion);
int MPI_Get_processor_name(char *name, int *resultlen);
int MPI_Errhandler_create(MPI_Handler_function *function, MPI_Errhandler *errhandler);
int MPI_Errhandler_set(MPI_Comm comm, MPI_Errhandler errhandler);
int MPI_Errhandler_get(MPI_Comm comm, MPI_Errhandler *errhandler);
int MPI_Errhandler_free(MPI_Errhandler *errhandler);
int MPI_Error_string(int errorcode, char *string, int *resultlen);
int MPI_Error_class(int errorcode, int *errorclass);
double MPI_Wtime(void);
double MPI_Wtick(void);
int MPI_Init(int *argc, char ***argv);
int MPI_Finalize(void);
int MPI_Initialized(int *flag);
int MPI_Abort(MPI_Comm comm, int errorcode);


  /*    C Bindings for Profiling  */

int MPI_Pcontrol(const int level,...);

    /*  C Bindings for Non-Blocking Collective Communication   */

int MPE_Ibarrier(MPI_Comm comm,MPI_Request *request);
int MPE_Ibcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm ,MPI_Request *request);
int MPE_Igather(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
              int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,MPI_Request *request) ;
int MPE_Igatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
              int *recvcounts, int *displs, MPI_Datatype recvtype, int root,
              MPI_Comm comm,MPI_Request *request);
int MPE_Iscatter(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
              int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,MPI_Request *request);
int MPE_Iscatterv(void* sendbuf, int *sendcounts, int *displs, MPI_Datatype sendtype,
              void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,MPI_Request *request);
int MPE_Iallgather(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
              int recvcount, MPI_Datatype recvtype, MPI_Comm comm,MPI_Request *request);
int MPE_Iallgatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
              int *recvcounts, int *displs, MPI_Datatype recvtype, MPI_Comm comm,MPI_Request *request);
int MPE_Ialltoall(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
              int recvcount, MPI_Datatype recvtype, MPI_Comm comm,MPI_Request *request);
int MPE_Ialltoallv(void* sendbuf, int *sendcounts, int *sdispls, MPI_Datatype sendtype,
              void* recvbuf, int *recvcounts, int *rdispls, MPI_Datatype recvtype,
              MPI_Comm comm,MPI_Request *request);
int MPE_Ireduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype,
              MPI_Op op, int root, MPI_Comm comm,MPI_Request *request);
int MPE_Iallreduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype,
              MPI_Op op, MPI_Comm comm,MPI_Request *request);
int MPE_Ireduce_scatter(void* sendbuf, void* recvbuf, int *recvcounts,
              MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,MPI_Request *request);
int MPE_Iscan(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype,
              MPI_Op op, MPI_Comm comm ,MPI_Request *request);

#ifdef _THREAD_SAFE
    /*  C Bindings for MPI 1-sided */

int MPI_Win_create(void *base, MPI_Aint size, int disp_unit,
                    MPI_Info info, MPI_Comm comm, MPI_Win *win);
int MPI_Win_free(MPI_Win *win);
int MPI_Win_get_group(MPI_Win win, MPI_Group *group);
int MPI_Put(void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
            int target_rank, MPI_Aint target_disp, int target_count,
            MPI_Datatype target_datatype, MPI_Win win);
int MPI_Get(void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
            int target_rank, MPI_Aint target_disp, int target_count,
            MPI_Datatype target_datatype, MPI_Win win);
int MPI_Accumulate(void *origin_addr, int origin_count, MPI_Datatype
                   origin_datatype, int target_rank, MPI_Aint target_disp,
                   int target_count, MPI_Datatype target_datatype,
                   MPI_Op op, MPI_Win win);
int MPI_Win_fence(int assert, MPI_Win win);
int MPI_Win_start(MPI_Group group, int assert, MPI_Win win);
int MPI_Win_complete(MPI_Win win);
int MPI_Win_post(MPI_Group group, int assert, MPI_Win win);
int MPI_Win_wait(MPI_Win win);
int MPI_Win_test(MPI_Win win, int *flag);
int MPI_Win_lock(int lock_type, int rank, int assert, MPI_Win win);
int MPI_Win_unlock(int rank, MPI_Win win);
int MPI_Win_get_attr(MPI_Win win, int win_keyval, void *attribute_val,
                       int *flag);
int MPI_Win_set_attr(MPI_Win win, int win_keyval, void *attribute_val);
int MPI_Win_create_keyval(MPI_Win_copy_attr_function *win_copy_attr_fn,
              MPI_Win_delete_attr_function *win_delete_attr_fn, int *win_keyval,
              void *extra_state);
int MPI_Win_delete_attr(MPI_Win win, int win_keyval);
int MPI_Win_free_keyval(int *win_keyval);
int MPI_Win_create_errhandler(MPI_Win_errhandler_fn *function, MPI_Errhandler *errhandler);
int MPI_Win_set_errhandler(MPI_Win win, MPI_Errhandler err);
int MPI_Win_get_errhandler(MPI_Win win, MPI_Errhandler *errhandler);

int MPI_Alloc_mem(MPI_Aint size, MPI_Info info, void *baseptr);
int MPI_Free_mem(void *base);

    /*  C Bindings for I/O  */

int MPI_File_open(MPI_Comm comm,char *filename,int amode,MPI_Info info,
                  MPI_File *fh);
int MPI_File_close(MPI_File *fh);
int MPI_File_delete(char *filename,MPI_Info info);
int MPI_File_set_size(MPI_File fh,MPI_Offset size);
int MPI_File_get_size(MPI_File fh,MPI_Offset *size);
int MPI_File_get_group(MPI_File fh,MPI_Group *group);
int MPI_File_get_amode(MPI_File fh,int *amode);
int MPI_File_set_info(MPI_File fh,MPI_Info info);
int MPI_File_get_info(MPI_File fh,MPI_Info *info_used);
int MPI_File_set_view(MPI_File fh,MPI_Offset disp,MPI_Datatype etype,
                      MPI_Datatype filetype,char *datarep,MPI_Info info);
int MPI_File_get_view(MPI_File fh,MPI_Offset *disp,MPI_Datatype *etype,
                      MPI_Datatype *filetype,char *datarep);
int MPI_File_read_at(MPI_File fh,MPI_Offset offset,void *buf,int count,
                     MPI_Datatype datatype,MPI_Status *status);
int MPI_File_read_at_all(MPI_File fh,MPI_Offset offset,void *buf,int count,
                         MPI_Datatype datatype,MPI_Status *status);
int MPI_File_write_at(MPI_File fh,MPI_Offset offset,void *buf,int count,
                      MPI_Datatype datatype,MPI_Status *status);
int MPI_File_write_at_all(MPI_File fh,MPI_Offset offset,void *buf,int count,
                          MPI_Datatype datatype,MPI_Status *status);
int MPI_File_iread_at(MPI_File fh,MPI_Offset offset,void *buf,int count,
                      MPI_Datatype datatype,MPI_Request *request);
int MPI_File_iwrite_at(MPI_File fh,MPI_Offset offset,void *buf,int count,
                       MPI_Datatype datatype,MPI_Request *request);
int MPI_File_get_atomicity(MPI_File fh,int *flag);
int MPI_File_sync(MPI_File fh);
int MPI_Type_create_subarray(int ndims,int array_of_sizes[],
                             int array_of_subsizes[],int array_of_starts[],
                             int order,MPI_Datatype oldtype,
                             MPI_Datatype *newtype);
int MPI_Type_create_darray(int size,int rank,int ndims,int array_of_gsizes[],
                           int array_of_distribs[],int array_of_dargs[],
                           int array_of_psizes[],int order,MPI_Datatype oldtype,
                           MPI_Datatype *newtype);
int MPI_File_create_errhandler(MPI_File_errhandler_fn *function,
                               MPI_Errhandler *errhandler);
int MPI_File_set_errhandler(MPI_File fh,MPI_Errhandler errhandler);
int MPI_File_get_errhandler(MPI_File fh,MPI_Errhandler *errhandler);
int MPI_File_create_errhandler(MPI_File_errhandler_fn *function, MPI_Errhandler *errhandler);
int MPI_File_get_byte_offset(MPI_File fh, MPI_Offset offset, MPI_Offset *disp);
int MPI_File_get_info(MPI_File fh, MPI_Info *info_used);
int MPI_File_get_position(MPI_File fh, MPI_Offset *offset);
int MPI_File_get_position_shared(MPI_File fh, MPI_Offset *offset);
int MPI_File_get_type_extent(MPI_File fh, MPI_Datatype datatype, MPI_Aint *extent);
int MPI_File_iread(MPI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Request *request);
int MPI_File_iread_shared(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
             MPI_Request *request);
int MPI_File_iwrite(MPI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Request *request);
int MPI_File_iwrite_shared(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
             MPI_Request *request);
int MPI_File_preallocate(MPI_File fh, MPI_Offset size);
int MPI_File_read(MPI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Status *status);
int MPI_File_read_all(MPI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Status *status);
int MPI_File_read_all_begin(MPI_File fh, void *buf, int count, MPI_Datatype datatype);
int MPI_File_read_all_end(MPI_File fh, void *buf, MPI_Status *status);
int MPI_File_read_at_all_begin(MPI_File fh, MPI_Offset offset, void *buf, int count,
             MPI_Datatype datatype);
int MPI_File_read_at_all_end(MPI_File fh, void *buf, MPI_Status *status);
int MPI_File_read_ordered(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
             MPI_Status *status);
int MPI_File_read_ordered_begin(MPI_File fh, void *buf, int count, MPI_Datatype datatype);
int MPI_File_read_ordered_end(MPI_File fh, void *buf, MPI_Status *status);
int MPI_File_read_shared(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
             MPI_Status *status);
int MPI_File_seek(MPI_File fh, MPI_Offset offset, int whence);
int MPI_File_seek_shared(MPI_File fh, MPI_Offset offset, int whence);
int MPI_File_set_atomicity(MPI_File fh, int flag);
int MPI_File_set_info(MPI_File fh, MPI_Info info);
int MPI_File_write(MPI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Status *status);
int MPI_File_write_all(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
             MPI_Status *status);
int MPI_File_write_all_begin(MPI_File fh, void *buf, int count, MPI_Datatype datatype);
int MPI_File_write_all_end(MPI_File fh, void *buf, MPI_Status *status);
int MPI_File_write_at_all_begin(MPI_File fh, MPI_Offset offset, void *buf, int count,
             MPI_Datatype datatype);
int MPI_File_write_at_all_end(MPI_File fh, void *buf, MPI_Status *status);
int MPI_File_write_ordered(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
             MPI_Status *status);
int MPI_File_write_ordered_begin(MPI_File fh, void *buf, int count, MPI_Datatype datatype);
int MPI_File_write_ordered_end(MPI_File fh, void *buf, MPI_Status *status);
int MPI_File_write_shared(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
             MPI_Status *status);
int MPI_Register_datarep(char *datarep, MPI_Datarep_conversion_function *read_conversion_fn,
             MPI_Datarep_conversion_function *write_conversion_fn,
             MPI_Datarep_extent_function *dtype_file_extent_fn, void *extra_state);
int MPI_Info_create(MPI_Info *info);
int MPI_Info_set(MPI_Info Info,char *key,char *value);
int MPI_Info_delete(MPI_Info info,char *key);
int MPI_Info_get(MPI_Info info,char *key,int valuelen,char *value,int *flag);
int MPI_Info_get_valuelen(MPI_Info info,char *key,int *valuelen,int *flag);
int MPI_Info_get_nkeys(MPI_Info info,int *nkeys);
int MPI_Info_get_nthkey(MPI_Info info,int n,char *key);
int MPI_Info_dup(MPI_Info info, MPI_Info *newinfo);
int MPI_Info_free(MPI_Info *info);
#endif

/****************************************************************************
 ****************  Profiling (PMPI) Entry Points for C  *********************
 ****************************************************************************/

   /* C Bindings for Point-to-Point Communication */

int PMPI_Send(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
             MPI_Comm comm);
int PMPI_Recv(void* buf, int count, MPI_Datatype datatype, int source, int tag,
             MPI_Comm comm, MPI_Status *status);
int PMPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int *count);
int PMPI_Bsend(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
             MPI_Comm comm);
int PMPI_Ssend(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
             MPI_Comm comm);
int PMPI_Rsend(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
             MPI_Comm comm);
int PMPI_Buffer_attach( void* buffer, int size);
int PMPI_Buffer_detach( void* buffer, int* size);
int PMPI_Isend(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request);
int PMPI_Ibsend(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request);
int PMPI_Issend(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request);
int PMPI_Irsend(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request);
int PMPI_Irecv(void* buf, int count, MPI_Datatype datatype, int source, int tag,
              MPI_Comm comm, MPI_Request *request);
int PMPI_Wait(MPI_Request *request, MPI_Status *status);
int PMPI_Test(MPI_Request *request, int *flag, MPI_Status *status);
int PMPI_Request_free(MPI_Request *request);
int PMPI_Waitany(int count, MPI_Request *array_of_requests, int *index,
              MPI_Status *status);
int PMPI_Testany(int count, MPI_Request *array_of_requests, int *index, int *flag,
              MPI_Status *status);
int PMPI_Waitall(int count, MPI_Request *array_of_requests,
              MPI_Status *array_of_statuses);
int PMPI_Testall(int count, MPI_Request *array_of_requests, int *flag,
              MPI_Status *array_of_statuses);
int PMPI_Waitsome(int incount, MPI_Request *array_of_requests, int *outcount,
              int *array_of_indices, MPI_Status *array_of_statuses);
int PMPI_Testsome(int incount, MPI_Request *array_of_requests, int *outcount,
              int *array_of_indices, MPI_Status *array_of_statuses);
int PMPI_Iprobe(int source, int tag, MPI_Comm comm, int *flag, MPI_Status *status);
int PMPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status *status);
int PMPI_Cancel(MPI_Request *request);
int PMPI_Test_cancelled(MPI_Status *status, int *flag);
int PMPI_Send_init(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request);
int PMPI_Bsend_init(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request);
int PMPI_Ssend_init(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request);
int PMPI_Rsend_init(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request);
int PMPI_Recv_init(void* buf, int count, MPI_Datatype datatype, int source, int tag,
              MPI_Comm comm, MPI_Request *request);
int PMPI_Start(MPI_Request *request);
int PMPI_Startall(int count, MPI_Request *array_of_requests);
int PMPI_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest,
              int sendtag, void *recvbuf, int recvcount, MPI_Datatype recvtype,
              int source, int recvtag, MPI_Comm comm, MPI_Status *status);
int PMPI_Sendrecv_replace(void* buf, int count, MPI_Datatype datatype, int dest,
              int sendtag, int source, int recvtag, MPI_Comm comm, MPI_Status *status);
int PMPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype);
int PMPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype,
              MPI_Datatype *newtype);
int PMPI_Type_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype,
              MPI_Datatype *newtype);
int PMPI_Type_indexed(int count, int *array_of_blocklengths, int *array_of_displacements,
              MPI_Datatype oldtype, MPI_Datatype *newtype);
int PMPI_Type_hindexed(int count, int *array_of_blocklengths, MPI_Aint *array_of_displacements,
              MPI_Datatype oldtype, MPI_Datatype *newtype);
int PMPI_Type_struct(int count, int *array_of_blocklengths, MPI_Aint *array_of_displacements,
              MPI_Datatype *array_of_types, MPI_Datatype *newtype);

#ifdef _THREAD_SAFE
int PMPI_Type_get_envelope(MPI_Datatype datatype, int *num_integers, int *num_addresses,
             int *num_datatypes, int *combiner );
int PMPI_Type_get_contents(MPI_Datatype datatype, int max_integers, int max_addresses,
             int max_datatypes, int array_of_integers[], MPI_Aint array_of_addresses[],
             MPI_Datatype array_of_datatypes[] );
int PMPI_Type_create_keyval(MPI_Type_copy_attr_function *type_copy_attr_fn,
             MPI_Type_delete_attr_function *type_delete_attr_fn, int *type_keyval,
             void *extra_state);
int PMPI_Type_delete_attr(MPI_Datatype type, int type_keyval);
int PMPI_Type_free_keyval(int *type_keyval);
int PMPI_Type_get_attr(MPI_Datatype type, int type_keyval, void *attribute_val, int *flag);
int PMPI_Type_set_attr(MPI_Datatype type, int type_keyval, void *attribute_val);
int PMPI_Type_dup(MPI_Datatype type, MPI_Datatype *newtype);
#endif

int PMPI_Address(void* location, MPI_Aint *address);
int PMPI_Type_extent(MPI_Datatype datatype, MPI_Aint *extent);
int PMPI_Type_size(MPI_Datatype datatype, int *size);
int PMPI_Type_lb(MPI_Datatype datatype, MPI_Aint *displacement);
int PMPI_Type_ub(MPI_Datatype datatype, MPI_Aint *displacement);
int PMPI_Type_commit(MPI_Datatype *datatype);
int PMPI_Type_free(MPI_Datatype *datatype);
int PMPI_Get_elements(MPI_Status *status, MPI_Datatype datatype, int *count);
int PMPI_Pack(void* inbuf, int incount, MPI_Datatype datatype, void *outbuf, int outsize,
              int *position,  MPI_Comm comm);
int PMPI_Unpack(void* inbuf, int insize, int *position, void *outbuf, int outcount,
              MPI_Datatype datatype, MPI_Comm comm);
int PMPI_Pack_size(int incount, MPI_Datatype datatype, MPI_Comm comm, int *size);

    /*  C Bindings for Collective Communication   */

int PMPI_Barrier(MPI_Comm comm );
int PMPI_Bcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm );
int PMPI_Gather(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
              int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) ;
int PMPI_Gatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
              int *recvcounts, int *displs, MPI_Datatype recvtype, int root,
              MPI_Comm comm);
int PMPI_Scatter(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
              int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
int PMPI_Scatterv(void* sendbuf, int *sendcounts, int *displs, MPI_Datatype sendtype,
              void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
int PMPI_Allgather(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
              int recvcount, MPI_Datatype recvtype, MPI_Comm comm);
int PMPI_Allgatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
              int *recvcounts, int *displs, MPI_Datatype recvtype, MPI_Comm comm);
int PMPI_Alltoall(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
              int recvcount, MPI_Datatype recvtype, MPI_Comm comm);
int PMPI_Alltoallv(void* sendbuf, int *sendcounts, int *sdispls, MPI_Datatype sendtype,
              void* recvbuf, int *recvcounts, int *rdispls, MPI_Datatype recvtype,
              MPI_Comm comm);
int PMPI_Reduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype,
              MPI_Op op, int root, MPI_Comm comm);
int PMPI_Op_create(MPI_User_function *function, int commute, MPI_Op *op);
int PMPI_Op_free( MPI_Op *op);
int PMPI_Allreduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype,
              MPI_Op op, MPI_Comm comm);
int PMPI_Reduce_scatter(void* sendbuf, void* recvbuf, int *recvcounts,
              MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int PMPI_Scan(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype,
              MPI_Op op, MPI_Comm comm );

    /* C  Bindings for Groups, Contexts, and Communicators  */

int PMPI_Group_size(MPI_Group group, int *size);
int PMPI_Group_rank(MPI_Group group, int *rank);
int PMPI_Group_translate_ranks (MPI_Group group1, int n, int *ranks1,
              MPI_Group group2, int *ranks2);
int PMPI_Group_compare(MPI_Group group1,MPI_Group group2, int *result);
int PMPI_Comm_group(MPI_Comm comm, MPI_Group *group);
int PMPI_Group_union(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup);
int PMPI_Group_intersection(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup);
int PMPI_Group_difference(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup);
int PMPI_Group_incl(MPI_Group group, int n, int *ranks, MPI_Group *newgroup);
int PMPI_Group_excl(MPI_Group group, int n, int *ranks, MPI_Group *newgroup);
int PMPI_Group_range_incl(MPI_Group group, int n, int ranges[][3], MPI_Group *newgroup);
int PMPI_Group_range_excl(MPI_Group group, int n, int ranges[][3], MPI_Group *newgroup);
int PMPI_Group_free(MPI_Group *group);
int PMPI_Comm_size(MPI_Comm comm, int *size);
int PMPI_Comm_rank(MPI_Comm comm, int *rank);
int PMPI_Comm_compare(MPI_Comm comm1, MPI_Comm comm2, int *result);
int PMPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm);
int PMPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm);
int PMPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm);
int PMPI_Comm_free(MPI_Comm *comm);
int PMPI_Comm_test_inter(MPI_Comm comm, int *flag);
int PMPI_Comm_remote_size(MPI_Comm comm, int *size);
int PMPI_Comm_remote_group(MPI_Comm comm, MPI_Group *group);
int PMPI_Intercomm_create(MPI_Comm local_comm, int local_leader, MPI_Comm peer_comm,
              int remote_leader, int tag, MPI_Comm *newintercomm);
int PMPI_Intercomm_merge(MPI_Comm intercomm, int high, MPI_Comm *newintracomm);
int PMPI_Keyval_create(MPI_Copy_function *copy_fn, MPI_Delete_function *delete_fn,
              int *keyval, void* extra_state);
int PMPI_Keyval_free(int *keyval);
int PMPI_Attr_put(MPI_Comm comm, int keyval, void *attribute_val);
int PMPI_Attr_get(MPI_Comm comm, int keyval, void *attribute_val, int *flag);
int PMPI_Attr_delete(MPI_Comm comm, int keyval);

#ifdef _THREAD_SAFE
int PMPI_Comm_create_errhandler(MPI_Comm_errhandler_fn *function, MPI_Errhandler *errhandler);
int PMPI_Comm_get_errhandler(MPI_Comm comm, MPI_Errhandler *errhandler);
int PMPI_Comm_set_errhandler(MPI_Comm comm, MPI_Errhandler errhandler);
int PMPI_Comm_create_keyval(MPI_Comm_copy_attr_function *comm_copy_attr_fn,
             MPI_Comm_delete_attr_function *comm_delete_attr_fn, int *comm_keyval,
             void *extra_state);
int PMPI_Comm_delete_attr(MPI_Comm comm, int comm_keyval);
int PMPI_Comm_free_keyval(int *comm_keyval);
int PMPI_Comm_get_attr(MPI_Comm comm, int comm_keyval, void *attribute_val, int *flag);
int PMPI_Comm_set_attr(MPI_Comm comm, int comm_keyval, void *attribute_val);
#endif


   /*   C Bindings for Process Topologies   */

int PMPI_Cart_create(MPI_Comm comm_old, int ndims, int *dims, int *periods,
              int reorder, MPI_Comm *comm_cart);
int PMPI_Dims_create(int nnodes, int ndims, int *dims);
int PMPI_Graph_create(MPI_Comm comm_old, int nnodes, int *index, int *edges,
              int reorder, MPI_Comm *comm_graph);
int PMPI_Topo_test(MPI_Comm comm, int *status);
int PMPI_Graphdims_get(MPI_Comm comm, int *nnodes, int *nedges);
int PMPI_Graph_get(MPI_Comm comm, int maxindex, int maxedges, int *index, int *edges);
int PMPI_Cartdim_get(MPI_Comm comm, int *ndims);
int PMPI_Cart_get(MPI_Comm comm, int maxdims, int *dims, int *periods, int *coords);
int PMPI_Cart_rank(MPI_Comm comm, int *coords, int *rank);
int PMPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int *coords);
int PMPI_Graph_neighbors_count(MPI_Comm comm, int rank, int *nneighbors);
int PMPI_Graph_neighbors(MPI_Comm comm, int rank, int maxneighbors, int *neighbors);
int PMPI_Cart_shift(MPI_Comm comm, int direction, int disp, int *rank_source, int *rank_dest);
int PMPI_Cart_sub(MPI_Comm comm, int *remain_dims, MPI_Comm *newcomm);
int PMPI_Cart_map(MPI_Comm comm, int ndims, int *dims, int *periods, int *newrank);
int PMPI_Graph_map(MPI_Comm comm, int nnodes, int *index, int *edges, int *newrank);


  /*    C bindings for Environmental Inquiry   */

int PMPI_Get_version(int *version, int *subversion);
int PMPI_Get_processor_name(char *name, int *resultlen);
int PMPI_Errhandler_create(MPI_Handler_function *function, MPI_Errhandler *errhandler);
int PMPI_Errhandler_set(MPI_Comm comm, MPI_Errhandler errhandler);
int PMPI_Errhandler_get(MPI_Comm comm, MPI_Errhandler *errhandler);
int PMPI_Errhandler_free(MPI_Errhandler *errhandler);
int PMPI_Error_string(int errorcode, char *string, int *resultlen);
int PMPI_Error_class(int errorcode, int *errorclass);
double PMPI_Wtime(void);
double PMPI_Wtick(void);
int PMPI_Init(int *argc, char ***argv);
int PMPI_Finalize(void);
int PMPI_Initialized(int *flag);
int PMPI_Abort(MPI_Comm comm, int errorcode);


  /*    C Bindings for Profiling  */

int PMPI_Pcontrol(const int level,...);

    /*  C Bindings for Non-Blocking Collective Communication   */

int PMPE_Ibarrier(MPI_Comm comm,MPI_Request *request);
int PMPE_Ibcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm ,MPI_Request *request);
int PMPE_Igather(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
              int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,MPI_Request *request) ;
int PMPE_Igatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
              int *recvcounts, int *displs, MPI_Datatype recvtype, int root,
              MPI_Comm comm,MPI_Request *request);
int PMPE_Iscatter(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
              int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,MPI_Request *request);
int PMPE_Iscatterv(void* sendbuf, int *sendcounts, int *displs, MPI_Datatype sendtype,
              void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,MPI_Request *request);
int PMPE_Iallgather(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
              int recvcount, MPI_Datatype recvtype, MPI_Comm comm,MPI_Request *request);
int PMPE_Iallgatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
              int *recvcounts, int *displs, MPI_Datatype recvtype, MPI_Comm comm,MPI_Request *request);
int PMPE_Ialltoall(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
              int recvcount, MPI_Datatype recvtype, MPI_Comm comm,MPI_Request *request);
int PMPE_Ialltoallv(void* sendbuf, int *sendcounts, int *sdispls, MPI_Datatype sendtype,
              void* recvbuf, int *recvcounts, int *rdispls, MPI_Datatype recvtype,
              MPI_Comm comm,MPI_Request *request);
int PMPE_Ireduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype,
              MPI_Op op, int root, MPI_Comm comm,MPI_Request *request);
int PMPE_Iallreduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype,
              MPI_Op op, MPI_Comm comm,MPI_Request *request);
int PMPE_Ireduce_scatter(void* sendbuf, void* recvbuf, int *recvcounts,
              MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,MPI_Request *request);
int PMPE_Iscan(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype,
              MPI_Op op, MPI_Comm comm ,MPI_Request *request);

#ifdef _THREAD_SAFE
    /*  C Bindings for MPI 1-sided */

int PMPI_Win_create(void *base, MPI_Aint size, int disp_unit,
                    MPI_Info info, MPI_Comm comm, MPI_Win *win);
int PMPI_Win_free(MPI_Win *win);
int PMPI_Win_get_group(MPI_Win win, MPI_Group *group);
int PMPI_Put(void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
            int target_rank, MPI_Aint target_disp, int target_count,
            MPI_Datatype target_datatype, MPI_Win win);
int PMPI_Get(void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
            int target_rank, MPI_Aint target_disp, int target_count,
            MPI_Datatype target_datatype, MPI_Win win);
int PMPI_Accumulate(void *origin_addr, int origin_count, MPI_Datatype
                   origin_datatype, int target_rank, MPI_Aint target_disp,
                   int target_count, MPI_Datatype target_datatype,
                   MPI_Op op, MPI_Win win);
int PMPI_Win_fence(int assert, MPI_Win win);
int PMPI_Win_start(MPI_Group group, int assert, MPI_Win win);
int PMPI_Win_complete(MPI_Win win);
int PMPI_Win_post(MPI_Group group, int assert, MPI_Win win);
int PMPI_Win_wait(MPI_Win win);
int PMPI_Win_test(MPI_Win win, int *flag);
int PMPI_Win_lock(int lock_type, int rank, int assert, MPI_Win win);
int PMPI_Win_unlock(int rank, MPI_Win win);
int PMPI_Win_get_attr(MPI_Win win, int win_keyval, void *attribute_val, int *flag);
int PMPI_Win_set_attr(MPI_Win win, int win_keyval, void *attribute_val);

int PMPI_Alloc_mem(MPI_Aint size, MPI_Info info, void *baseptr);
int PMPI_Free_mem(void *base);

    /*  C Bindings for I/O  */

int PMPI_File_open(MPI_Comm comm,char *filename,int amode,MPI_Info info,
                   MPI_File *fh);
int PMPI_File_close(MPI_File *fh);
int PMPI_File_delete(char *filename,MPI_Info info);
int PMPI_File_set_size(MPI_File fh,MPI_Offset size);
int PMPI_File_get_size(MPI_File fh,MPI_Offset *size);
int PMPI_File_get_group(MPI_File fh,MPI_Group *group);
int PMPI_File_get_amode(MPI_File fh,int *amode);
int PMPI_File_set_info(MPI_File fh,MPI_Info info);
int PMPI_File_get_info(MPI_File fh,MPI_Info *info_used);
int PMPI_File_set_view(MPI_File fh,MPI_Offset disp,MPI_Datatype etype,
                       MPI_Datatype filetype,char *datarep,MPI_Info info);
int PMPI_File_get_view(MPI_File fh,MPI_Offset *disp,MPI_Datatype *etype,
                       MPI_Datatype *filetype,char *datarep);
int PMPI_File_read_at(MPI_File fh,MPI_Offset offset,void *buf,int count,
                      MPI_Datatype datatype,MPI_Status *status);
int PMPI_File_read_at_all(MPI_File fh,MPI_Offset offset,void *buf,int count,
                          MPI_Datatype datatype,MPI_Status *status);
int PMPI_File_write_at(MPI_File fh,MPI_Offset offset,void *buf,int count,
                       MPI_Datatype datatype,MPI_Status *status);
int PMPI_File_write_at_all(MPI_File fh,MPI_Offset offset,void *buf,int count,
                           MPI_Datatype datatype,MPI_Status *status);
int PMPI_File_iread_at(MPI_File fh,MPI_Offset offset,void *buf,int count,
                       MPI_Datatype datatype,MPI_Request *request);
int PMPI_File_iwrite_at(MPI_File fh,MPI_Offset offset,void *buf,int count,
                        MPI_Datatype datatype,MPI_Request *request);
int PMPI_File_get_atomicity(MPI_File fh,int *flag);
int PMPI_File_sync(MPI_File fh);
int PMPI_Type_create_subarray(int ndims,int array_of_sizes[],
                              int array_of_subsizes[],int array_of_starts[],
                              int order,MPI_Datatype oldtype,
                              MPI_Datatype *newtype);
int PMPI_Type_create_darray(int size,int rank,int ndims,int array_of_gsizes[],
                            int array_of_distribs[],int array_of_dargs[],
                            int array_of_psizes[],int order,MPI_Datatype oldtype,
                            MPI_Datatype *newtype);
int PMPI_File_create_errhandler(MPI_File_errhandler_fn *function,
                                MPI_Errhandler *errhandler);
int PMPI_File_set_errhandler(MPI_File fh,MPI_Errhandler errhandler);
int PMPI_File_get_errhandler(MPI_File fh,MPI_Errhandler *errhandler);
int PMPI_File_create_errhandler(MPI_File_errhandler_fn *function, MPI_Errhandler *errhandler);
int PMPI_File_get_byte_offset(MPI_File fh, MPI_Offset offset, MPI_Offset *disp);
int PMPI_File_get_info(MPI_File fh, MPI_Info *info_used);
int PMPI_File_get_position(MPI_File fh, MPI_Offset *offset);
int PMPI_File_get_position_shared(MPI_File fh, MPI_Offset *offset);
int PMPI_File_get_type_extent(MPI_File fh, MPI_Datatype datatype, MPI_Aint *extent);
int PMPI_File_iread(MPI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Request *request);
int PMPI_File_iread_shared(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
             MPI_Request *request);
int PMPI_File_iwrite(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
             MPI_Request *request);
int PMPI_File_iwrite_shared(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
             MPI_Request *request);
int PMPI_File_preallocate(MPI_File fh, MPI_Offset size);
int PMPI_File_read(MPI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Status *status);
int PMPI_File_read_all(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
             MPI_Status *status);
int PMPI_File_read_all_begin(MPI_File fh, void *buf, int count, MPI_Datatype datatype);
int PMPI_File_read_all_end(MPI_File fh, void *buf, MPI_Status *status);
int PMPI_File_read_at_all_begin(MPI_File fh, MPI_Offset offset, void *buf, int count,
             MPI_Datatype datatype);
int PMPI_File_read_at_all_end(MPI_File fh, void *buf, MPI_Status *status);
int PMPI_File_read_ordered(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
             MPI_Status *status);
int PMPI_File_read_ordered_begin(MPI_File fh, void *buf, int count, MPI_Datatype datatype);
int PMPI_File_read_ordered_end(MPI_File fh, void *buf, MPI_Status *status);
int PMPI_File_read_shared(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
             MPI_Status *status);
int PMPI_File_seek(MPI_File fh, MPI_Offset offset, int whence);
int PMPI_File_seek_shared(MPI_File fh, MPI_Offset offset, int whence);
int PMPI_File_set_atomicity(MPI_File fh, int flag);
int PMPI_File_set_info(MPI_File fh, MPI_Info info);
int PMPI_File_write(MPI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Status *status);
int PMPI_File_write_all(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
             MPI_Status *status);
int PMPI_File_write_all_begin(MPI_File fh, void *buf, int count, MPI_Datatype datatype);
int PMPI_File_write_all_end(MPI_File fh, void *buf, MPI_Status *status);
int PMPI_File_write_at_all_begin(MPI_File fh, MPI_Offset offset, void *buf, int count,
             MPI_Datatype datatype);
int PMPI_File_write_at_all_end(MPI_File fh, void *buf, MPI_Status *status);
int PMPI_File_write_ordered(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
             MPI_Status *status);
int PMPI_File_write_ordered_begin(MPI_File fh, void *buf, int count, MPI_Datatype datatype);
int PMPI_File_write_ordered_end(MPI_File fh, void *buf, MPI_Status *status);
int PMPI_File_write_shared(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
             MPI_Status *status);
int PMPI_Register_datarep(char *datarep, MPI_Datarep_conversion_function *read_conversion_fn,
             MPI_Datarep_conversion_function *write_conversion_fn,
             MPI_Datarep_extent_function *dtype_file_extent_fn, void *extra_state);
int PMPI_Info_create(MPI_Info *info);
int PMPI_Info_set(MPI_Info Info,char *key,char *value);
int PMPI_Info_delete(MPI_Info info,char *key);
int PMPI_Info_get(MPI_Info info,char *key,int valuelen,char *value,int *flag);
int PMPI_Info_get_valuelen(MPI_Info info,char *key,int *valuelen,int *flag);
int PMPI_Info_get_nkeys(MPI_Info info,int *nkeys);
int PMPI_Info_get_nthkey(MPI_Info info,int n,char *key);
int PMPI_Info_dup(MPI_Info info, MPI_Info *newinfo);
int PMPI_Info_free(MPI_Info *info);
#endif

#if defined(__cplusplus)
   }
#endif
#endif /* _H_MPI */
