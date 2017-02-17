#include <gbpLib.h>
#include <mpi.h>

void __SID_Bcast__(void *buffer, int data_size, int source_rank, SID_Comm *comm);
void __SID_Reduce__(void *sendbuf,void *recvbuf,int count,SID_Datatype datatype,SID_Op op,int root,SID_Comm *comm);
void __SID_Allreduce__(void *sendbuf,void *recvbuf,int count,SID_Datatype datatype,SID_Op op,SID_Comm *comm);
void __SID_Send__(void        *sendbuf,
                  int          sendcount,
                  SID_Datatype sendtype,
                  int          dest,
                  int          sendtag,
                  SID_Comm    *comm);
void __SID_Recv__(void        *recvbuf,
                  int          recvcount,
                  SID_Datatype recvtype,
                  int          source,
                  int          recvtag,
                  SID_Comm    *comm);
void __SID_Sendrecv__(void        *sendbuf,
                      int          sendcount,
                      SID_Datatype sendtype,
                      int          dest,
                      int          sendtag,
                      void        *recvbuf,
                      int          recvcount,
                      SID_Datatype recvtype,
                      int          source,
                      int          recvtag,
                      SID_Comm    *comm);

OMPI_DECLSPEC int __MPI_Bcast__(void *buffer, int count, MPI_Datatype datatype,
                                int root, MPI_Comm comm);
OMPI_DECLSPEC int __MPI_Send__(const void *buf, int count, MPI_Datatype datatype, int dest,
                               int tag, MPI_Comm comm);
OMPI_DECLSPEC int __MPI_Sendrecv__(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                   int dest, int sendtag, void *recvbuf, int recvcount,
                                   MPI_Datatype recvtype, int source, int recvtag,
                                   MPI_Comm comm,  MPI_Status *status);
OMPI_DECLSPEC int __MPI_Recv__(void *buf, int count, MPI_Datatype datatype, int source,
                               int tag, MPI_Comm comm, MPI_Status *status);
OMPI_DECLSPEC int __MPI_Reduce__(const void *sendbuf, void *recvbuf, int count,
                                 MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
OMPI_DECLSPEC int __MPI_Allreduce__(const void *sendbuf, void *recvbuf, int count,
                                    MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
OMPI_DECLSPEC int __MPI_Allgather__(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                    void *recvbuf, int recvcount,
                                    MPI_Datatype recvtype, MPI_Comm comm);
