#include "logger.h"
#include <zlog.h>

void __SID_Bcast__(void *buffer, int data_size, int source_rank, SID_Comm *comm)
{
  zlog_category_t *log_cat = zlog_get_category("mpi");
  zlog_info(log_cat, "Bcast : %d bytes from rank %d", data_size, source_rank);
  SID_Bcast(buffer, data_size, source_rank, comm);
}

void __SID_Reduce__(void *sendbuf, void *recvbuf, int count, SID_Datatype datatype, SID_Op op, int root, SID_Comm *comm)
{
  zlog_category_t *log_cat = zlog_get_category("mpi");
  zlog_info(log_cat, "Reduce : %zu bytes to rank %d", count*sizeof(datatype), root);
  SID_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
}

void __SID_Allreduce__(void *sendbuf,void *recvbuf,int count,SID_Datatype datatype,SID_Op op,SID_Comm *comm)
{
  zlog_category_t *log_cat = zlog_get_category("mpi");
  zlog_info(log_cat, "Allreduce : %zu bytes", count*sizeof(datatype));
  SID_Allreduce(sendbuf,recvbuf,count,datatype,op,comm);
}

void __SID_Send__(void *sendbuf, int sendcount, SID_Datatype sendtype, int dest, int sendtag, SID_Comm *comm)
{
  zlog_category_t *log_cat = zlog_get_category("mpi");
  zlog_info(log_cat, "Send : %zu bytes to %d", sendcount*sizeof(sendtype), dest);
  SID_Send(sendbuf, sendcount, sendtype, dest, sendtag, comm);
}


void __SID_Recv__(void         *recvbuf,
    int           recvcount,
    SID_Datatype  recvtype,
    int           source,
    int           recvtag,
    SID_Comm     *comm)
{
  zlog_category_t *log_cat = zlog_get_category("mpi");
  zlog_info(log_cat, "Recv : %zu bytes from %d", recvcount*sizeof(recvtype), source);
  SID_Recv(recvbuf, recvcount, recvtype, source, recvtag, comm);
}


void __SID_Sendrecv__(void         *sendbuf,
    int           sendcount,
    SID_Datatype  sendtype,
    int           dest,
    int           sendtag,
    void         *recvbuf,
    int           recvcount,
    SID_Datatype  recvtype,
    int           source,
    int           recvtag,
    SID_Comm     *comm)
{
  zlog_category_t *log_cat = zlog_get_category("mpi");
  zlog_info(log_cat, "Sendrecv : %zu bytes to %d, %zu bytes from %d", sendcount*sizeof(sendtype), dest, recvcount*sizeof(recvtype), source);
  SID_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm);
}

OMPI_DECLSPEC  int __MPI_Bcast__(void *buffer, int count, MPI_Datatype datatype,
    int root, MPI_Comm comm)
{
  zlog_category_t *log_cat = zlog_get_category("mpi");
  zlog_info(log_cat, "Bcast : %zu bytes from rank %d", count*sizeof(datatype), root);
  MPI_Bcast(buffer, count, datatype, root, comm);
}

OMPI_DECLSPEC  int __MPI_Send__(const void *buf, int count, MPI_Datatype datatype, int dest,
    int tag, MPI_Comm comm)
{
  zlog_category_t *log_cat = zlog_get_category("mpi");
  zlog_info(log_cat, "Send : %zu bytes to %d", count*sizeof(datatype), dest);
  MPI_Send(buf, count, datatype, dest, tag, comm);
}

OMPI_DECLSPEC  int __MPI_Sendrecv__(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
    int dest, int sendtag, void *recvbuf, int recvcount,
    MPI_Datatype recvtype, int source, int recvtag,
    MPI_Comm comm,  MPI_Status *status)
{
  zlog_category_t *log_cat = zlog_get_category("mpi");
  zlog_info(log_cat, "Sendrecv : %zu bytes to %d, %zu bytes from %d", sendcount*sizeof(sendtype), dest, recvcount*sizeof(recvtype), source);
  MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status);
}

OMPI_DECLSPEC  int __MPI_Recv__(void *buf, int count, MPI_Datatype datatype, int source,
    int tag, MPI_Comm comm, MPI_Status *status)
{
  zlog_category_t *log_cat = zlog_get_category("mpi");
  zlog_info(log_cat, "Recv : %zu bytes from %d", count*sizeof(datatype), source);
  MPI_Recv(buf, count, datatype, source, tag, comm, status);
}

OMPI_DECLSPEC  int __MPI_Reduce__(const void *sendbuf, void *recvbuf, int count,
    MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
{
  zlog_category_t *log_cat = zlog_get_category("mpi");
  zlog_info(log_cat, "Reduce : %zu bytes to rank %d", count*sizeof(datatype), root);
  MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
}

OMPI_DECLSPEC  int __MPI_Allreduce__(const void *sendbuf, void *recvbuf, int count,
    MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  zlog_category_t *log_cat = zlog_get_category("mpi");
  zlog_info(log_cat, "Allreduce : %zu bytes", count*sizeof(datatype));
  MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
}

OMPI_DECLSPEC  int __MPI_Allgather__(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
    void *recvbuf, int recvcount,
    MPI_Datatype recvtype, MPI_Comm comm)
{
  zlog_category_t *log_cat = zlog_get_category("mpi");
  zlog_info(log_cat, "Allgather : send %zu bytes, recv %zu bytes", sendcount*sizeof(recvtype), sendcount*sizeof(recvtype));
  MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
}
