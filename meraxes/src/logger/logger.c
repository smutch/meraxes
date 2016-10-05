#include "logger.h"
#include <zlog.h>

void __SID_Bcast__(void *buffer, int data_size, int source_rank, SID_Comm *comm)
{
  zlog_category_t *log_cat = zlog_get_category("mpi");
  zlog_info(log_cat, "Bcast : %d bytes from rank %d", data_size, source_rank);
  SID_Bcast(buffer, data_size, source_rank, comm);
}
