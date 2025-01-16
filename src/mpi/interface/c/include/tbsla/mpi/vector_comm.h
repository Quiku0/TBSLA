#ifndef TBSLA_CINTERFACE_MPI_vector_comm
#define TBSLA_CINTERFACE_MPI_vector_comm
#include <stdbool.h>
#include <tbsla/cpp/Vector.h>
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

C_CPP_Vector_t *C_MPI_allgatherv(MPI_Comm comm, C_CPP_Vector_t *v, long long int bn_row, long long int lgr);
C_CPP_Vector_t *C_MPI_reduce_sum(MPI_Comm comm, C_CPP_Vector_t *v, long long int n);
C_CPP_Vector_t *C_MPI_reduce_gather(MPI_Comm comm, C_CPP_Vector_t *v, long long int bn_row, long long int bn_col, long long int lpr, long long int lpc, long long int lgr, long long int lgc);
C_CPP_Vector_t *C_MPI_redistribute(MPI_Comm comm, C_CPP_Vector_t *v, long long int bn_row, long long int bn_col, long long int lpr, long long int lpc, long long int lgr, long long int lgc);

#ifdef __cplusplus
}
#endif

#endif /* TBSLA_CINTERFACE_MPI_vector_comm */
