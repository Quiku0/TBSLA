#ifndef TBSLA_CINTERFACE_MPI_MatrixDENSE
#define TBSLA_CINTERFACE_MPI_MatrixDENSE
#include <stdbool.h>
#include <tbsla/cpp/Vector.h>
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

struct C_MPI_MatrixDENSE;
typedef struct C_MPI_MatrixDENSE C_MPI_MatrixDENSE_t;

C_MPI_MatrixDENSE_t *C_MPI_MatrixDENSE_create();
void C_MPI_MatrixDENSE_destroy(C_MPI_MatrixDENSE_t *m);

void C_MPI_MatrixDENSE_fill_cdiag(C_MPI_MatrixDENSE_t *m, long long int n_row, long long int n_col, long long int cdiag, long long int pr, long long int pc, long long int NR, long long int NC);
void C_MPI_MatrixDENSE_fill_cqmat(C_MPI_MatrixDENSE_t *m, long long int n_row, long long int n_col, long long int c, double q, unsigned long long int seed_mult, long long int pr, long long int pc, long long int NR, long long int NC);
bool C_MPI_MatrixDENSE_read(C_MPI_MatrixDENSE_t *m, char *filename, long long int seek);
bool C_MPI_MatrixDENSE_write(C_MPI_MatrixDENSE_t *m, char *filename);
void C_MPI_MatrixDENSE_print(C_MPI_MatrixDENSE_t *m);
C_CPP_Vector_t *C_MPI_MatrixDENSE_spmv(C_MPI_MatrixDENSE_t *m, MPI_Comm comm, C_CPP_Vector_t *v);

#ifdef __cplusplus
}
#endif

#endif /* TBSLA_C_MatrixDENSE */
