#ifndef TBSLA_PETSC_Matrix
#define TBSLA_PETSC_Matrix

#include <vector>

#include <mpi.h>
#include <petscmat.h>
#include <petscvec.h>

namespace tbsla { namespace petsc {

class Matrix {
  public:
    void fill_cdiag(MPI_Comm comm, long long int nr, long long int nc, long long int cdiag);
    void fill_cqmat(MPI_Comm comm, long long int n_row, long long int n_col, long long int c, double q, unsigned long long int seed_mult);
    Vec spmv(MPI_Comm comm, Vec &v);
    Vec a_axpx_(MPI_Comm comm, Vec &v);
    long long int get_n_row() { return n_row; }
    long long int get_n_col() { return n_col; }

  protected:
    Mat m;
    long long int n_row;
    long long int n_col;
};

}}

#endif
