#include <tbsla/petsc/Matrix.hpp>
#include <tbsla/cpp/utils/range.hpp>
#include <tbsla/cpp/utils/values_generation.hpp>
#include <tbsla/cpp/utils/vector.hpp>

void tbsla::petsc::Matrix::fill_cdiag(MPI_Comm comm, long long int n_row, long long int n_col, long long int cdiag) {
  long long int world, rank;
  MPI_Comm_size(comm, &world);
  MPI_Comm_rank(comm, &rank);

  this->n_row = n_row;
  this->n_col = n_col;

  long long int nv = std::max(std::min(n_row, n_col - cdiag), 0) + std::max(std::min(n_row - cdiag, n_col), 0);
  if(cdiag == 0)
    nv /= 2;

  MatCreate(comm, &m);
  MatSetSizes(m, PETSC_DECIDE, PETSC_DECIDE, n_row, n_col);
  MatSetUp(m);

  if(nv != 0) {
    long long int s = tbsla::utils::range::pflv(nv, rank, world);
    long long int n = tbsla::utils::range::lnv(nv, rank, world);
    for(long long int k = s; k < s + n; k++) {
      auto tuple = tbsla::utils::values_generation::cdiag_value(k, nv, n_row, n_col, cdiag);
      long long int i = std::get<0>(tuple);
      long long int j = std::get<1>(tuple);
      double v = std::get<2>(tuple);
      MatSetValues(m, 1, &i, 1, &j, &v, ADD_VALUES);
    }
  }
  MatAssemblyBegin(m, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(m, MAT_FINAL_ASSEMBLY);
}

void tbsla::petsc::Matrix::fill_cqmat(MPI_Comm comm, long long int n_row, long long int n_col, long long int c, double q, unsigned long long int seed_mult) {
  long long int world, rank;
  MPI_Comm_size(comm, &world);
  MPI_Comm_rank(comm, &rank);

  this->n_row = n_row;
  this->n_col = n_col;

  long long int nv = 0;
  for(long long int i = 0; i < std::min(n_col - std::min(c, n_col) + 1, n_row); i++) {
    nv += std::min(c, n_col);
  }
  for(long long int i = 0; i < std::min(n_row, n_col) - std::min(n_col - std::min(c, n_col) + 1, n_row); i++) {
    nv += std::min(c, n_col) - i - 1;
  }

  MatCreate(comm, &m);
  MatSetSizes(m, PETSC_DECIDE, PETSC_DECIDE, n_row, n_col);
  MatSetUp(m);

  if(nv != 0) {
    long long int s = tbsla::utils::range::pflv(nv, rank, world);
    long long int n = tbsla::utils::range::lnv(nv, rank, world);

    for(long long int k = s; k < s + n; k++) {
      auto tuple = tbsla::utils::values_generation::cqmat_value(k, n_row, n_col, c, q, seed_mult);
      long long int i = std::get<0>(tuple);
      long long int j = std::get<1>(tuple);
      double v = std::get<2>(tuple);
      MatSetValues(m, 1, &i, 1, &j, &v, ADD_VALUES);
    }
  }

  MatAssemblyBegin(m, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(m, MAT_FINAL_ASSEMBLY);
}

Vec tbsla::petsc::Matrix::spmv(MPI_Comm comm, Vec &v) {
  Vec r;
  VecCreate(comm, &r);
  VecSetSizes(r, PETSC_DECIDE, this->n_row);
  VecSetUp(r);
  MatMult(m, v, r);
  return r;
}

Vec tbsla::petsc::Matrix::a_axpx_(MPI_Comm comm, Vec &v) {
  Vec r = spmv(comm, v);
  VecAXPY(r, 1, v);
  return spmv(comm, r);
}
