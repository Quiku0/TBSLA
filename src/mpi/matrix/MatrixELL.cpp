#include <tbsla/mpi/MatrixELL.hpp>
#include <tbsla/cpp/utils/range.hpp>
#include <algorithm>
#include <vector>
#include <mpi.h>

#define TBSLA_MATRIX_ELL_READLINES 2048

long long int tbsla::mpi::MatrixELL::read_bin_mpiio(MPI_Comm comm, std::string filename, long long int pr, long long int pc, long long int NR, long long int NC) {
  int world, rank;
  MPI_Comm_size(comm, &world);
  MPI_Comm_rank(comm, &rank);

  MPI_File fh;
  MPI_Status status;
  MPI_File_open(comm, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

  MPI_File_read_all(fh, &this->n_row, 1, MPI_INT, &status);
  MPI_File_read_all(fh, &this->n_col, 1, MPI_INT, &status);
  MPI_File_read_at_all(fh, 6 * sizeof(int), &this->gnnz, 1, MPI_LONG_LONG, &status);
  MPI_File_read_at_all(fh, 10 * sizeof(int) + sizeof(long long int), &this->max_col, 1, MPI_INT, &status);

  size_t vec_size, depla_general, values_start, columns_start;
  depla_general = 11 * sizeof(int) + sizeof(long long int);

  this->pr = pr;
  this->pc = pc;
  this->NR = NR;
  this->NC = NC;

  this->ln_row = tbsla::utils::range::lnv(n_row, pr, NR);
  this->f_row = tbsla::utils::range::pflv(n_row, pr, NR);
  this->ln_col = tbsla::utils::range::lnv(n_col, pc, NC);
  this->f_col = tbsla::utils::range::pflv(n_col, pc, NC);

  MPI_File_read_at_all(fh, depla_general, &vec_size, 1, MPI_UNSIGNED_LONG, &status);
  depla_general += sizeof(size_t);
  this->nnz = 0;
  values_start = depla_general;
  depla_general += vec_size * sizeof(double);

  MPI_File_read_at_all(fh, depla_general, &vec_size, 1, MPI_UNSIGNED_LONG, &status);
  depla_general += sizeof(size_t);
  columns_start = depla_general;

  if(this->values)
    delete[] this->values;
  if(this->columns)
    delete[] this->columns;
  this->values = new double[this->ln_row * this->max_col];
  this->columns = new int[this->ln_row * this->max_col];

  long long int mod = this->ln_row % TBSLA_MATRIX_ELL_READLINES;
  tbsla::mpi::MatrixELL::mpiio_read_lines(fh, 0, mod, columns_start, values_start);
  for(long long int i = mod; i < this->ln_row; i += TBSLA_MATRIX_ELL_READLINES) {
    tbsla::mpi::MatrixELL::mpiio_read_lines(fh, i, TBSLA_MATRIX_ELL_READLINES, columns_start, values_start);
  }

  MPI_File_close(&fh);
  return 0;
}

void tbsla::mpi::MatrixELL::mpiio_read_lines(MPI_File &fh, long long int s, long long int n, long long int columns_start, long long int values_start) {
  MPI_Status status;
  std::vector<int> ctmp(n * this->max_col);
  std::vector<double> vtmp(n * this->max_col);
  MPI_File_read_at(fh, columns_start + (this->f_row + s) * this->max_col * sizeof(int), ctmp.data(), n * this->max_col, MPI_INT, &status);
  MPI_File_read_at(fh, values_start + (this->f_row + s) * this->max_col * sizeof(double), vtmp.data(), n * this->max_col, MPI_DOUBLE, &status);
  for(long long int i = 0; i < n; i++) {
    long long int incr = 0;
    for(long long int j = 0; j < this->max_col; j++) {
      long long int idx = ctmp[i * this->max_col + j];
      if(idx >= this->f_col && idx < this->f_col + this->ln_col) {
        this->columns[(i + s) * this->max_col + incr] = idx;
        this->values[(i + s) * this->max_col + incr] = vtmp[i * this->max_col + j];
        incr++;
      }
    }
    this->nnz += incr;
    for(long long int j = incr; j < this->max_col; j++) {
      this->columns[(i + s) * this->max_col + j] = 0;
      this->values[(i + s) * this->max_col + j] = 0;
    }
  }
}
