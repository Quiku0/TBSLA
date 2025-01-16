#ifndef TBSLA_MPI_MatrixELL
#define TBSLA_MPI_MatrixELL

#include <tbsla/cpp/MatrixELL.hpp>
#include <tbsla/mpi/Matrix.hpp>
#include <iostream>

#include <mpi.h>

namespace tbsla { namespace mpi {

class MatrixELL : public tbsla::cpp::MatrixELL, public tbsla::mpi::Matrix {
  public:
    long long int read_bin_mpiio(MPI_Comm comm, std::string filename, long long int pr, long long int pc, long long int NR, long long int NC);
    void mpiio_read_lines(MPI_File &fh, long long int s, long long int n, long long int columns_start, long long int values_start);
    void fill_cdiag(MPI_Comm comm, long long int nr, long long int nc, long long int cdiag);
    void fill_cqmat(MPI_Comm comm, long long int n_row, long long int n_col, long long int c, double q, unsigned long long int seed_mult, long long int pr, long long int pc, long long int NR, long long int NC);
    using tbsla::cpp::MatrixELL::spmv;
    using tbsla::cpp::MatrixELL::Ax;
    using tbsla::cpp::MatrixELL::fill_cdiag;
    using tbsla::cpp::MatrixELL::fill_cqmat;
    using tbsla::cpp::MatrixELL::read;
    using tbsla::cpp::MatrixELL::write;
    using tbsla::mpi::Matrix::spmv_no_redist;
    using tbsla::mpi::Matrix::spmv;
    using tbsla::mpi::Matrix::Ax;
    using tbsla::mpi::Matrix::Ax_;
    using tbsla::mpi::Matrix::a_axpx_;
};

}}

#endif
