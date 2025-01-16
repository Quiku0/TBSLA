#ifndef TBSLA_MPI_MatrixCSR
#define TBSLA_MPI_MatrixCSR

#include <tbsla/cpp/MatrixCSR.hpp>
#include <tbsla/mpi/Matrix.hpp>
#include <iostream>

#include <mpi.h>

namespace tbsla { namespace mpi {

class MatrixCSR : public tbsla::cpp::MatrixCSR, public tbsla::mpi::Matrix {
  public:
    long long int read_bin_mpiio(MPI_Comm comm, std::string filename, long long int pr, long long int pc, long long int NR, long long int NC);
    void fill_cdiag(MPI_Comm comm, long long int nr, long long int nc, long long int cdiag);
    void fill_cqmat(MPI_Comm comm, long long int n_row, long long int n_col, long long int c, double q, unsigned long long int seed_mult, long long int pr, long long int pc, long long int NR, long long int NC);
    using tbsla::cpp::MatrixCSR::spmv;
    using tbsla::cpp::MatrixCSR::Ax;
    using tbsla::cpp::MatrixCSR::fill_cdiag;
    using tbsla::cpp::MatrixCSR::fill_cqmat;
    using tbsla::cpp::MatrixCSR::read;
    using tbsla::cpp::MatrixCSR::write;
    using tbsla::mpi::Matrix::spmv_no_redist;
    using tbsla::mpi::Matrix::spmv;
    using tbsla::mpi::Matrix::Ax;
    using tbsla::mpi::Matrix::Ax_;
    using tbsla::mpi::Matrix::a_axpx_;
private:
    void mpiio_read_lines(MPI_File &fh, long long int s, long long int n, long long int rowptr_start, long long int colidx_start, long long int values_start, size_t& mem_alloc);
};

}}

#endif
