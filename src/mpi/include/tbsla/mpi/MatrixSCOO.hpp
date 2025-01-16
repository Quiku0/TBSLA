#ifndef TBSLA_MPI_MatrixSCOO
#define TBSLA_MPI_MatrixSCOO

#include <tbsla/mpi/Matrix.hpp>
#include <tbsla/cpp/MatrixSCOO.hpp>
#include <iostream>
#include <fstream>

#include <mpi.h>

namespace tbsla { namespace mpi {

class MatrixSCOO : public tbsla::cpp::MatrixSCOO, public virtual tbsla::mpi::Matrix {
  public:
    long long int read_bin_mpiio(MPI_Comm comm, std::string filename, long long int pr, long long int pc, long long int NR, long long int NC);
    void fill_cdiag(MPI_Comm comm, long long int nr, long long int nc, long long int cdiag);
    void fill_cqmat(MPI_Comm comm, long long int n_row, long long int n_col, long long int c, double q, unsigned long long int seed_mult, long long int pr, long long int pc, long long int NR, long long int NC);
    using tbsla::cpp::MatrixSCOO::spmv;
    using tbsla::cpp::MatrixSCOO::Ax;
    using tbsla::cpp::MatrixSCOO::fill_cdiag;
    using tbsla::cpp::MatrixSCOO::fill_cqmat;
    using tbsla::cpp::MatrixSCOO::read;
    using tbsla::cpp::MatrixSCOO::write;
    using tbsla::mpi::Matrix::spmv_no_redist;
    using tbsla::mpi::Matrix::spmv;
    using tbsla::mpi::Matrix::Ax;
    using tbsla::mpi::Matrix::Ax_;
    using tbsla::mpi::Matrix::a_axpx_;
};

}}

#endif
