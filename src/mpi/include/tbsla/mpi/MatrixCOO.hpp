#ifndef TBSLA_MPI_MatrixCOO
#define TBSLA_MPI_MatrixCOO

#include <tbsla/mpi/Matrix.hpp>
#include <tbsla/cpp/MatrixCOO.hpp>
#include <iostream>
#include <fstream>

#include <mpi.h>

namespace tbsla { namespace mpi {

class MatrixCOO : public tbsla::cpp::MatrixCOO, public virtual tbsla::mpi::Matrix {
  public:
    long long int read_bin_mpiio(MPI_Comm comm, std::string filename, long long int pr, long long int pc, long long int NR, long long int NC);
    void fill_cdiag(MPI_Comm comm, long long int nr, long long int nc, long long int cdiag);
    void fill_cqmat(MPI_Comm comm, long long int n_row, long long int n_col, long long int c, double q, unsigned long long int seed_mult, long long int pr, long long int pc, long long int NR, long long int NC);
    double* spmv(MPI_Comm comm, const double* v, long long int vect_incr = 0);
    double* a_axpx_(MPI_Comm comm, const double* v, long long int vect_incr = 0);
    inline void Ax(MPI_Comm comm, double* r, const double* v, double* buffer, double* buffer2, long long int vect_incr = 0);
    double* page_rank(MPI_Comm comm, double beta, double epsilon, long long int max_iterations, long long int &nb_iterations_done);
    using tbsla::cpp::MatrixCOO::spmv;
    using tbsla::cpp::MatrixCOO::Ax;
    using tbsla::cpp::MatrixCOO::fill_cdiag;
    using tbsla::cpp::MatrixCOO::fill_cqmat;
    using tbsla::cpp::MatrixCOO::read;
    using tbsla::cpp::MatrixCOO::write;
    using tbsla::mpi::Matrix::spmv_no_redist;
};

}}

#endif
