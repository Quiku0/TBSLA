#ifndef TBSLA_CPP_MatrixSCOO
#define TBSLA_CPP_MatrixSCOO

#include <tbsla/cpp/Matrix.hpp>
#include <tbsla/cpp/MatrixCOO.hpp>
#include <iostream>
#include <fstream>

namespace tbsla { namespace cpp {

class MatrixSCOO : public virtual Matrix {
  public:
    friend std::ostream & operator<<( std::ostream &os, const MatrixSCOO &m);
    MatrixSCOO(long long int n_row, long long int n_col, double* values, long long int* row,long long int* col);
    MatrixSCOO(long long int n_row, long long int n_col, long long int n_values);
    MatrixSCOO(long long int n_row, long long int n_col);
    MatrixSCOO() : values(0), row(0), col(0) {};
    MatrixSCOO(const tbsla::cpp::MatrixCOO & m);
    ~MatrixSCOO();
    double* spmv(const double* v, long long int vect_incr = 0) const;
    inline void Ax(double* r, const double* v, long long int vect_incr = 0) const;
    using tbsla::cpp::Matrix::a_axpx_;
    using tbsla::cpp::Matrix::AAxpAx;
    using tbsla::cpp::Matrix::AAxpAxpx;
    void push_back(long long int r, long long int c, double v);
    std::ostream & print_as_dense(std::ostream &os);
    std::ostream & write(std::ostream &os);
    std::istream & read(std::istream &is, std::size_t pos = 0, std::size_t n = 1);
    std::ostream& print(std::ostream& os) const;
    void NUMAinit();

    void fill_cdiag(long long int n_row, long long int n_col, long long int cdiag, long long int pr = 0, long long int pc = 0, long long int NR = 1, long long int NC = 1);
    void fill_cqmat(long long int n_row, long long int n_col, long long int c, double q, unsigned long long int seed_mult = 1, long long int pr = 0, long long int pc = 0, long long int NR = 1, long long int NC = 1);
    void fill_random(long long int n_row, long long int n_col, double nnz_ratio, unsigned long long int seed_mult = 1, long long int pr = 0, long long int pc = 0, long long int NR = 1, long long int NC = 1);
    void fill_cdistrib(long long int n_row, long long int n_col, long long int nnz, long long int pr = 0, long long int pc = 0, long long int NR = 1, long long int NC = 1);
    void fill_brain(long long int n_row, long long int n_col, int* neuron_type, std::vector<std::vector<double> > proba_conn, std::vector<std::unordered_map<int,std::vector<int> > > brain_struct, unsigned long long int seed_mult, long long int pr = 0, long long int pc = 0, long long int NR = 1, long long int NC = 1);
	
	void get_row_sums(double* buffer);
	void normalize_rows(double* buffer);
	void get_col_sums(double* buffer);
	void normalize_cols(double* buffer);
    void set_diag(double* s);

  protected:
    double* values;
    long long int* row;
    long long int* col;
};

}}

#endif
