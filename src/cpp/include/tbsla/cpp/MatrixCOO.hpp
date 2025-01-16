#ifndef TBSLA_CPP_MatrixCOO
#define TBSLA_CPP_MatrixCOO

#include <tbsla/cpp/Matrix.hpp>
#include <iostream>
#include <fstream>

namespace tbsla { namespace cpp {

class MatrixCOO : public virtual Matrix {
  public:
    friend std::ostream & operator<<( std::ostream &os, const MatrixCOO &m);
    MatrixCOO(long long int n_row, long long int n_col, double* values, long long int* row, long long int* col);
    MatrixCOO(long long int n_row, long long int n_col, long long int n_values);
    MatrixCOO(long long int n_row, long long int n_col);
    MatrixCOO() : values(0), row(0), col(0) {};
    ~MatrixCOO();
    const double* get_values() const { return values; }
    const long long int* get_row() const { return row; }
    const long long int* get_col() const { return col; }
    double* spmv(const double* v, long long int vect_incr = 0) const;
    inline void Ax(double* r, const double* v, long long int vect_incr = 0) const;
    using tbsla::cpp::Matrix::a_axpx_;
    using tbsla::cpp::Matrix::AAxpAx;
    using tbsla::cpp::Matrix::AAxpAxpx;
    std::ostream & print_as_dense(std::ostream &os);
    std::ostream & write(std::ostream &os);
    std::istream & read(std::istream &is, std::size_t pos = 0, std::size_t n = 1);
    std::ostream& print(std::ostream& os) const;
    void NUMAinit();

    void readMM(std::string name);
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

  private:
    void init(long long int n_row, long long int n_col, long long int n_values);
};

}}

#endif
