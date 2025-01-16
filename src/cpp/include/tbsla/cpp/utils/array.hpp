#ifndef TBSLA_CPP_ARRAY
#define TBSLA_CPP_ARRAY

#include <iostream>

namespace tbsla { namespace utils { namespace array {

template <class myType>
void stream(std::ostream &os, const std::string name, myType* v, long long int size) {
  os << name << " : ";
  if (size > 0) {
    for (long long int i = 0; i < size - 1; i++) {
      os << v[i] << ", ";
    }
    os << v[size - 1];
  }
}

long long int test_spmv_cdiag(long long int nr, long long int nc, long long int c, double* v, double* r, bool debug);
long long int test_a_axpx__cdiag(long long int nr, long long int nc, long long int c, double* v, double* r, bool debug);
void print_dense_matrix(long long int nr, long long int nc, const double* m, std::ostream& os);
long long int compare_arrays(double* v1, double* v2, long long int size);
long long int check(long long int i, double v, double exp, long long int return_value, bool debug);

}}}

#endif
