#ifndef TBSLA_CPP_UTILS_CSR
#define TBSLA_CPP_UTILS_CSR

namespace tbsla { namespace cpp { namespace utils { namespace csr {

template <typename T>
T* applyPermutation(const long long int* order, const T* t, long long int size) {
  T* st = new T[size];
  for(long long int i = 0; i < size; i++) {
    st[i] = t[order[i]];
  }
  return st;
}

bool compare_row(const long long int* row, const long long int* col, unsigned long i, unsigned long j);

}}}}

#endif
