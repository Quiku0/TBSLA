#include <tbsla/cpp/utils/array.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>

long long int tbsla::utils::array::check(long long int i, double v, double exp, long long int return_value, bool debug) {
  if(debug)
    std::cout << "i = " << i << " r[i] = " << v << " exp " << exp << std::endl;
  if(v != exp) {
    std::cout << "error --> i = " << i << " r[i] = " << v << " exp " << exp << std::endl;
    return return_value;
  }
  return 0;
}

double* compute_spmv_res_cdiag(long long int nr, long long int nc, long long int c, double* v) {
  double* r = new double[nr]();
  long long int i;
  for(i = 0; i < std::min(c, nr); i++) {
    r[i] = i < nc - c ? v[i + c] : 0;
  }
  for(; i < std::min(nr, nc - c); i++) {
    r[i] = c == 0 ? v[i] : v[i + c] + v[i - c];
  }
  for(; i < nr; i++) {
    r[i] = i < nc + c ? v[i - c] : 0;
  }
  return r;
}

long long int tbsla::utils::array::test_spmv_cdiag(long long int nr, long long int nc, long long int c, double* v, double* r, bool debug) {
  double* rc = compute_spmv_res_cdiag(nr, nc, c, v);
  long long int i;
  for(i = 0; i < std::min(c, nr); i++) {
    long long int rv = tbsla::utils::array::check(i, r[i], rc[i], 10, debug);
    if(rv) return rv;
  }
  for(; i < std::min(nr, nc - c); i++) {
    long long int rv = tbsla::utils::array::check(i, r[i], rc[i], 12, debug);
    if(rv) return rv;
  }
  for(; i < nr; i++) {
    long long int rv = tbsla::utils::array::check(i, r[i], rc[i], 14, debug);
    if(rv) return rv;
  }
  delete [] rc;
  if(i < nr) {
    return 222;
  }
  return 0;
}

long long int tbsla::utils::array::test_a_axpx__cdiag(long long int nr, long long int nc, long long int c, double* v, double* r, bool debug) {
  if(nr != nc) {
    return 111;
  }
  double* rc = compute_spmv_res_cdiag(nr, nc, c, v);
  for (long long int i = 0; i < nr; i++) {
    rc[i] += v[i];
  }
  double * tmp = rc;
  rc = compute_spmv_res_cdiag(nr, nc, c, rc);
  delete[] tmp;
  long long int i;
  for(i = 0; i < std::min(c, nr); i++) {
    long long int rv = tbsla::utils::array::check(i, r[i], rc[i], 10, debug);
    if(rv) return rv;
  }
  for(; i < std::min(nr, nc - c); i++) {
    long long int rv = tbsla::utils::array::check(i, r[i], rc[i], 12, debug);
    if(rv) return rv;
  }
  for(; i < nr; i++) {
    long long int rv = tbsla::utils::array::check(i, r[i], rc[i], 14, debug);
    if(rv) return rv;
  }
  delete [] rc;
  if(i < nr) {
    return 222;
  }
  return 0;
}

void tbsla::utils::array::print_dense_matrix(long long int nr, long long int nc, const double* m, std::ostream& os) {
  if (m == NULL) {
    return;
  }
  long long int i, j, max_p = 1, log_r;
  for (i = 0; i < nr * nc; i++) {
    log_r = log10(fabs(m[i])) + 1;
    if (log_r > max_p) {
      max_p = log_r;
    }
  }
  os << std::fixed;
  os.precision(max_p);
  for(i = 0; i < nr; i++) {
    for(j = 0; j < nc; j++) {
      os << std::setw(max_p + 3) << m[i * nc + j] << "  ";
    }
    os << std::endl;
  }
}

long long int tbsla::utils::array::compare_arrays(double* v1, double* v2, long long int size) {
  long long int r = 0;
  for(long long int i = 0; i < size; i++) {
    long long int s = (v1[i] - v2[i]) * (v1[i] - v2[i]);
    if(r < s) r = s;
  }
  return r;
}
