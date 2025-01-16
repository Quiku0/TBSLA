#include <tbsla/cpp/utils/vector.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>

long long int check(long long int i, double v, double exp, long long int return_value, bool debug) {
  if(debug)
    std::cout << "i = " << i << " r[i] = " << v << " exp " << exp << std::endl;
  if(v != exp) {
    return return_value;
  }
  return 0;
}

std::vector<double> compute_spmv_res_cdiag(long long int nr, long long int nc, long long int c, std::vector<double> & v) {
  std::vector<double> r(nr);
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

long long int tbsla::utils::vector::test_spmv_cdiag(long long int nr, long long int nc, long long int c, std::vector<double> & v, std::vector<double> & r, bool debug) {
  std::vector<double> rc = compute_spmv_res_cdiag(nr, nc, c, v);
  long long int i;
  for(i = 0; i < std::min(c, nr); i++) {
    long long int rv = check(i, r[i], rc[i], 10, debug);
    if(rv) return rv;
  }
  for(; i < std::min(nr, nc - c); i++) {
    long long int rv = check(i, r[i], rc[i], 12, debug);
    if(rv) return rv;
  }
  for(; i < nr; i++) {
    long long int rv = check(i, r[i], rc[i], 14, debug);
    if(rv) return rv;
  }
  if(i < nr) {
    return 222;
  }
  return 0;
}

long long int tbsla::utils::vector::test_a_axpx__cdiag(long long int nr, long long int nc, long long int c, std::vector<double> & v, std::vector<double> & r, bool debug) {
  if(nr != nc) {
    return 111;
  }
  std::vector<double> rc = compute_spmv_res_cdiag(nr, nc, c, v);
  std::transform (rc.begin(), rc.end(), v.begin(), rc.begin(), std::plus<double>());
  rc = compute_spmv_res_cdiag(nr, nc, c, rc);
  long long int i;
  for(i = 0; i < std::min(c, nr); i++) {
    long long int rv = check(i, r[i], rc[i], 10, debug);
    if(rv) return rv;
  }
  for(; i < std::min(nr, nc - c); i++) {
    long long int rv = check(i, r[i], rc[i], 12, debug);
    if(rv) return rv;
  }
  for(; i < nr; i++) {
    long long int rv = check(i, r[i], rc[i], 14, debug);
    if(rv) return rv;
  }
  if(i < nr) {
    return 222;
  }
  return 0;
}

void tbsla::utils::vector::print_dense_matrix(long long int nr, long long int nc, const std::vector<double> & m, std::ostream& os) {
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

long long int tbsla::utils::vector::compare_vectors(std::vector<double> v1, std::vector<double> v2) {
  if(v1.size() != v2.size()) {
    return -1;
  }
  long long int r = 0;
  for(long long int i = 0; i < v1.size(); i++) {
    long long int s = (v1[i] - v2[i]) * (v1[i] - v2[i]);
    if(r < s) r = s;
  }
  return r;
}
