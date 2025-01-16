#include <tbsla/cpp/MatrixCOO.hpp>
#include <tbsla/cpp/MatrixSCOO.hpp>
#include <tbsla/cpp/MatrixCSR.hpp>
#include <tbsla/cpp/MatrixELL.hpp>
#include <tbsla/cpp/MatrixDENSE.hpp>
#include <tbsla/cpp/Matrix.hpp>
#include <tbsla/cpp/utils/array.hpp>

#include <iostream>

void print(tbsla::cpp::Matrix & m) {
  std::cout << "--------" << std::endl;
  std::cout << m << std::endl;
  std::cout << "--------" << std::endl;
}

void test_matrix(tbsla::cpp::Matrix & m, long long int c) {
  long long int nc, nr;
  nc = m.get_n_col();
  nr = m.get_n_row();
  double* v = new double[nc]();
  for(long long int i = 0; i < nc; i++) {
    v[i] = 2 * i + 1;
  }
  double* r = m.a_axpx_(v);
  long long int res;
  res = tbsla::utils::array::test_a_axpx__cdiag(nr, nc, c, v, r, false);
  if (res) {
    tbsla::utils::array::stream<double>(std::cout, "v", v, nc);
    std::cout << std::endl;
    tbsla::utils::array::stream<double>(std::cout, "r", r, nr);
    std::cout << std::endl;
    print(m);
    res = tbsla::utils::array::test_a_axpx__cdiag(nr, nc, c, v, r, true);
    exit(res);
  }
  delete[] v;
  delete[] r;
}

void test_cdiag(long long int nr, long long int nc, long long int c) {
  std::cout << "---- nr : " << nr << "; nc : " << nc << "; c : " << c << " ----  " << std::endl;
  tbsla::cpp::MatrixCOO mcoo;
  mcoo.fill_cdiag(nr, nc, c);
  test_matrix(mcoo, c);

  tbsla::cpp::MatrixSCOO mscoo(mcoo);
  test_matrix(mscoo, c);

  tbsla::cpp::MatrixSCOO mscoo2;
  mscoo2.fill_cdiag(nr, nc, c);
  test_matrix(mscoo2, c);

  tbsla::cpp::MatrixCSR mcsr(mcoo);
  test_matrix(mcsr, c);

  tbsla::cpp::MatrixCSR mcsr2;
  mcsr2.fill_cdiag(nr, nc, c);
  test_matrix(mcsr2, c);

  tbsla::cpp::MatrixELL mell(mcoo);
  test_matrix(mell, c);

  tbsla::cpp::MatrixELL mell2;
  mell2.fill_cdiag(nr, nc, c);
  test_matrix(mell2, c);

  tbsla::cpp::MatrixDENSE mdense(mcoo);
  test_matrix(mdense, c);

  tbsla::cpp::MatrixDENSE mdense2;
  mdense2.fill_cdiag(nr, nc, c);
  test_matrix(mdense2, c);
}

int main(int argc, char** argv) {

  long long int t = 0;
  for(long long int i = 0; i <= 12; i++) {
    std::cout << "=== test " << t++ << " ===" << std::endl;
    test_cdiag(10, 10, i);
  }
  for(long long int i = 0; i <= 12; i++) {
    std::cout << "=== test " << t++ << " ===" << std::endl;
    test_cdiag(30, 30, 2 * i);
  }

  std::cout << "=== finished without error === " << std::endl;
}
