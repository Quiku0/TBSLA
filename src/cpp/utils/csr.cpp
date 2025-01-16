#include <tbsla/cpp/utils/csr.hpp>

bool tbsla::cpp::utils::csr::compare_row(const long long int* row, const long long int* col, unsigned long i, unsigned long j) {
  if (row[i] == row[j]) {
    return col[i] < col[j];
  }
  return row[i] < row[j];
}
