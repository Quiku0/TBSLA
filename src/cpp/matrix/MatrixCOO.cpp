#include <tbsla/cpp/MatrixCOO.hpp>
#include <tbsla/cpp/utils/vector.hpp>
#include <iostream>
#include <vector>
#include <string>

MatrixCOO::MatrixCOO(int n_row, int n_col, std::vector<double> & values, std::vector<int> & row,  std::vector<int> & col) {
  this->n_row = n_row;
  this->n_col = n_col;
  this->values = values;
  this->row = row;
  this->col = col;
}

MatrixCOO::MatrixCOO(int n_row, int n_col, int n_values) {
  this->n_row = n_row;
  this->n_col = n_col;
  this->values.reserve(n_values);
  this->row.reserve(n_values);
  this->col.reserve(n_values);
}

MatrixCOO::MatrixCOO(int n_row, int n_col) {
  this->n_row = n_row;
  this->n_col = n_col;
}

std::ostream & operator<<( std::ostream &os, const MatrixCOO &m) {
  os << "n_row : " << m.n_row << std::endl;
  os << "n_col : " << m.n_col << std::endl;
  os << "n_values : " << m.values.size() << std::endl;
  tbsla::utils::vector::streamvector<double>(os, "val", m.values);
  os << std::endl;
  tbsla::utils::vector::streamvector<int>(os, "row", m.row);
  os << std::endl;
  tbsla::utils::vector::streamvector<int>(os, "col", m.col);
  return os;
}

std::vector<double> MatrixCOO::spmv(const std::vector<double> &v) {
  std::vector<double> r (this->n_row, 0);
  for (int i = 0; i < this->values.size(); i++) {
     r[this->row[i]] += this->values[i] * v[this->col[i]];
  }
  return r;
}

void MatrixCOO::push_back(int r, int c, double v) {
  this->values.push_back(v);
  this->row.push_back(r);
  this->col.push_back(c);
}
