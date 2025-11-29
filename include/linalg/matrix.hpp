#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <map>
#include <ostream>
#include <utility>
#include <vector>

#include "linalg/types.hpp"
#include "linalg/vector.hpp"

namespace linalg {

class matrix {
public:
  using container_type = std::vector<std::vector<real>>;

  matrix(container_type &&__container, bool augmented = false);
  matrix(size_t M, size_t N, bool augmented = false);

  static matrix rotation(real theta);

  template <size_t N> static real angle(const matrix &mat) {
    std::cout << "error: no implementation for " << N << "-by-" << N
              << " matrices\n";
    exit(0);
  }

  void row(size_t i, std::vector<real> &&row);
  void row(size_t i, vector &&row);
  void column(size_t j, std::vector<real> &&col);
  void column(size_t j, vector &&col);

  size_t M() const;
  size_t N() const;

  real entry(size_t i, size_t j) const;
  real &entry(size_t i, size_t j);

  vector row_vec(size_t i) const;
  vector col_vec(size_t j) const;

  bool is_square() const;
  bool is_orthogonal() const;

  matrix transpose();
  matrix augment(vector &&b) const;

  void add(const matrix &other);
  matrix add(const matrix &other) const;

  void mult(real scalar);
  matrix mult(real scalar) const;

  // m x p * p x n = m x n
  template <size_t _N_> matrix mult(const matrix &other) const {
    container_type underlying_arr;

    for (size_t k = 0; k < M_; k++) {
      auto row = row_vec(k + 1);

      for (size_t l = 0; l < _N_; ++l) {
        auto col = other.col_vec(l + 1);
        underlying_arr[k][l] = row.dot(col);
      }
    }

    return matrix(std::move(underlying_arr));
  }

  /**
   * Type I - E_ij
   *
   * Swaps row i and row j
   */
  template <size_t Ty, typename = std::enable_if_t<Ty == 1>>
  void elementary(size_t i, size_t j) {
    std::swap(container_[i - 1], container_[j - 1]);
  }

  /**
   * Type II - E_i(c)
   *
   * M_ultiplies row i by the scalar c
   */
  template <size_t Ty, typename = std::enable_if_t<Ty == 2>>
  void elementary(size_t i, real scalar) {
    for (auto &comp : container_[i - 1]) {
      if (comp == 0)
        continue;

      comp *= scalar;
    }
  }

  /**
   * Type III - E_ij(c)
   *
   * M_ultiplies row j by scalar c and adds the product
   * to row i and replaces row i
   */
  template <size_t Ty, typename = std::enable_if_t<Ty == 3>>
  void elementary(size_t i, size_t j, real scalar) {
    for (size_t k = 0; k < N_; ++k)
      container_[i - 1][k] += scalar * container_[j - 1][k];
  }

  matrix invert();

  /**
   * Returns the vector of pairs of pivots and their indices
   */
  std::map<size_t, real> gaussian();
  std::map<size_t, real> gauss_jordan();

  std::vector<vector> solve_linear(vector &&b);

  std::vector<vector> column_space();
  std::vector<vector> null_space();

  real cofactor_expansion();
  vector transform(vector &&x) const;

  friend std::ostream &operator<<(std::ostream &os, const matrix &mat);

protected:
  size_t M_, N_;
  container_type container_{};

  bool augmented_ = false;
};

} // namespace linalg
