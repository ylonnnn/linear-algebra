#pragma once

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <utility>

#include "types.hpp"
#include "vector.hpp"

namespace linear_algebra {
template <size_t M, size_t N> class matrix {
public:
  template <size_t _M = M, size_t _N = N>
  using container_type = std::array<std::array<real, _N>, _M>;
  using row_type = vector<N>;
  using column_type = vector<M>;

  matrix() = default;
  matrix(container_type<> &&__container) : container_(std::move(__container)) {}

  constexpr size_t n() const { return M * N; }

  real entry(size_t i, size_t j) const {
    assert(i <= M);
    assert(j <= N);

    return container_[i - 1][j - 1];
  }

  real &entry(size_t i, size_t j) {
    assert(i <= M);
    assert(j <= N);

    return container_[i - 1][j - 1];
  }

  matrix<N, M> transpose() {
    return matrix<N, M>(container_type(this->container_));
  }

  row_type row_vec(size_t i) const {
    return row_type(typename row_type::container_type(container_[i - 1]));
  }

  column_type col_vec(size_t j) const {
    typename column_type::container_type underlying_arr;

    // Iterator over all the rows to retrieve the column within that row
    for (size_t i = 0; i < M; ++i)
      underlying_arr[i] = container_[i][j - 1];

    return column_type(std::move(underlying_arr));
  }

  void add(const matrix<M, N> &other) {
    for (size_t i = 0; i < M; ++i)
      for (size_t j = 0; j < N; ++j)
        container_[i][j] += other.container_[i][j];
  }

  matrix<M, N> add(const matrix<M, N> &other) const {
    matrix<M, N> sum(container_type<>(this->container_));
    sum.add(other);

    return sum;
  }

  void mult(real scalar) {
    for (size_t i = 0; i < M; ++i)
      for (size_t j = 0; j < N; ++j)
        container_[i][j] *= scalar;
  }

  matrix<M, N> mult(real scalar) const {
    matrix<M, N> prod(container_type(this->container_));
    prod.mult(scalar);

    return prod;
  }

  // m x p * p x n = m x n
  template <size_t _N> matrix<M, _N> mult(const matrix<N, _N> &other) const {
    container_type<M, _N> underlying_arr;

    for (size_t k = 0; k < M; k++) {
      auto row = row_vec(k + 1);

      for (size_t l = 0; l < _N; ++l) {
        auto col = other.col_vec(l + 1);
        underlying_arr[k][l] = row.dot(col);
      }
    }

    return matrix<M, _N>(std::move(underlying_arr));
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
   * Multiplies row i by the scalar c
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
   * Multiplies row j by scalar c and adds the product
   * to row i and replaces row i
   */
  template <size_t Ty, typename = std::enable_if_t<Ty == 3>>
  void elementary(size_t i, size_t j, real scalar) {
    for (size_t k = 0; k < N; ++k)
      container_[i - 1][k] += scalar * container_[j - 1][k];
  }

  matrix<M, N> invert() {
    // Inversion can only be done on square matrices
    static_assert(M == N, "inversion can only be done on square matrices");

    // [A | In];
    matrix<M, N * 2> temp;

    // Copy the current contents to the left side (A)
    for (size_t i = 0; i < M; ++i)
      for (size_t j = 0; j < N; ++j) {
        real &entry = temp.entry(i + 1, j + 1);
        entry = container_[i][j];
      }

    // Initialize the identity matrix (In) on the right-side
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j) {
        real &entry = temp.entry(i + 1, N + j + 1);
        entry = i == j;
      }
    }

    // Inverse through Gauss-Jordan Elimination
    temp.gauss_jordan(); // Results to [In | A^-1]

    // Retrieve the right-side of the temp matrix
    matrix<M, N> inverted;

    // Retrieve the right-side (A^-1) of the matrix
    for (size_t i = 1; i <= M; ++i)
      for (size_t j = 1; j <= N; ++j) {
        real &entry = inverted.entry(i, j);
        entry = temp.entry(i, N + j);
      }

    return inverted;
  }

  /**
   * Returns the vector of pivot indices
   */
  std::vector<size_t> gaussian() {
    std::vector<size_t> p_idcs;
    p_idcs.reserve(N);

    size_t p_idx = 0;

    // Row Echelon Form
    for (size_t i = 0; i < M; ++i) {
      size_t max_row = i;

      // Partial pivoting
      for (size_t j = i + 1; j < M; ++j)
        if (std::abs(container_[j][p_idx]) >
            std::abs(container_[max_row][p_idx]))
          max_row = j;

      elementary<1>(i + 1, max_row + 1);

      std::array<real, N> &row = container_[i];

      // Attempt to create the pivot on the current row
      real pivot = row[p_idx];
      assert(std::abs(pivot) >= std::numeric_limits<real>::epsilon());

      if (pivot != 1)
        elementary<2>(i + 1, 1.0 / pivot);

      // Eliminate the nonzero numbers below the current pivot
      for (size_t j = i + 1; j < M; ++j) {
        std::array<real, N> &curr_row = container_[j]; // Below current pivot
        real b_pivot = curr_row[p_idx];

        if (b_pivot == 0)
          continue;

        elementary<3>(j + 1, i + 1, -b_pivot);
      }

      p_idcs.push_back(p_idx);
      ++p_idx;
    }

    return p_idcs;
  }

  void gauss_jordan() {
    std::vector<size_t> p_idcs = gaussian();

    // Reduced Row Echelon Form
    for (auto it = p_idcs.rbegin(); it != p_idcs.rend(); ++it) {
      size_t p_idx = *it;

      // If the pivot index is 0, there are no other rows above it
      if (p_idx == 0)
        continue;

      for (size_t i = p_idx; i-- > 0;) {
        std::array<real, N> &curr_row = container_[i];
        real a_pivot = curr_row[p_idx];

        if (a_pivot == 0)
          continue;

        elementary<3>(i + 1, p_idx + 1, -a_pivot);
      }
    }
  }

  real cofactor_expansion() {
    static_assert(M == N, "determinant requires a square matrix");

    if constexpr (M == 1)
      return entry(1, 1);

    // Edge case - 2 x 2
    else if constexpr (M == 2)
      return (entry(1, 1) * entry(2, 2)) - ((entry(1, 2) * entry(2, 1)));

    else {
      row_type top = row_vec(1);
      real determinant = 0;

      for (size_t r = 0; r < M; ++r) {
        container_type<M - 1, N - 1> submat_arr;
        real cofactor = ((r % 2 == 0) ? 1 : -1) * top[r];

        for (size_t i = 1; i < M; ++i) {
          size_t col = 0;
          for (size_t j = 0; j < N; ++j) {
            if (j == r)
              continue;

            submat_arr[i - 1][col] = container_[i][j];
            ++col;
          }
        }

        matrix<M - 1, N - 1> submatrix(std::move(submat_arr));
        determinant += cofactor * submatrix.cofactor_expansion();
      }

      return determinant;
    }
  }

  friend std::ostream &operator<<(std::ostream &os, const matrix<M, N> &mat) {
    os << "[";

    for (size_t i = 0; i < M; ++i) {
      os << "\n";
      for (size_t j = 0; j < N; ++j) {
        os << std::setw(4) << mat.entry(i + 1, j + 1)
           << (j == N - 1 ? "" : " ");
      }
    }

    return os << "\n]";
  }

protected:
  container_type<> container_;
};

} // namespace linear_algebra
