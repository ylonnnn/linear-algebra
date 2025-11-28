#pragma once

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <ostream>
#include <utility>

#include "types.hpp"
#include "vector.hpp"

namespace linear_algebra {
template <size_t M, size_t N> class matrix {
public:
  template <size_t _M = M, size_t _N = N>
  using container_type = std::array<std::array<real, _N>, _M>;

  matrix(bool augmented) : augmented_(augmented) {}
  matrix(container_type<> &&__container, bool augmented = false)
      : container_(std::move(__container)), augmented_(augmented) {}

  constexpr size_t n() const { return M * N; }

  real entry(size_t i, size_t j) const {
    return const_cast<matrix *>(this)->entry(i, j);
  }

  real &entry(size_t i, size_t j) {
    assert(i <= M);
    assert(j <= N);

    return container_[i - 1][j - 1];
  }

  matrix<N, M> transpose() {
    return matrix<N, M>(container_type(this->container_));
  }

  matrix<M, N + 1> augment(vector &&b) const {
    matrix<M, N + 1> augmented(true);

    for (size_t i = 0; i < M; ++i)
      for (size_t j = 0; j < N; ++j)
        augmented.entry(i + 1, j + 1) = entry(i + 1, j + 1);

    for (size_t j = 0; j < M; ++j)
      augmented.entry(j + 1, N + 1) = b[j];

    return augmented;
  }

  vector row_vec(size_t i) const {
    return vector(typename vector::container_type(container_[i - 1]));
  }

  vector col_vec(size_t j) const {
    typename vector::container_type underlying_arr(M);

    // Iterator over all the rows to retrieve the column within that row
    for (size_t i = 0; i < M; ++i)
      underlying_arr[i] = container_[i][j - 1];

    return vector(std::move(underlying_arr));
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
   * Returns the vector of pairs of pivots and their indices
   */
  std::map<size_t, real> gaussian() {
    size_t _N = N - augmented_;

    std::map<size_t, real> pivots;
    size_t p_idx = 0;

    // Row Echelon Form
    for (size_t i = 0; i < M; ++i) {
      size_t max_row = i;

      // Partial pivoting
      for (size_t j = i + 1; j < M; ++j)
        if (std::abs(container_[j][p_idx]) >
            std::abs(container_[max_row][p_idx]))
          max_row = j;

      if (max_row != i)
        elementary<1>(i + 1, max_row + 1);

      std::array<real, N> &row = container_[i];

      // Attempt to create the pivot on the current row
      real pivot = row[p_idx];
      for (; pivot == 0 && p_idx < _N; pivot = row[++p_idx])
        ;

      if (p_idx >= _N)
        continue;

      // assert(std::abs(pivot) >= std::numeric_limits<real>::epsilon());

      // Eliminate the nonzero numbers below the current pivot
      for (size_t j = i + 1; j < M; ++j) {
        std::array<real, N> &curr_row = container_[j]; // Below current pivot
        real b_pivot = curr_row[p_idx];

        if (b_pivot == 0)
          continue;

        elementary<3>(j + 1, i + 1, (1.0 / pivot) * -b_pivot);
      }

      pivots.try_emplace(p_idx, pivot);
      ++p_idx;
    }

    if (augmented_) {
      // Invalid Case: Zero Row = Non-Zero Number
      for (std::array<real, N> &a_row : container_) {
        vector row(vector::container_type(a_row.begin(), a_row.end() - 1));
        real b = a_row[N - 1];

        // Program termination due to system inconsisntency
        if (row.is_zero() && b != 0) {
          std::cout << "error: the system is inconsisntent\n";
          exit(0);
        }
      }
    }

    return pivots;
  }

  std::map<size_t, real> gauss_jordan() {
    std::map<size_t, real> piv_map = gaussian();
    std::vector<std::pair<size_t, real>> pivots(piv_map.begin(), piv_map.end());

    // Reduced Row Echelon Form
    for (size_t i = pivots.size(); i-- > 0;) {
      auto &[p_idx, _] = pivots[i];

      real pivot = container_[i][p_idx];
      elementary<2>(i + 1, 1.0 / pivot);

      piv_map[p_idx] = container_[i][p_idx];

      // Nonzero elimination above the pivot
      for (size_t j = i; j-- > 0;) {
        auto &row_above = container_[j];
        real above = row_above[p_idx];

        if (above == 0)
          continue;

        elementary<3>(j + 1, i + 1, -above);
      }
    }

    return piv_map;
  }

  std::vector<vector> solve_linear(vector &&b) {
    assert(b.n() == M);

    constexpr size_t _N = N + 1;
    matrix<M, _N> augmented = augment(std::move(b));

    std::map<size_t, real> pivots = augmented.gauss_jordan();
    size_t pn = pivots.size();

    std::cout << "A = " << augmented << "\n";

    std::vector<vector> solution;

    bool homogeneous = b.is_zero();
    size_t sol_n = (N - pn) + !homogeneous;

    solution.reserve(sol_n);

    for (size_t k = 0; k < (N + !homogeneous); ++k) {
      auto it = pivots.find(k);
      if (it != pivots.end())
        continue;

      vector sln(N);
      vector col = augmented.col_vec(k + 1);

      for (size_t i = 0; i < pn; ++i)
        sln[i] = -col[i];

      // Populate free variables
      sln[pn + solution.size()] = 1;

      solution.push_back(std::move(sln));
    }

    return solution;
  }

  std::vector<vector> column_space() {
    auto copy = matrix<M, N>(container_type<>(container_));
    std::map<size_t, real> pivots = copy.gauss_jordan();

    std::vector<vector> columns;
    columns.reserve(N);

    for (auto it = pivots.begin(); it != pivots.end(); ++it)
      columns.push_back(std::move(col_vec(it->first + 1)));

    return columns;
  }

  std::vector<vector> null_space() { return solve_linear(vector::zero<M>()); }

  real cofactor_expansion() {
    static_assert(M == N, "determinant requires a square matrix");

    if constexpr (M == 1)
      return entry(1, 1);

    // Edge case - 2 x 2
    else if constexpr (M == 2)
      return (entry(1, 1) * entry(2, 2)) - ((entry(1, 2) * entry(2, 1)));

    else {
      vector top = row_vec(1);
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

  vector transform(vector &&x) const {
    assert(N == x.n());

    vector y(M);

    for (size_t i = 0; i < N; ++i) {
      const std::array<real, N> &row = container_[i];
      vector row_vec(vector::container_type{row.begin(), row.end()});

      y[i] = row_vec.dot(x);
    }

    return y;
  }

  friend std::ostream &operator<<(std::ostream &os, const matrix<M, N> &mat) {
    os << "[";

    for (size_t i = 0; i < M; ++i) {
      os << "\n";
      for (size_t j = 0; j < N; ++j) {
        os << std::setw(8) << std::setprecision(4) << mat.entry(i + 1, j + 1);
        // << (j == N - 1 ? "" : " ");
      }
    }

    return os << "\n]";
  }

protected:
  container_type<> container_{};
  bool augmented_;
};

} // namespace linear_algebra
