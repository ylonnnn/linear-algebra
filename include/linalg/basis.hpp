#pragma once

#include "linalg/matrix.hpp"
#include "linalg/vector.hpp"

namespace linalg::basis {
template <size_t N> bool verify(const std::vector<vector> &set) {
  size_t n = set.size();

  // Dimension and column count must be the same
  if (n != N)
    return false;

  matrix Bmat(set[0].n(), n);

  for (size_t i = 0; i < n; ++i)
    Bmat.column(i + 1, vector(set[i]));

  return Bmat.gauss_jordan().size() == n;
}

template <size_t N> bool is_orthogonal_basis(const std::vector<vector> &set) {
  if (!verify<N>(set))
    return false;

  return vector::is_orthogonal_set(set);
}

template <size_t N> bool is_orthonormal_basis(const std::vector<vector> &set) {
  if (!verify<N>(set))
    return false;

  return vector::is_orthonormal_set(set);
}

} // namespace linalg::basis
