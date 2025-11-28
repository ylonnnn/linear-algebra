#pragma once

#include "linalg/matrix.hpp"
#include "linalg/vector.hpp"

namespace linalg::basis {
template <size_t N> bool verify(const std::vector<vector> &basis) {
  size_t n = basis.size();

  // Dimension and column count must be the same
  if (n != N)
    return false;

  matrix Bmat(basis[0].n(), n);

  for (size_t i = 0; i < n; ++i)
    Bmat.column(i + 1, vector(basis[i]));

  return Bmat.gauss_jordan().size() == n;
}

bool is_orthogonal_basis(const std::vector<vector> &basis);
bool is_orthonormal_basis(const std::vector<vector> &basis);

} // namespace linalg::basis
