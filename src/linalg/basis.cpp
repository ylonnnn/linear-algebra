#include "linalg/basis.hpp"

namespace linalg::basis {
bool is_orthogonal_basis(const std::vector<vector> &basis) {
  size_t n = basis.size();

  for (size_t i = 0; i < n; ++i) {
    const vector &curr = basis[i];

    for (size_t j = i + 1; j < n; ++j) {
      const vector &c_pair = basis[j];
      if (!curr.is_orthogonal(c_pair))
        return false;
    }
  }

  return true;
}

bool is_orthonormal_basis(const std::vector<vector> &basis) {
  size_t n = basis.size();

  for (size_t i = 0; i < n; ++i) {
    const vector &curr = basis[i];

    for (size_t j = i + 1; j < n; ++j) {
      const vector &c_pair = basis[j];
      if (!curr.is_orthonormal(c_pair))
        return false;
    }
  }

  return true;
}

} // namespace linalg::basis
