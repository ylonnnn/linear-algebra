#include <cmath>
#include <iostream>

#include "linalg/angle.hpp"
#include "linalg/matrix.hpp"
#include "linalg/vector.hpp"

void matrix_tests() {
  //
  using namespace linalg;

  matrix A{{
      {1, 1, -1},
      {-1, -1, -1},
      {1, -1, -1},
  }};

  std::cout << A.cofactor_expansion() << "\n";
  std::cout << A.gaussian_diagonal_det() << "\n";
}

void vector_tests() {
  //
  using namespace linalg;

  std::vector<vector> set{vector{{1, 2, 2}}, {{2, 1, 0}}, {{0, 1, 2}}};
  std::vector<vector> orthogonal_set = vector::gram_schmidt(set);
  std::vector<vector> orthonormal_set = vector::gram_schmidt(set, true);

  for (auto &v : orthogonal_set)
    std::cout << v << "\n";

  for (auto &v : orthonormal_set)
    std::cout << v << "\n";
}

int main() {
  matrix_tests();

  return 0;
}
