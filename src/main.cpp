#include <cmath>
#include <iostream>

#include "linalg/angle.hpp"
#include "linalg/matrix.hpp"
#include "linalg/vector.hpp"

void matrix_tests() {
  //
  using namespace linalg;

  matrix T = matrix::formula(
      {
          {{-1, 1, 1}},
          {{1, -1, 1}},
          {{1, 1, -1}},
      },
      {
          {{-2, 0}},
          {{-4, 2}},
          {{6, 2}},
      });

  std::cout << T << "\n";
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
  vector_tests();

  return 0;
}
