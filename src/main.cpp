#include <cmath>
#include <iostream>

#include "linalg/matrix.hpp"
#include "linalg/vector.hpp"

void matrix_tests() {
  using namespace linalg;

  matrix A{{
      {1, 2, 1, -4},
      {0, 0, 3, 0},
      {1, 2, 4, -4},
  }};

  std::cout << A << "\n";

  std::cout << "det A = " << A.cofactor_expansion() << "\n";
}

void vector_tests() {
  using namespace linalg;

  std::vector<vector> S = {
      vector{{1, 0, 1}},
      vector{{2, 0, 2}},
      vector{{1, 3, 4}},
      vector{{-4, 0, -4}},
  };

  std::vector<vector> redundancies = vector::redundant(S);

  for (auto &redundancy : redundancies)
    std::cout << redundancy << "\n";
}

int main() {
  vector_tests();

  return 0;
}
