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

  vector u = {{1, 1}}, v = {{3, 2}}, proju_v = v.projection(u);

  std::cout << u << "\n";
  std::cout << v << "\n";

  std::cout << proju_v << "\n";

  std::cout << (v - proju_v) << "\n";

  std::cout << u.is_orthogonal(v - proju_v) << "\n";
  std::cout << u.is_orthogonal(proju_v) << "\n";
}

int main() {
  vector_tests();

  return 0;
}
