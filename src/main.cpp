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
}

int main() {
  matrix_tests();

  return 0;
}
