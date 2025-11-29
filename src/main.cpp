#include <cmath>
#include <iostream>

#include "linalg/angle.hpp"
#include "linalg/matrix.hpp"
#include "linalg/vector.hpp"

void matrix_tests() {
  using namespace linalg;

  auto R_ = matrix{{
      {-sqrt(3) / 2, -1.0 / 2},
      {1.0 / 2, -sqrt(3) / 2},
  }};

  std::cout << std::boolalpha;

  std::cout << R_.is_square() << "\n";
  std::cout << R_.is_orthogonal() << "\n";

  std::cout << "Angle of Rotation (Degrees): "
            << angle::rad_to_deg(matrix::angle<2>(R_)) << "\n";
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
  matrix_tests();

  return 0;
}
