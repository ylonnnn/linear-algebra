#include <cmath>
#include <iostream>

#include "matrix.hpp"
#include "vector.hpp"

void matrix_tests() {
  using namespace linear_algebra;

  matrix<2, 3> A({{
      {2, 3, -1},
      {1, 0, 2},
  }});

  std::cout << A.transform(vector::standard<3>(1)) << "\n";
  std::cout << A.transform(vector::standard<3>(2)) << "\n";
  std::cout << A.transform(vector::standard<3>(3)) << "\n";
}

void vector_tests() {
  using namespace linear_algebra;

  vector a{vector::container_type{3, 0}}, b{vector::container_type{0, 2}};
  vector x{vector::container_type{1, 0}}, y{vector::container_type{0, 1}};

  std::cout << "ab: " << a.dot(b) << "\n";
  std::cout << "xy: " << x.dot(y) << "\n";

  std::cout << "is_orthogonal: " << std::boolalpha << a.is_orthogonal(b)
            << "\n";
  std::cout << "is_orthonormal: " << std::boolalpha << a.is_orthonormal(b)
            << "\n";

  std::cout << "is_orthogonal: " << std::boolalpha << x.is_orthogonal(y)
            << "\n";
  std::cout << "is_orthonormal: " << std::boolalpha << x.is_orthonormal(y)
            << "\n";
}

int main() {
  matrix_tests();

  return 0;
}
