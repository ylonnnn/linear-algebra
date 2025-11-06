#include <cmath>
#include <iostream>

#include "matrix.hpp"
#include "vector.hpp"

int main() {
  using namespace linear_algebra;

  matrix<3, 4> equation({{
      {1, 3, -2, 4},
      {-1, -3, 5, 2},
      {0, 0, 0, 0},
  }});

  equation.gauss_jordan();

  std::cout << equation << "\n";

  return 0;
}
