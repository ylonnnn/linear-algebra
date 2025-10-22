#include <cmath>
#include <iostream>

#include "vector.hpp"

int main() {
  using namespace linear_algebra;

  vector<3> x{{1, 2, 3}};
  vector<3> y{{4, 5, 6}};

  std::cout << x << "\n";
  std::cout << y << "\n";

  real angle = x.angle_between(y);

  std::cout << "degrees: " << angle * (180 / M_PI) << "\n";
  std::cout << "radians: " << angle << "\n";

  return 0;
}
