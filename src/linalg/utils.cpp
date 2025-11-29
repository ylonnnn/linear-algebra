#include <cmath>
#include <limits>

#include "linalg/utils.hpp"

namespace linalg::utils {
real sign(real n) { return n < 0 ? -1 : 1; }

bool approx_eq(real a, real b) {
  return std::abs(a - b) <= std::numeric_limits<real>::epsilon();
}

real cos(real x) {
  real cosine = std::cos(x);
  return std::abs(cosine) <= std::numeric_limits<real>::epsilon() ? 0 : cosine;
}

real sin(real x) {
  real sine = std::sin(x);
  return std::abs(sine) <= std::numeric_limits<real>::epsilon() ? 0 : sine;
}

} // namespace linalg::utils
