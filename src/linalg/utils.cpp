#include <cstdlib>
#include <limits>

#include "linalg/utils.hpp"

namespace linalg::utils {
bool approx_eq(real a, real b) {
  return std::abs(a - b) <= std::numeric_limits<real>::epsilon();
}

} // namespace linalg::utils
