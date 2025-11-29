#include <cmath>

#include "linalg/angle.hpp"

namespace linalg::angle {
real to_radian(real degree) { return degree * M_PI / 180; }

real from_radian(real radian) { return radian * 180 / M_PI; }

} // namespace linalg::angle
