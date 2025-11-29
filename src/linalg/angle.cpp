#include <cmath>

#include "linalg/angle.hpp"

namespace linalg::angle {
real deg_to_rad(real degree) { return degree * M_PI / 180; }

real rad_to_deg(real radian) { return radian * 180 / M_PI; }

} // namespace linalg::angle
