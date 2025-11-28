#include <iomanip>

#include "linalg/matrix.hpp"
#include "linalg/vector.hpp"

namespace linalg {
vector::vector(size_t n) { container_ = container_type(n); }

vector::vector(typename vector_container::const_iterator first,
               typename vector_container::const_iterator last) {
  std::copy(first, last, container_.begin());
}

vector::vector(vector_container &&container)
    : container_(std::move(container)) {}

vector vector::zero(size_t N) { return container_type(N); }

std::vector<vector> vector::redundant(const std::vector<vector> &set) {
  size_t n = set.size();
  matrix mat(set[0].n(), n);

  for (size_t i = 0; i < n; ++i)
    mat.column(i + 1, vector::container_type{set[i].container_});

  std::map<size_t, real> pivots = mat.gauss_jordan();
  std::vector<vector> redundancies;

  size_t N = mat.N();
  redundancies.reserve(N);

  for (size_t k = 0; k < N; ++k) {
    auto it = pivots.find(k);
    if (it != pivots.end())
      continue;

    redundancies.push_back(set[k]);
  }

  return redundancies;
}

bool vector::is_orthogonal(const vector &vec) const { return dot(vec) == 0; }

bool vector::is_orthonormal(const vector &vec) const {
  return is_orthogonal(vec) && is_unit();
}

vector vector::projection(const vector &u) const {
  if (n() != u.n()) {
    std::cout << "error: dimension mismatch\n";
    exit(0);
  }

  return u.mult(dot(u) / u.dot(u));
}

size_t vector::n() const { return container_.size(); }

bool vector::is_zero() const {
  for (const real &x : container_)
    if (x != 0)
      return false;

  return true;
}

bool vector::is_unit() const { return norm() == 1; }

/**
 * Returns the angle in radians
 */
real vector::angle_between(const vector &other) const {
  real y = std::clamp(unit().dot(other.unit()), static_cast<real>(-1),
                      static_cast<real>(1));

  return std::acos(y);
}

vector vector::unit() const {
  vector res(*this);
  res.mult(1.0 / res.magnitude());

  return res;
}

void vector::add(const vector &other) {
  iterate([&](real &comp, size_t i) { comp += other[i]; });
}

vector vector::add(const vector &other) const {
  vector sum(*this);
  sum.add(other);

  return sum;
}

void vector::mult(real scalar) {
  iterate([&](real &comp, size_t) { comp *= scalar; });
}

vector vector::mult(real scalar) const {
  vector prod(*this);
  prod.mult(scalar);

  return prod;
}

real vector::dot(const vector &other) const {
  auto prod = 0.0;
  iterate([&](real comp, size_t i) { prod += comp * other[i]; });

  return prod;
}

real vector::magnitude() const { return std::sqrt(dot(*this)); }

real vector::norm() const { return magnitude(); }

void vector::normalize() { mult(1.0 / magnitude()); }

vector vector::normalize() const {
  vector normalized(*this);
  normalized.normalize();

  return normalized;
}

real &vector::operator[](size_t i) { return container_[i]; }
const real &vector::operator[](size_t i) const { return container_[i]; }

vector vector::operator+(const vector &b) const { return add(b); }

vector vector::operator-(const vector &b) const { return add(b * -1); }

vector vector::operator*(const vector &b) const { return dot(b); }
vector vector::operator*(real b) const { return mult(b); }

vector vector::operator/(real b) const { return mult(1.0 / b); }

std::ostream &operator<<(std::ostream &os, const vector &vec) {
  os << "[\n";

  vec.iterate([&](real comp, size_t i) { os << std::setw(8) << comp << "\n"; });

  return os << "]";
}

} // namespace linalg
