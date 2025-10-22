#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <ostream>

#include "types.hpp"

namespace linear_algebra {
template <size_t N> class vector {
  using vector_container = std::array<real, N>;

  vector_container container;

public:
  using container_type = vector_container;

  vector(typename vector_container::const_iterator first,
         typename vector_container::const_iterator last) {
    std::copy(first, last, container.begin());
  }

  vector(vector_container &&container) : container(std::move(container)) {}

  size_t n() const { return container.size(); }

  /**
   * Returns the angle in radians
   */
  real angle_between(const vector<N> &other) const {
    real y = std::clamp(unit().dot(other.unit()), static_cast<real>(-1),
                        static_cast<real>(1));

    return std::acos(y);
  }

  vector<N> unit() const {
    vector<N> res(*this);
    res.mult(1.0 / res.magnitude());

    return res;
  }

  void add(const vector<N> &other) {
    iterate([&](real &comp, size_t i) { comp += other[i]; });
  }

  vector<N> add(const vector<N> &other) const {
    vector<N> sum(*this);
    sum.add(other);

    return sum;
  }

  void mult(real scalar) {
    iterate([&](real &comp, size_t i) { comp *= scalar; });
  }

  vector<N> mult(real scalar) const {
    vector<N> prod(*this);
    prod.mult(scalar);

    return prod;
  }

  real dot(const vector<N> &other) const {
    auto prod = 0.0;
    iterate([&](real comp, size_t i) { prod += comp * other[i]; });

    return prod;
  }

  real magnitude() const {
    auto product = 0.0;
    iterate([&](real comp, size_t) { product += std::pow(comp, 2); });

    return std::sqrt(product);
  }

  real norm() const { return magnitude(); }

  real &operator[](size_t i) { return container[i]; }
  const real &operator[](size_t i) const { return container[i]; }

  friend std::ostream &operator<<(std::ostream &os, const vector<N> &vec) {
    os << "<";
    const char *sep = vec.n() > 4 ? ",\n" : ", ";

    vec.iterate([&](real comp, size_t i) {
      os << comp << (i == vec.n() - 1 ? "" : sep);
    });

    return os << ">";
  }

protected:
  using vector_iterator_fn = std::function<void(real &, size_t)>;
  using vector_const_iterator_fn = std::function<void(real, size_t)>;

  void iterate(const vector_iterator_fn &fn) {
    for (size_t i = 0; i < container.size(); ++i)
      fn(container[i], i);
  }

  void iterate(const vector_const_iterator_fn &fn) const {
    for (size_t i = 0; i < container.size(); ++i)
      fn(container[i], i);
  }
};

} // namespace linear_algebra
