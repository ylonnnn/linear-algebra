#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iomanip>
#include <iostream>
#include <ostream>

#include "types.hpp"

namespace linear_algebra {
class vector {
  using vector_container = std::vector<real>;

public:
  using container_type = vector_container;

  vector(size_t n) { container_ = container_type(n); }

  vector(typename vector_container::const_iterator first,
         typename vector_container::const_iterator last) {
    std::copy(first, last, container_.begin());
  }

  vector(vector_container &&container) : container_(std::move(container)) {}

  template <size_t N> static vector zero() { return container_type(N); }

  template <size_t N> static vector standard(size_t k) {
    assert(k <= N);

    vector e = vector_container(N);
    e[k - 1] = 1;

    return e;
  }

  static bool is_orthogonal_set(const std::vector<vector> &set) {
    size_t n = set.size();

    for (size_t i = 0; i < n; ++i) {
      const vector &curr = set[i];

      for (size_t j = i + 1; j < n; ++j) {
        const vector &c_pair = set[j];
        if (!curr.is_orthogonal(c_pair))
          return false;
      }
    }

    return true;
  }

  static bool is_orthonormal_set(const std::vector<vector> &set) {
    size_t n = set.size();

    for (size_t i = 0; i < n; ++i) {
      const vector &curr = set[i];

      for (size_t j = i + 1; j < n; ++j) {
        const vector &c_pair = set[j];
        if (!curr.is_orthonormal(c_pair))
          return false;
      }
    }

    return true;
  }

  bool is_orthogonal(const vector &vec) const { return dot(vec) == 0; }

  bool is_orthonormal(const vector &vec) const {
    return is_orthogonal(vec) && is_unit();
  }

  size_t n() const { return container_.size(); }

  bool is_zero() const {
    for (const real &x : container_)
      if (x != 0)
        return false;

    return true;
  }

  bool is_unit() const { return norm() == 1; }

  /**
   * Returns the angle in radians
   */
  real angle_between(const vector &other) const {
    real y = std::clamp(unit().dot(other.unit()), static_cast<real>(-1),
                        static_cast<real>(1));

    return std::acos(y);
  }

  vector unit() const {
    vector res(*this);
    res.mult(1.0 / res.magnitude());

    return res;
  }

  void add(const vector &other) {
    iterate([&](real &comp, size_t i) { comp += other[i]; });
  }

  vector add(const vector &other) const {
    vector sum(*this);
    sum.add(other);

    return sum;
  }

  void mult(real scalar) {
    iterate([&](real &comp, size_t) { comp *= scalar; });
  }

  vector mult(real scalar) const {
    vector prod(*this);
    prod.mult(scalar);

    return prod;
  }

  real dot(const vector &other) const {
    auto prod = 0.0;
    iterate([&](real comp, size_t i) { prod += comp * other[i]; });

    return prod;
  }

  real magnitude() const { return std::sqrt(dot(*this)); }

  real norm() const { return magnitude(); }

  void normalize() { mult(1.0 / magnitude()); }

  vector normalize() const {
    vector normalized(*this);
    normalized.normalize();

    return normalized;
  }

  real &operator[](size_t i) { return container_[i]; }
  const real &operator[](size_t i) const { return container_[i]; }

  friend std::ostream &operator<<(std::ostream &os, const vector &vec) {
    os << "[\n";

    vec.iterate(
        [&](real comp, size_t i) { os << std::setw(8) << comp << "\n"; });

    return os << "]";
  }

protected:
  using vector_iterator_fn = std::function<void(real &, size_t)>;
  using vector_const_iterator_fn = std::function<void(real, size_t)>;

  void iterate(const vector_iterator_fn &fn) {
    for (size_t i = 0; i < container_.size(); ++i)
      fn(container_[i], i);
  }

  void iterate(const vector_const_iterator_fn &fn) const {
    for (size_t i = 0; i < container_.size(); ++i)
      fn(container_[i], i);
  }

private:
  vector_container container_{};
};

} // namespace linear_algebra
