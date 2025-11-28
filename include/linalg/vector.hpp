#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <ostream>

#include "linalg/types.hpp"

namespace linalg {
class matrix;

class vector {
  using vector_container = std::vector<real>;

public:
  using container_type = vector_container;

  vector(size_t n);
  vector(typename vector_container::const_iterator first,
         typename vector_container::const_iterator last);

  vector(vector_container &&container);

  static vector zero(size_t N);
  template <size_t N> static vector standard(size_t k) {
    assert(k <= N);

    vector e = vector_container(N);
    e[k - 1] = 1;

    return e;
  }

  static std::vector<vector> redundant(const std::vector<vector> &set);

  static bool is_orthogonal_set(const std::vector<vector> &set);
  static bool is_orthonormal_set(const std::vector<vector> &set);

  bool is_orthogonal(const vector &vec) const;
  bool is_orthonormal(const vector &vec) const;

  size_t n() const;

  bool is_zero() const;
  bool is_unit() const;

  /**
   * Returns the angle in radians
   */
  real angle_between(const vector &other) const;

  vector unit() const;
  void add(const vector &other);
  vector add(const vector &other) const;

  void mult(real scalar);
  vector mult(real scalar) const;

  real dot(const vector &other) const;
  real magnitude() const;

  real norm() const;
  void normalize();

  vector normalize() const;

  real &operator[](size_t i);
  const real &operator[](size_t i) const;

  friend std::ostream &operator<<(std::ostream &os, const vector &vec);

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

  friend matrix;
};

} // namespace linalg
