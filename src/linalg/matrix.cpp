#include <iomanip>
#include <limits>

#include "linalg/matrix.hpp"
#include "linalg/utils.hpp"
#include "linalg/vector.hpp"

namespace linalg {

matrix::matrix(container_type &&__container, bool augmented)
    : M_(__container.size()), N_(__container[0].size() + augmented),
      container_(std::move(__container)), augmented_(augmented) {}

matrix::matrix(size_t M, size_t N, bool augmented)
    : M_(M), N_(N + augmented), container_({}), augmented_(augmented) {
  container_.resize(M_);

  for (auto &row : container_)
    row.resize(N_);
}

template <> matrix matrix::rotation<2>(real theta) {
  real cosine = utils::cos(theta), sine = utils::sin(theta);

  return matrix({
      {cosine, -sine}, // [cos, -sin]
      {sine, cosine},  // [sin,  cos]
  });
}

template <> real matrix::angle<2>(const matrix &mat) {
  assert(mat.is_orthogonal());
  return utils::sign(mat.entry(2, 1)) * acos(mat.entry(1, 1));
}

matrix matrix::formula(const std::vector<vector> &basis,
                       const std::vector<vector> &output) {
  size_t n = basis.size(), m = output.size();
  matrix b_mat{basis[0].n(), n}, o_mat{output[0].n(), m};

  b_mat.columns(basis);
  o_mat.columns(output);

  return o_mat.mult(b_mat.pseudo_inverse());
}

void matrix::row(size_t i, std::vector<real> &&row) {
  assert(i <= M_);

  container_[i - 1] = std::move(row);
}

void matrix::row(size_t i, vector &&row) {
  this->row(i, std::move(row.container_));
}

void matrix::column(size_t j, std::vector<real> &&col) {
  assert(j <= N_);

  for (size_t i = 0; i < M_; ++i)
    container_[i][j - 1] = col[i];
}

void matrix::column(size_t j, vector &&col) {
  column(j, std::move(col.container_));
}

void matrix::rows(const std::vector<vector> &rows) {
  for (size_t i = 0; i < rows.size(); ++i)
    row(i + 1, vector::container_type(rows[i].container_));
}

void matrix::columns(const std::vector<vector> &cols) {
  for (size_t j = 0; j < cols.size(); ++j)
    column(j + 1, vector::container_type(cols[j].container_));
}

size_t matrix::M() const { return M_; }

size_t matrix::N() const { return N_; }

real matrix::entry(size_t i, size_t j) const {
  return const_cast<matrix *>(this)->entry(i, j);
}

real &matrix::entry(size_t i, size_t j) {
  assert(i <= M_);
  assert(j <= N_);

  return container_[i - 1][j - 1];
}

vector matrix::row_vec(size_t i) const {
  return vector(typename vector::container_type(container_[i - 1]));
}

vector matrix::col_vec(size_t j) const {
  typename vector::container_type underlying_arr(M_);

  // Iterator over all the rows to retrieve the column within that row
  for (size_t i = 0; i < M_; ++i)
    underlying_arr[i] = container_[i][j - 1];

  return vector(std::move(underlying_arr));
}

bool matrix::is_square() const { return M_ == N_; }

bool matrix::is_orthogonal() const {
  if (!is_square()) // Orthogonal matrices must always be square
    return false;

  std::vector<vector> set;
  set.reserve(M_);

  for (auto &row : container_)
    set.emplace_back(vector::container_type(row));

  return vector::is_orthonormal_set(set);
}

matrix matrix::transpose() const {
  matrix mat(N_, M_);

  for (size_t i = 1; i <= M_; ++i)
    for (size_t j = 1; j <= N_; ++j)
      mat.entry(j, i) = entry(i, j);

  return mat;
}

matrix matrix::T() const { return transpose(); }

matrix matrix::pseudo_inverse() const {
  matrix _T = T();
  return (_T * *this).invert() * _T;
}

matrix matrix::augment(vector &&b) const {
  matrix augmented(M_, N_, true);

  for (size_t i = 0; i < M_; ++i)
    for (size_t j = 0; j < N_; ++j)
      augmented.entry(i + 1, j + 1) = entry(i + 1, j + 1);

  for (size_t j = 0; j < M_; ++j)
    augmented.entry(j + 1, N_ + 1) = b[j];

  return augmented;
}

void matrix::add(const matrix &other) {
  for (size_t i = 0; i < M_; ++i)
    for (size_t j = 0; j < N_; ++j)
      container_[i][j] += other.container_[i][j];
}

matrix matrix::add(const matrix &other) const {
  matrix sum(container_type(this->container_));
  sum.add(other);

  return sum;
}

void matrix::mult(real scalar) {
  for (size_t i = 0; i < M_; ++i)
    for (size_t j = 0; j < N_; ++j)
      container_[i][j] *= scalar;
}

matrix matrix::mult(real scalar) const {
  matrix prod(container_type(this->container_));
  prod.mult(scalar);

  return prod;
}
matrix matrix::mult(const matrix &other) const {
  matrix mat(M_, other.N_);

  for (size_t k = 0; k < M_; k++) {
    auto row = row_vec(k + 1);

    for (size_t l = 0; l < other.N_; ++l) {
      auto col = other.col_vec(l + 1);
      mat.entry(k + 1, l + 1) = row.dot(col);
    }
  }

  return mat;
}

matrix matrix::invert() {
  // Inversion can only be done on square matrices
  assert(is_square());

  // [A | In];
  matrix temp(M_, N_ * 2);

  // Copy the current contents to the left side (A)
  for (size_t i = 0; i < M_; ++i)
    for (size_t j = 0; j < N_; ++j) {
      real &entry = temp.entry(i + 1, j + 1);
      entry = container_[i][j];
    }

  // Initialize the identity matrix (In) on the right-side
  for (size_t i = 0; i < N_; ++i) {
    for (size_t j = 0; j < N_; ++j) {
      real &entry = temp.entry(i + 1, N_ + j + 1);
      entry = i == j;
    }
  }

  // Inverse through Gauss-Jordan Elimination
  temp.gauss_jordan(); // Results to [In | A^-1]

  // Retrieve the right-side of the temp matrix
  matrix inverted(M_, N_);

  // Retrieve the right-side (A^-1) of the matrix
  for (size_t i = 1; i <= M_; ++i)
    for (size_t j = 1; j <= N_; ++j) {
      real &entry = inverted.entry(i, j);
      entry = temp.entry(i, N_ + j);
    }

  return inverted;
}

/**
 * Returns the vector of pairs of pivots and their indices
 */
std::map<size_t, real> matrix::gaussian() {
  size_t _N_ = N_ - augmented_;

  std::map<size_t, real> pivots;
  size_t p_idx = 0;

  // Row Echelon Form
  for (size_t i = 0; i < M_; ++i) {
    size_t max_row = i;

    // Partial pivoting
    for (size_t j = i + 1; j < M_; ++j)
      if (std::abs(container_[j][p_idx]) > std::abs(container_[max_row][p_idx]))
        max_row = j;

    if (max_row != i)
      elementary<1>(i + 1, max_row + 1);

    std::vector<real> &row = container_[i];

    // Attempt to create the pivot on the current row
    real pivot = row[p_idx];
    for (; pivot == 0 && p_idx < _N_; pivot = row[++p_idx])
      ;

    if (p_idx >= _N_)
      continue;

    // assert(std::abs(pivot) >= std::numeric_limits<real>::epsilon());

    // Eliminate the nonzero numbers below the current pivot
    for (size_t j = i + 1; j < M_; ++j) {
      std::vector<real> &curr_row = container_[j]; // Below current pivot
      real b_pivot = curr_row[p_idx];

      if (b_pivot == 0)
        continue;

      elementary<3>(j + 1, i + 1, (1.0 / pivot) * -b_pivot);
    }

    pivots.try_emplace(p_idx, pivot);
    ++p_idx;
  }

  if (augmented_) {
    // Invalid Case: Zero Row = N_on-Zero N_umber
    for (std::vector<real> &a_row : container_) {
      vector row(vector::container_type(a_row.begin(), a_row.end() - 1));
      real b = a_row[N_ - 1];

      // Program termination due to system inconsisntency
      if (row.is_zero() && b != 0) {
        std::cout << "error: the system is inconsisntent\n";
        exit(0);
      }
    }
  }

  return pivots;
}

std::map<size_t, real> matrix::gauss_jordan() {
  std::map<size_t, real> piv_map = gaussian();
  std::vector<std::pair<size_t, real>> pivots(piv_map.begin(), piv_map.end());

  // Reduced Row Echelon Form
  for (size_t i = pivots.size(); i-- > 0;) {
    auto &[p_idx, _] = pivots[i];

    real pivot = container_[i][p_idx];
    elementary<2>(i + 1, 1.0 / pivot);

    piv_map[p_idx] = container_[i][p_idx];

    // N_onzero elimination above the pivot
    for (size_t j = i; j-- > 0;) {
      auto &row_above = container_[j];
      real above = row_above[p_idx];

      if (above == 0)
        continue;

      elementary<3>(j + 1, i + 1, -above);
    }
  }

  return piv_map;
}

std::vector<vector> matrix::solve_linear(vector &&b) {
  assert(b.n() == M_);

  matrix augmented = augment(std::move(b));

  std::map<size_t, real> pivots = augmented.gauss_jordan();
  size_t pn = pivots.size();

  std::cout << "augmented: " << augmented << "\n";

  std::vector<vector> solution;

  bool homogeneous = b.is_zero();
  size_t sol_n = (N_ - pn) + !homogeneous;

  solution.reserve(sol_n);

  for (size_t k = 0; k < (N_ + !homogeneous); ++k) {
    auto it = pivots.find(k);
    if (it != pivots.end())
      continue;

    vector sln(N_);
    vector col = augmented.col_vec(k + 1);

    for (size_t i = 0; i < pn; ++i)
      sln[i] = -col[i];

    // Populate free variables
    sln[pn + solution.size()] = 1;

    solution.push_back(std::move(sln));
  }

  return solution;
}

std::vector<vector> matrix::column_space() {
  auto copy = matrix(container_type(container_));
  std::map<size_t, real> pivots = copy.gauss_jordan();

  std::vector<vector> columns;
  columns.reserve(N_);

  for (auto it = pivots.begin(); it != pivots.end(); ++it)
    columns.push_back(std::move(col_vec(it->first + 1)));

  return columns;
}

std::vector<vector> matrix::null_space() {
  return solve_linear(vector::zero(M_));
}

real matrix::cofactor_expansion() {
  assert(is_square());

  if (M_ == 1)
    return entry(1, 1);

  // Edge case - 2 x 2
  else if (M_ == 2)
    return (entry(1, 1) * entry(2, 2)) - ((entry(1, 2) * entry(2, 1)));

  else {
    vector top = row_vec(1);
    real determinant = 0;

    for (size_t r = 0; r < M_; ++r) {
      matrix submatrix(M_ - 1, N_ - 1);
      real cofactor = ((r % 2 == 0) ? 1 : -1) * top[r];

      for (size_t i = 1; i < M_; ++i) {
        size_t col = 0;
        for (size_t j = 0; j < N_; ++j) {
          if (j == r)
            continue;

          submatrix.container_[i - 1][col] = container_[i][j];
          ++col;
        }
      }

      determinant += cofactor * submatrix.cofactor_expansion();
    }

    return determinant;
  }
}

vector matrix::transform(const vector &x) const {
  assert(N_ == x.n());

  vector y(M_);

  for (size_t i = 0; i < N_; ++i) {
    const std::vector<real> &row = container_[i];
    vector row_vec(vector::container_type{row.begin(), row.end()});

    y[i] = row_vec.dot(x);
  }

  return y;
}

matrix matrix::operator+(const matrix &b) const { return add(b); }

matrix matrix::operator-() const { return *this * -1; }
matrix matrix::operator-(const matrix &b) const { return add(-b); }

matrix matrix::operator*(real b) const { return mult(b); }

matrix matrix::operator*(const matrix &b) const { return mult(b); }

std::ostream &operator<<(std::ostream &os, const matrix &mat) {
  os << "[";

  for (size_t i = 0; i < mat.M_; ++i) {
    os << "\n";
    for (size_t j = 0; j < mat.N_; ++j) {
      os << std::setw(8) << std::setprecision(4) << mat.entry(i + 1, j + 1);
      // << (j == N_ - 1 ? "" : " ");
    }
  }

  return os << "\n]";
}

} // namespace linalg
