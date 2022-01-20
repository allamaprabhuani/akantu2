using size_type = Index;

template <bool _is_vector = IsVectorAtCompileTime,
          std::enable_if_t<_is_vector> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
MatrixBase(std::initializer_list<Scalar> list)
    : MatrixBase() {
  Index i = 0;
  for (auto val : list) {
    this->operator()(i++) = val;
  }
}

template <bool _is_vector = IsVectorAtCompileTime,
          std::enable_if_t<not _is_vector> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
MatrixBase(std::initializer_list<std::initializer_list<Scalar>> list)
    : MatrixBase() {
  static_assert(std::is_trivially_copyable<Scalar>{},
                "Cannot create a tensor on non trivial types");
  Index m = list.size();
  Index n = 0;
  for (auto row : list) {
    n = std::max(n, row.size());
  }

  if (RowsAtCompileTime != -1 and RowsAtCompileTime != m) {
    AKANTU_EXCEPTION(
        "The size of the matrix does not correspond to the initializer_list");
  }

  if (ColsAtCompileTime != -1 and ColsAtCompileTime != n) {
    AKANTU_EXCEPTION(
        "The size of the matrix does not correspond to the initializer_list");
  }

  if (RowsAtCompileTime != -1 or ColsAtCompileTime != -1) {
    this->resize(m, n);
  }

  Index i = 0, j = 0;
  for (auto & row : list) {
    for (auto & val : row) {
      this->operator()(i, j++) = val;
    }
    ++i;
    j = 0;
  }
}

EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void zero();

template <bool _is_vector = IsVectorAtCompileTime,
          std::enable_if_t<not _is_vector> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE auto operator()(Index c) {
  auto & d = this->derived();
  return Map<Matrix<Scalar, RowsAtCompileTime, 1>>(
      d.data() + c * d.outerStride(), d.outerStride());
}

template <bool _is_vector = IsVectorAtCompileTime,
          std::enable_if_t<not _is_vector> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE auto operator()(Index c) const {
  const auto & d = this->derived();
  return d.col(c);
}

template <bool _is_vector = IsVectorAtCompileTime,
          std::enable_if_t<_is_vector> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto) operator()(Index c) {
  return Base::operator()(c);
}

template <bool _is_vector = IsVectorAtCompileTime,
          std::enable_if_t<_is_vector> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto) operator()(Index c) const {
  return Base::operator()(c);
}

EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto) operator()(Index i,
                                                                Index j) {
  return Base::operator()(i, j);
}

EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto) operator()(Index i,
                                                                Index j) const {
  return Base::operator()(i, j);
}

template <
    typename ED = Derived,
    typename T = std::remove_reference_t<decltype(*std::declval<ED>().data())>,
    std::enable_if_t<not std::is_const<T>::value and
                     ED::IsVectorAtCompileTime> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto) begin();
template <
    typename ED = Derived,
    typename T = std::remove_reference_t<decltype(*std::declval<ED>().data())>,
    std::enable_if_t<not std::is_const<T>::value and
                     ED::IsVectorAtCompileTime> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto) end();

template <typename ED = Derived,
          std::enable_if_t<ED::IsVectorAtCompileTime> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto) begin() const;
template <typename ED = Derived,
          std::enable_if_t<ED::IsVectorAtCompileTime> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto) end() const;

template <
    typename ED = Derived,
    typename T = std::remove_reference_t<decltype(*std::declval<ED>().data())>,
    std::enable_if_t<not std::is_const<T>::value and
                     not ED::IsVectorAtCompileTime> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto) begin();
template <
    typename ED = Derived,
    typename T = std::remove_reference_t<decltype(*std::declval<ED>().data())>,
    std::enable_if_t<not std::is_const<T>::value and
                     not ED::IsVectorAtCompileTime> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto) end();

template <typename ED = Derived,
          std::enable_if_t<not ED::IsVectorAtCompileTime> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto) begin() const;
template <typename ED = Derived,
          std::enable_if_t<not ED::IsVectorAtCompileTime> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto) end() const;

// clang-format off
[[deprecated("use data instead to be stl compatible")]]
Scalar * storage() {
  return this->data();
}

[[deprecated("use data instead to be stl compatible")]]
const Scalar * storage() const {
  return this->data();
}
// clang-format on

template <bool _is_vector = IsVectorAtCompileTime,
          std::enable_if_t<_is_vector> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE auto size() const {
  return this->cols();
}

EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE auto size(Index i) const {
  AKANTU_DEBUG_ASSERT(i < 1, "This tensor has only " << 1 << " dimensions, not "
                                                     << (i + 1));
  return this->cols();
}

EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void set(const Scalar & t) {
  this->fill(t);
}

EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void clear() { this->fill(Scalar()); };

template <typename OtherDerived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE auto
distance(const MatrixBase<OtherDerived> & other) const {
  return (*this - other).norm();
}

template <typename OtherDerived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Scalar
doubleDot(const MatrixBase<OtherDerived> & other) const {
  Scalar sum = 0;
  eigen_assert(rows() == cols() and rows() == other.rows() and
               cols() == other.cols());

  for (Index i = 0; i < rows(); ++i) {
    for (Index j = 0; j < cols(); ++j) {
      sum += coeff(i, j) * other.coeff(i, j);
    }
  }
}

template <typename OtherDerived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void
eig(const MatrixBase<OtherDerived> & other) const;

template <typename D1, typename D2>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void
eig(const MatrixBase<D1> & values, const MatrixBase<D2> & vectors) const;

template <typename OtherDerived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void
eigh(const MatrixBase<OtherDerived> & other) const;

template <typename D1, typename D2>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void
eigh(const MatrixBase<D1> & values, const MatrixBase<D2> & vectors) const;
/*
public:
*/
// template <bool ip, typename R = T,
//           std::enable_if_t<std::is_floating_point<R>::value, int> = 0>
// inline bool equal(const VectorBase<R, ip> & v,
//                   R tolerance = Math::getTolerance()) const {
//   T * a = this->data();
//   T * b = v.data();
//   idx_t i = 0;
//   while (i < this->_size && (std::abs(*(a++) - *(b++)) < tolerance))
//     ++i;
//   return i == this->_size;
// }

// template <bool ip, typename R = T,
//           std::enable_if_t<std::is_floating_point<R>::value, int> = 0>
// inline short compare(const TensorBase<R, 1, ip> & v,
//                      Real tolerance = Math::getTolerance()) const {
//   T * a = this->data();
//   T * b = v.data();
//   for (UInt i(0); i < this->_size; ++i, ++a, ++b) {
//     if (std::abs(*a - *b) > tolerance)
//       return (((*a - *b) > tolerance) ? 1 : -1);
//   }
//   return 0;
// }

// template <bool ip, typename R = T,
//           std::enable_if_t<not std::is_floating_point<R>::value, int> = 0>
// inline bool equal(const TensorBase<R, 1, ip> & v) const {
//   return std::equal(this->values, this->values + this->_size, v.data());
// }

// template <bool ip, typename R = T,
//           std::enable_if_t<not std::is_floating_point<R>::value, int> = 0>
// inline short compare(const TensorBase<R, 1, ip> & v) const {
//   T * a = this->data();
//   T * b = v.data();
//   for (idx_t i(0); i < this->_size; ++i, ++a, ++b) {
//     if (*a > *b)
//       return 1;
//     else if (*a < *b)
//       return -1;
//   }
//   return 0;
// }

template <typename OtherDerived>
inline bool operator<=(const MatrixBase<OtherDerived> & v) const {
  return this->isMuchSmallerThan(v);
}

template <typename OtherDerived>
inline bool operator>=(const MatrixBase<OtherDerived> & v) const {
  return v.isMuchSmallerThan(*this);
}

template <typename OtherDerived>
inline bool operator<(const MatrixBase<OtherDerived> & v) const {
  return (*this <= v) and (*this != v);
}

template <typename OtherDerived>
inline bool operator>(const MatrixBase<OtherDerived> & v) const {
  return (*this >= v) and (*this != v);
}
