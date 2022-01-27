using size_type = Index;

EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void zero();

template <bool _is_vector = IsVectorAtCompileTime,
          std::enable_if_t<not _is_vector> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE auto operator()(Index c) {
  auto & d = this->derived();
  return d.col(c);
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

EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Index size() const {
  return this->rows() * this->cols();
}

template <bool _is_vector = IsVectorAtCompileTime,
          std::enable_if_t<not _is_vector> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Index size(Index i) const {
  AKANTU_DEBUG_ASSERT(i < 2, "This tensor has only " << 2 << " dimensions, not "
                                                     << (i + 1));
  return (i == 0) ? this->rows() : this->cols();
}

EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void set(const Scalar & t) {
  this->fill(t);
}

EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void eye(const Scalar & t = Scalar()) {
  (*this) = t * Matrix<Scalar, Derived::RowsAtCompileTime,
                       Derived::ColsAtCompileTime>::Identity(this->rows(),
                                                             this->cols());
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

  this->cwiseProduct(other).sum();
}

template <typename OtherDerived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void
eig(const MatrixBase<OtherDerived> & other) const;

template <typename D1, typename D2>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void
eig(const MatrixBase<D1> & values, const MatrixBase<D2> & vectors,
    bool sort = std::is_floating_point<typename D1::Scalar>::value) const;

template <typename OtherDerived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void
eigh(const MatrixBase<OtherDerived> & other) const;

template <typename D1, typename D2>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void eigh(const MatrixBase<D1> & values,
                                                const MatrixBase<D2> & vectors,
                                                bool sort = true) const;

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
