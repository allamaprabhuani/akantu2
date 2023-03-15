/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 * 
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 * 
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

using size_type = Index;

EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void zero();

template <bool _is_vector = IsVectorAtCompileTime,
          std::enable_if_t<not _is_vector> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto) operator()(Index c) {
  auto & d = this->derived();
  return d.col(c);
}

template <bool _is_vector = IsVectorAtCompileTime,
          std::enable_if_t<not _is_vector> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto) operator()(Index c) const {
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

EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void eye(const Scalar & t = 1.) {
  (*this).noalias() =
      t *
      Matrix<Scalar, Derived::RowsAtCompileTime,
             Derived::ColsAtCompileTime>::Identity(this->rows(), this->cols());
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
  eigen_assert(rows() == cols() and rows() == other.rows() and
               cols() == other.cols());

  return this->cwiseProduct(other).sum();
}

template <typename OtherDerived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void
eig(MatrixBase<OtherDerived> & other) const;

template <typename D1, typename D2>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void
eig(MatrixBase<D1> & values, MatrixBase<D2> & vectors,
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
