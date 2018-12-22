/**
 * @file   aka_types.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Feb 17 2011
 * @date last modification: Wed Dec 09 2020
 *
 * @brief  description of the "simple" types
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_error.hh"
#include "aka_compatibilty_with_cpp_standard.hh"
#include "aka_fwd.hh"
#include "aka_math.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <array>
#include <initializer_list>
#include <numeric>
#include <type_traits>
/* -------------------------------------------------------------------------- */
#define EIGEN_MATRIXBASE_PLUGIN "aka_types_eigen_matrix_base_plugin.hh"
#define EIGEN_MATRIX_PLUGIN "aka_types_eigen_matrix_plugin.hh"
#define EIGEN_PLAINOBJECTBASE_PLUGIN "aka_types_eigen_plain_object_base_plugin.hh"
#include <Eigen/Dense>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_AKA_TYPES_HH_
#define AKANTU_AKA_TYPES_HH_

namespace akantu {

using Eigen::Ref;

template <typename T, Eigen::Index n = Eigen::Dynamic>
using Vector = Eigen::Matrix<T, n, 1>;

template <typename T, Eigen::Index m = Eigen::Dynamic , Eigen::Index n = Eigen::Dynamic>
using Matrix = Eigen::Matrix<T, m, n,
                             Eigen::AutoAlign | Eigen::ColMajor>;

template <typename T, Eigen::Index n = Eigen::Dynamic>
using VectorProxy = Eigen::Map<Vector<T, n>>;

template <typename T, Eigen::Index m = Eigen::Dynamic , Eigen::Index n = Eigen::Dynamic>
using MatrixProxy = Eigen::Map<Matrix<T, m, n>>;

using VectorXr = Vector<Real>;
using MatrixXr = Matrix<Real>;

enum NormType : int8_t { L_1 = 1, L_2 = 2, L_inf = -1 };

} // namespace akantu

namespace aka {
template <typename Derived>
using is_matrice = aka::bool_constant<not Derived::IsVectorAtCompileTime>;

template <typename Derived>
using is_vector = aka::bool_constant<Derived::IsVectorAtCompileTime>;

template <typename ... Ds>
using are_matrices = aka::conjunction<is_matrice<Ds>...>;

template <typename ... Ds>
using enable_if_matrices_t = std::enable_if_t<are_matrices<Ds...>::value>;
} // namespace aka

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename T, std::size_t ndim, bool is_proxy>
class TensorBase : public TensorTrait<ndim> {
  using RetType = TensorBase<T, ndim, is_proxy>;
  static_assert(ndim > 2, "TensorBase cannot by used for dimensions < 3");

protected:
  using idx_t = std::size_t;

public:
  using proxy = TensorBase<T, ndim, true>;

  template <size_t _ndim = ndim,
            std::enable_if_t<_ndim == 1 or _ndim == 2, int> = 0>
  TensorBase() {
    n.fill(0);
  }

  TensorBase() { n.fill(0); }

  // proxy constructor
  template <typename... Args,
            std::enable_if_t<
                aka::conjunction<aka::disjunction<
                    std::is_integral<Args>, std::is_enum<Args>>...>::value and
                    ndim == sizeof...(Args),
                int> = 0>
  constexpr TensorBase(T * data, Args... args)
      : n{idx_t(args)...}, _size(detail::product_all(args...)), values(data),
        wrapped(true) {}

  // tensor constructor
  template <typename... Args,
            std::enable_if_t<
                aka::conjunction<aka::disjunction<
                    std::is_integral<Args>, std::is_enum<Args>>...>::value and
                    ndim == sizeof...(Args) and not is_proxy,
                int> = 0>
  constexpr TensorBase(Args... args)
      : n{idx_t(args)...}, _size(detail::product_all(args...)),
        values(new T[_size]) {
    static_assert(
        std::is_trivially_constructible<T>{},
        "Cannot create a tensor on non trivially constructible types");
  }

  /* ------------------------------------------------------------------------ */
  virtual ~TensorBase() {
    if (not this->wrapped and not is_proxy)
      delete[] this->values;
  }

  // copy constructors ---------------------------------------------------------
  // non proxy -> non proxy
  constexpr TensorBase(const TensorBase & other)
      : n(other.n), _size(other._size), wrapped(false) {
    if (is_proxy) {
      values = other.values;
    } else {
      std::cout << "Warning: copy" << std::endl;
      values = new T[_size];
      std::copy(other.values, other.values + _size, values);
    }
  }

  // non proxy -> proxy
  template <bool is_proxy_other,
            std::enable_if_t<is_proxy_other != is_proxy and is_proxy, int> = 0>
  constexpr TensorBase(const TensorBase<T, ndim, is_proxy_other> & other)
      : n(other.n), _size(other._size), values(other.values), wrapped(true) {}

  // proxy -> non proxy (wrapped)
  template <
      bool is_proxy_other,
      std::enable_if_t<is_proxy_other != is_proxy and not is_proxy, int> = 0>
  explicit constexpr TensorBase(
      const TensorBase<T, ndim, is_proxy_other> & other)
      : n(other.n), _size(other._size), values(other.values), wrapped(true) {}

  // move constructors ---------------------------------------------------------
  // proxy -> proxy, non proxy -> non proxy
  TensorBase(TensorBase && other)
      : n(std::move(other.n)), _size(std::exchange(other._size, 0)),
        values(std::exchange(other.values, nullptr)),
        wrapped(std::move(other.wrapped)) {}

  // proxy -> non proxy (wrapped)
  template <bool is_proxy_other,
            std::enable_if_t<is_proxy_other and not is_proxy, int> = 0>
  TensorBase(TensorBase<T, ndim, is_proxy_other> && other)
      : n(std::move(other.n)), _size(std::exchange(other._size, 0)),
        values(std::exchange(other.values, nullptr)), wrapped(true) {}

  // non proxy -> proxy
  template <bool is_proxy_other,
            std::enable_if_t<not is_proxy_other and is_proxy, int> = 0>
  TensorBase(TensorBase<T, ndim, is_proxy_other> && /*other*/) {
    static_assert(not(not is_proxy_other and is_proxy),
                  "This is not a valid constructor, cannot create a proxy on "
                  "temporary memory");
  }

  // copy operator -------------------------------------------------------------
  /// operator= copy-and-swap
  TensorBase & operator=(const TensorBase & other) {
    if (&other == this)
      return *this;

    if (is_proxy or this->wrapped) {
      AKANTU_DEBUG_ASSERT(
          other.size() == this->size(),
          "You are trying to copy two tensors with different sizes "
              << _size << " != " << other._size);
    } else {
      std::cout << "Warning: operator= delete data" << std::endl;
      delete[] values;
      n = other.n;
      _size = other._size;
      static_assert(
          std::is_trivially_constructible<T>{},
          "Cannot create a tensor on non trivially constructible types");
      values = new T[_size];
    }

    static_assert(std::is_trivially_copyable<T>{},
                  "Cannot copy a tensor on non trivial types");

    std::copy(other.values, other.values + _size, values);

    return *this;
  }

protected:
  template <typename Array, std::size_t... I>
  constexpr auto check_indices(const Array & idx,
                               std::index_sequence<I...>) const {
    bool result = true;
    (void)std::initializer_list<int>{(result &= idx[I] < n[I], 0)...};
    return result;
  }

  template <typename... Args> constexpr auto compute_index(Args... args) const {
    std::array<idx_t, sizeof...(Args)> idx{idx_t(args)...};
    AKANTU_DEBUG_ASSERT(
        check_indices(idx, std::make_index_sequence<sizeof...(Args)>{}),
        "The indexes are out of bound");

    idx_t index = 0, i = (sizeof...(Args)) - 1;
    for (; i > 0; i--) {
      index += idx[i];
      if (i > 0)
        index *= n[i - 1];
    }
    return index + idx[0];
  }

  template <typename S, std::size_t... I>
  constexpr auto get_slice(idx_t s, std::index_sequence<I...>) {
    return S(this->values + s * detail::product_all(n[I]...), n[I]...);
  }

  template <typename S, std::size_t... I>
  constexpr auto get_slice(idx_t s, std::index_sequence<I...>) const {
    return S(this->values + s * detail::product_all(n[I]...), n[I]...);
  }

public:
  template <typename... Args,
            std::enable_if_t<
                aka::conjunction<aka::disjunction<
                    std::is_integral<Args>, std::is_enum<Args>>...>::value and
                    ndim == sizeof...(Args),
                int> = 0>
  inline T & operator()(Args... args) {
    return *(this->values + compute_index(std::move(args)...));
  }

  template <typename... Args,
            std::enable_if_t<
                aka::conjunction<aka::disjunction<
                    std::is_integral<Args>, std::is_enum<Args>>...>::value and
                    ndim == sizeof...(Args),
                int> = 0>
  inline const T & operator()(Args... args) const {
    return *(this->values + compute_index(std::move(args)...));
  }

  template <idx_t _ndim = ndim, std::enable_if_t<(_ndim > 3), int> = 0>
  inline auto operator()(idx_t s) {
    return get_slice<TensorBase<T, ndim - 1, true>>(
        std::move(s), std::make_index_sequence<ndim - 1>());
  }

  template <idx_t _ndim = ndim, std::enable_if_t<(_ndim > 3), int> = 0>
  inline auto operator()(idx_t s) const {
    return get_slice<TensorBase<T, ndim - 1, true>>(
        std::move(s), std::make_index_sequence<ndim - 1>());
  }

  template <idx_t _ndim = ndim, std::enable_if_t<_ndim == 3, int> = 0>
  inline auto operator()(idx_t s) {
    return get_slice<
        Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>>(
        std::move(s), std::make_index_sequence<ndim - 1>());
  }

  template <idx_t _ndim = ndim, std::enable_if_t<_ndim == 3, int> = 0>
  inline auto operator()(idx_t s) const {
    return get_slice<
        Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>>(
        std::move(s), std::make_index_sequence<ndim - 1>());
  }

protected:
  template <class Operator> RetType & transform(Operator && op) {
    std::transform(this->values, this->values + this->_size, this->values,
                   std::forward<Operator>(op));
    return *(static_cast<RetType *>(this));
  }

  template <class Other, class Operator>
  RetType & transform(Other && other, Operator && op) {
    AKANTU_DEBUG_ASSERT(_size == other.size(),
                        "The two tensors do not have the same size "
                            << this->_size << " != " << other._size);

    std::transform(this->values, this->values + this->_size, other.values,
                   this->values, std::forward<Operator>(op));
    return *(static_cast<RetType *>(this));
  }

  template <class Operator> T accumulate(T init, Operator && op) {
    return std::accumulate(this->values, this->values + this->_size,
                           std::move(init), std::forward<Operator>(op));
  }

  template <class Other, class Init, class Accumulate, class Operator>
  T transform_reduce(Other && other, T init, Accumulate && acc,
                     Operator && op) {
    return std::inner_product(
        this->values, this->values + this->_size, other.data(), std::move(init),
        std::forward<Accumulate>(acc), std::forward<Operator>(op));
  }

  // arithmetic operators ------------------------------------------------------
  /* ------------------------------------------------------------------------ */
public:
  template <bool ip>
  inline RetType & operator+=(const TensorBase<T, ndim, ip> & other) {
    return transform(other, [](auto && a, auto && b) { return a + b; });
  }

  /* ------------------------------------------------------------------------ */
  template <bool ip>
  inline RetType & operator-=(const TensorBase<T, ndim, ip> & other) {
    return transform(other, [](auto && a, auto && b) { return a - b; });
  }

  /* ------------------------------------------------------------------------ */
  inline RetType & operator+=(const T & x) {
    return transform([&x](auto && a) { return a + x; });
  }

  /* ------------------------------------------------------------------------ */
  inline RetType & operator-=(const T & x) {
    return transform([&x](auto && a) { return a - x; });
  }

  /* ------------------------------------------------------------------------ */
  inline RetType & operator*=(const T & x) {
    return transform([&x](auto && a) { return a * x; });
  }

  /* ---------------------------------------------------------------------- */
  inline RetType & operator/=(const T & x) {
    return transform([&x](auto && a) { return a / x; });
  }

  /// Y = \alpha X + Y
  template <bool ip>
  inline RetType & aXplusY(const TensorBase<T, ndim, ip> & other,
                           const T alpha = 1.) {
    return transform(other,
                     [&alpha](auto && a, auto && b) { return alpha * a + b; });
  }

  /* ------------------------------------------------------------------------ */

  template <typename TO, std::size_t ndim_o, bool is_proxy_o>
  friend class TensorBase;

  T * data() { return values; }
  const T * data() const { return values; }

  // clang-format off
  [[deprecated("use data instead to be stl compatible")]]
  T * storage() {
    return values;
  }

  [[deprecated("use data instead to be stl compatible")]]
  const T * storage() const {
    return values;
  }
  // clang-format on

  auto size() const { return _size; }
  auto size(idx_t i) const {
    AKANTU_DEBUG_ASSERT(i < ndim, "This tensor has only " << ndim
                                                          << " dimensions, not "
                                                          << (i + 1));
    return n[i];
  };

  inline void set(const T & t) { std::fill_n(values, _size, t); };
  inline void clear() { set(T()); };

  bool isWrapped() const { return is_proxy or this->wrapped; }

public:
  /// "Entrywise" norm norm<L_p> @f[ \|\boldsymbol{T}\|_p = \left(
  /// \sum_i^{n[0]}\sum_j^{n[1]}\sum_k^{n[2]} |T_{ijk}|^p \right)^{\frac{1}{p}}
  /// @f]
  template <NormType norm_type>
  auto norm() -> std::enable_if_t<norm_type == L_inf, T> const {
    return accumulate(
        T(), [](auto && init, auto && a) { return init + std::abs(a); });
  }

  template <NormType norm_type>
  auto norm() -> std::enable_if_t<norm_type == L_1, T> const {
    return accumulate(T(), [](auto && init, auto && a) {
      return std::max(init, std::abs(a));
    });
  }

  template <NormType norm_type>
  auto norm() -> std::enable_if_t<norm_type == L_2, T> const {
    return std::sqrt(
        accumulate(T(), [](auto && init, auto && a) { return init + a * a; }));
  }

  template <NormType norm_type>
  auto norm() -> std::enable_if_t<(norm_type > L_2), T> const {
    return std::pow(accumulate(T(),
                               [](auto && init, auto && a) {
                                 return init + std::pow(a, norm_type);
                               }),
                    1. / norm_type);
  }

  T norm() const { return norm<L_2>(); }

protected:
  template <std::size_t N, typename... Args,
            std::enable_if_t<(sizeof...(Args) == ndim), int> = 0>
  void serialize(std::ostream & stream, Args... args) const {
    stream << this->operator()(std::move(args)...);
  }

  template <std::size_t N, typename... Args,
            std::enable_if_t<(sizeof...(Args) < ndim), int> = 0>
  void serialize(std::ostream & stream, Args... args) const {
    stream << "[";
    for (idx_t i = 0; i < n[N]; ++i) {
      if (i != 0)
        stream << ",";
      serialize<N + 1>(stream, std::move(args)..., i);
    }
    stream << "]";
  }

public:
  void printself(std::ostream & stream) const { serialize<0>(stream); };

protected:
  // size per dimension
  std::array<idx_t, ndim> n;

  // total storage size
  idx_t _size{0};

  // actual data location
  T * values{nullptr};

  // keep for backward compatibility
  bool wrapped{false};
};

/* -------------------------------------------------------------------------- */
template<typename T>
using Tensor3 = TensorBase<T, 3, false>;
template<typename T>
using Tensor3Proxy = TensorBase<T, 3, true>;

/* -------------------------------------------------------------------------- */

} // namespace akantu

#include <iterator>

// namespace std {
// template <typename R>
// struct iterator_traits<::akantu::types::details::vector_iterator<R>> {
// protected:
//   using iterator = ::akantu::types::details::vector_iterator<R>;

// public:
//   using iterator_category = typename iterator::iterator_category;
//   using value_type = typename iterator::value_type;
//   using difference_type = typename iterator::difference_type;
//   using pointer = typename iterator::pointer;
//   using reference = typename iterator::reference;
// };
// } // namespace std
/* -------------------------------------------------------------------------- */
#include "aka_view_iterators.hh"

namespace Eigen {
namespace {
template <typename T> struct SliceMap {};

template <typename Derived> struct SliceMap<MatrixBase<Derived>> {
  using type = std::conditional_t<
      Derived::IsVectorAtCompileTime,
      std::conditional_t<std::is_const<Derived>::value,
                         const typename Derived::Scalar,
                         typename Derived::Scalar>,
      Map<std::conditional_t<
          std::is_const<Derived>::value,
          const Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, 1>,
          Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, 1>>>>;
};

template <typename T> using SliceMap_t = typename SliceMap<T>::type;
}

template <typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto)
MatrixBase<Derived>::begin() {
  return ::akantu::view_iterator<SliceMap_t<MatrixBase<Derived>>>(this->derived().data());
}

template <typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto)
MatrixBase<Derived>::end() {
  return ::akantu::view_iterator<SliceMap_t<MatrixBase<Derived>>>(this->derived().data() +
                                                                  this->size());
}

template <typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto)
MatrixBase<Derived>::begin() const {
  return ::akantu::const_view_iterator<SliceMap_t<MatrixBase<Derived>>>(
      this->derived().data());
}

template <typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto)
MatrixBase<Derived>::end() const {
  return ::akantu::const_view_iterator<SliceMap_t<MatrixBase<Derived>>>(
      this->derived().data() + this->size());
}
} // namespace Eigen

#endif /* AKANTU_AKA_TYPES_HH_ */
