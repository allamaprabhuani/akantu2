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
#include "aka_common.hh"
#include "aka_math.hh"
#include "aka_view_iterators.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <array>
#include <initializer_list>
//#include <iomanip>
#include <Eigen/Dense>
#include <numeric>
#include <type_traits>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_AKA_TYPES_HH_
#define AKANTU_AKA_TYPES_HH_

namespace akantu {

using Eigen::Ref;

template <typename T> using VectorXt = Eigen::Matrix<T, Eigen::Dynamic, 1>;
template <typename T>
using MatrixXt = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic,
                               Eigen::AutoAlign | Eigen::ColMajor>;

using VectorXr = Eigen::Matrix<Real, Eigen::Dynamic, 1>;
using MatrixXr = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic,
                               Eigen::AutoAlign | Eigen::ColMajor>;

enum NormType : int8_t { L_1 = 1, L_2 = 2, L_inf = -1 };

template <typename T, bool is_proxy> class VectorBase;
template <typename T, bool is_proxy> class MatrixBase;

namespace detail {
  struct NotEigen {};
  template <typename T, std::size_t ndim, bool is_proxy>
  using EigenParent = std::conditional_t<
      ndim == 2,
      std::conditional_t<
          is_proxy,
          Eigen::Map<std::conditional_t<std::is_const<T>::value,
                                        const MatrixXt<std::remove_cv_t<T>>,
                                        MatrixXt<std::remove_cv_t<T>>>>,
          MatrixXt<T>>,
      std::conditional_t<
          ndim == 1,
          std::conditional_t<
              is_proxy,
              Eigen::Map<std::conditional_t<std::is_const<T>::value,
                                            const VectorXt<std::remove_cv_t<T>>,
                                            VectorXt<std::remove_cv_t<T>>>>,
              VectorXt<T>>,
          NotEigen>>;
} // namespace detail

/* -------------------------------------------------------------------------- */
template <typename T, std::size_t ndim, bool is_proxy>
class TensorBase : public TensorTrait<ndim>,
                   public detail::EigenParent<T, ndim, is_proxy> {
  using RetType = TensorBase<T, ndim, is_proxy>;
  using EigenParent = detail::EigenParent<T, ndim, is_proxy>;

protected:
  using idx_t = std::size_t;

public:
  using proxy = TensorBase<T, ndim, true>;

  template <size_t _ndim = ndim,
            std::enable_if_t<_ndim == 1 or _ndim == 2, int> = 0>
  TensorBase() : EigenParent() {
    n.fill(0);
  }

  TensorBase() { n.fill(0); }

  // proxy constructor
  template <typename... Args,
            std::enable_if_t<
                aka::conjunction<aka::disjunction<
                    std::is_integral<Args>, std::is_enum<Args>>...>::value and
                    ndim == sizeof...(Args) and (ndim > 2),
                int> = 0>
  constexpr TensorBase(T * data, Args... args)
      : n{idx_t(args)...}, _size(detail::product_all(args...)), values(data),
        wrapped(true) {}

  template <typename... Args,
            std::enable_if_t<
                aka::conjunction<aka::disjunction<
                    std::is_integral<Args>, std::is_enum<Args>>...>::value and
                        ndim == sizeof...(Args) and ndim == 1 or
                    ndim == 2 and is_proxy,
                int> = 0>
  constexpr TensorBase(T * data, Args... args)
      : EigenParent(data, args...), n{idx_t(args)...},
        _size(detail::product_all(args...)), values(data), wrapped(true) {}

  // tensor constructor
  template <typename... Args,
            std::enable_if_t<
                aka::conjunction<aka::disjunction<
                    std::is_integral<Args>, std::is_enum<Args>>...>::value and
                    ndim == sizeof...(Args) and (ndim > 2) and not is_proxy,
                int> = 0>
  constexpr TensorBase(Args... args)
      : n{idx_t(args)...}, _size(detail::product_all(args...)),
        values(new T[_size]) {
    static_assert(
        std::is_trivially_constructible<T>{},
        "Cannot create a tensor on non trivially constructible types");
  }

  template <typename... Args,
            std::enable_if_t<
                aka::conjunction<aka::disjunction<
                    std::is_integral<Args>, std::is_enum<Args>>...>::value and
                        ndim == sizeof...(Args) and ndim == 2 or
                    ndim == 1 and not is_proxy,
                int> = 0>
  constexpr TensorBase(Args... args)
      : EigenParent(args...), n{idx_t(args)...},
        _size(detail::product_all(args...)), values(data()), wrapped(true) {
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
  template <size_t ndim_ = ndim, std::enable_if_t<(ndim_ > 2), int> = 0>
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
  template <bool is_proxy_other, std::enable_if_t<is_proxy_other != is_proxy and
                                                      is_proxy and (ndim > 2),
                                                  int> = 0>
  constexpr TensorBase(const TensorBase<T, ndim, is_proxy_other> & other)
      : n(other.n), _size(other._size), values(other.values), wrapped(true) {}

  // proxy -> non proxy (wrapped)
  template <
      bool is_proxy_other,
      std::enable_if_t<
          is_proxy_other != is_proxy and not is_proxy and (ndim > 2), int> = 0>
  explicit constexpr TensorBase(
      const TensorBase<T, ndim, is_proxy_other> & other)
      : n(other.n), _size(other._size), values(other.values), wrapped(true) {}

  // move constructors ---------------------------------------------------------
  // proxy -> proxy, non proxy -> non proxy
  template <size_t ndim_ = ndim, std::enable_if_t<(ndim_ > 2), int> = 0>
  TensorBase(TensorBase && other)
      : n(std::move(other.n)), _size(std::exchange(other._size, 0)),
        values(std::exchange(other.values, nullptr)),
        wrapped(std::move(other.wrapped)) {}

  // proxy -> non proxy (wrapped)
  template <
      bool is_proxy_other,
      std::enable_if_t<is_proxy_other and not is_proxy and (ndim > 2), int> = 0>
  TensorBase(TensorBase<T, ndim, is_proxy_other> && other)
      : n(std::move(other.n)), _size(std::exchange(other._size, 0)),
        values(std::exchange(other.values, nullptr)), wrapped(true) {}

  // non proxy -> proxy
  template <
      bool is_proxy_other,
      std::enable_if_t<not is_proxy_other and is_proxy and (ndim > 2), int> = 0>
  TensorBase(TensorBase<T, ndim, is_proxy_other> && /*other*/) {
    static_assert(not(not is_proxy_other and is_proxy),
                  "This is not a valid constructor, cannot create a proxy on "
                  "temporary memory");
  }

  // Eigen functionalities
  // This constructor allows you to construct TensorBase from Eigen
  // expressions

  template <typename OtherDerived, size_t _ndim = ndim,
            std::enable_if_t<_ndim == 1 or _ndim == 2, int> = 0>
  TensorBase(const Eigen::MatrixBase<OtherDerived> & other)
      : EigenParent(other) {
    this->_size = this->n[0] = this->rows();
    if (ndim == 2) {
      this->n[1] = this->cols();
      this->size_ *= this->n[1];
    }
    this->values = this->data();
  }

  // This method allows you to assign Eigen expressions to TensorBase
  template <typename OtherDerived, size_t _ndim = ndim,
            std::enable_if_t<_ndim == 1 or _ndim == 2, int> = 0>
  TensorBase & operator=(const Eigen::MatrixBase<OtherDerived> & other) {
    EigenParent::operator=(other);
    this->_size = this->n[0] = this->rows();
    if (ndim == 2) {
      this->n[1] = this->cols();
      this->size_ *= this->n[1];
    }
    this->values = this->data();
    return *this;
  }

  // copy operator -------------------------------------------------------------
  /// operator= copy-and-swap
  template <size_t ndim_ = ndim, std::enable_if_t<(ndim_ > 2), int> = 0>
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
                    ndim == sizeof...(Args) and (ndim > 2),
                int> = 0>
  inline T & operator()(Args... args) {
    return *(this->values + compute_index(std::move(args)...));
  }

  template <typename... Args,
            std::enable_if_t<
                aka::conjunction<aka::disjunction<
                    std::is_integral<Args>, std::is_enum<Args>>...>::value and
                    ndim == sizeof...(Args) and (ndim > 2),
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

  template <idx_t _ndim = ndim,
            std::enable_if_t<_ndim == 3 or _ndim == 2, int> = 0>
  inline auto operator()(idx_t s) {
    return get_slice<TensorBase<T, ndim - 1, true>>(
        std::move(s), std::make_index_sequence<ndim - 1>());
  }

  template <idx_t _ndim = ndim,
            std::enable_if_t<_ndim == 3 or _ndim == 2, int> = 0>
  inline auto operator()(idx_t s) const {
    return get_slice<TensorBase<T, ndim - 1, true>>(
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
template <typename T, bool is_proxy>
class TensorBase<T, 1, is_proxy> : public TensorTrait<1>,
                                   public detail::EigenParent<T, 1, is_proxy> {
  using RetType = TensorBase<T, 1, is_proxy>;
  using EigenParent = detail::EigenParent<T, 1, is_proxy>;

protected:
  using idx_t = std::size_t;

public:
  using proxy = TensorBase<T, 1, true>;

public:
  TensorBase() : EigenParent() {}

  using EigenParent::EigenParent;

  // Eigen functionalities
  // This constructor allows you to construct TensorBase from Eigen
  // expressions
  template <typename OtherDerived>
  TensorBase(const Eigen::MatrixBase<OtherDerived> & other)
      : EigenParent(other) {}

  // This method allows you to assign Eigen expressions to TensorBase
  template <typename OtherDerived>
  TensorBase & operator=(const Eigen::MatrixBase<OtherDerived> & other) {
    EigenParent::operator=(other);
    return *this;
  }

public:
  inline T & operator()(idx_t i) {
    return EigenParent::operator()(std::move(i));
  }

  inline const T & operator()(idx_t i) const {
    return EigenParent::operator()(std::move(i));
  }

public:
  using iterator = view_iterator<T>;
  using const_iterator = const_view_iterator<T>;

  iterator begin() { return view_iterator<T>(this->data()); }
  iterator end() { return view_iterator<T>(this->data() + this->size()); }

  const_iterator begin() const { return const_view_iterator<T>(this->data()); }
  const_iterator end() const {
    return const_view_iterator<T>(this->data() + this->size());
  }

public:
  // clang-format off
  [[deprecated("use data instead to be stl compatible")]]
  T * storage() {
    return this->data();
  }

  [[deprecated("use data instead to be stl compatible")]]
  const T * storage() const {
    return this->data();
  }
  // clang-format on

  auto size() const { return this->cols(); }
  auto size(idx_t i) const {
    AKANTU_DEBUG_ASSERT(i < 1, "This tensor has only "
                                   << 1 << " dimensions, not " << (i + 1));
    return this->cols();
  };

  inline void set(const T & t) { this->fill(t); };
  inline void clear() { set(T()); };

public:
  template <class Tensor> T distance(const Tensor & other) const {
    return (*this - other).norm();
  }

public:
  /* ------------------------------------------------------------------------ */
  template <bool ip, typename R = T,
            std::enable_if_t<std::is_floating_point<R>::value, int> = 0>
  inline bool equal(const VectorBase<R, ip> & v,
                    R tolerance = Math::getTolerance()) const {
    T * a = this->data();
    T * b = v.data();
    idx_t i = 0;
    while (i < this->_size && (std::abs(*(a++) - *(b++)) < tolerance))
      ++i;
    }
    return i == this->_size;
  }

  /* ------------------------------------------------------------------------ */
  template <bool ip, typename R = T,
            std::enable_if_t<std::is_floating_point<R>::value, int> = 0>
  inline short compare(const TensorBase<R, 1, ip> & v,
                       Real tolerance = Math::getTolerance()) const {
    T * a = this->data();
    T * b = v.data();
    for (UInt i(0); i < this->_size; ++i, ++a, ++b) {
      if (std::abs(*a - *b) > tolerance) {
        return (((*a - *b) > tolerance) ? 1 : -1);
      }
    }
    return 0;
  }

  template <bool ip, typename R = T,
            std::enable_if_t<not std::is_floating_point<R>::value, int> = 0>
  inline bool equal(const TensorBase<R, 1, ip> & v) const {
    return std::equal(this->values, this->values + this->_size, v.data());
  }

  /* ------------------------------------------------------------------------ */
  template <bool ip, typename R = T,
            std::enable_if_t<not std::is_floating_point<R>::value, int> = 0>
  inline short compare(const TensorBase<R, 1, ip> & v) const {
    T * a = this->data();
    T * b = v.data();
    for (idx_t i(0); i < this->_size; ++i, ++a, ++b) {
      if (*a > *b)
        return 1;
      else if (*a < *b)
        return -1;
    }
    return 0;
  }

  /* ------------------------------------------------------------------------ */
  template <bool ip>
  inline bool operator==(const TensorBase<T, 1, ip> & v) const {
    return equal(v);
  }
  template <bool ip>
  inline bool operator!=(const TensorBase<T, 1, ip> & v) const {
    return !operator==(v);
  }
  template <bool ip>
  inline bool operator<(const TensorBase<T, 1, ip> & v) const {
    return compare(v) == -1;
  }
  template <bool ip>
  inline bool operator>(const TensorBase<T, 1, ip> & v) const {
    return compare(v) == 1;
  }

  template <bool ip>
  inline bool operator<=(const TensorBase<T, 1, ip> & v) const {
    return compare(v) <= 0;
  }

  template <bool ip>
  inline bool operator>=(const TensorBase<T, 1, ip> & v) const {
    return compare(v) >= 0;
  }

public:
  void printself(std::ostream & stream) const {
    stream << "[";
    for (idx_t i = 0; i < this->cols(); ++i) {
      if (i != 0)
        stream << ",";
      stream << this->operator()(i);
    }
    stream << "]";
  };
};

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
template <typename T, bool is_proxy>
class TensorBase<T, 2, is_proxy> : public TensorTrait<2>,
                                   public detail::EigenParent<T, 2, is_proxy> {
  using EigenParent = detail::EigenParent<T, 2, is_proxy>;

protected:
  using idx_t = std::size_t;

public:
  using proxy = TensorBase<T, 2, true>;

public:
  TensorBase() : EigenParent() {}

  using EigenParent::EigenParent;

  // Eigen functionalities
  // This constructor allows you to construct TensorBase from Eigen
  // expressions
  template <typename OtherDerived>
  TensorBase(const Eigen::MatrixBase<OtherDerived> & other)
      : EigenParent(other) {}

  // This method allows you to assign Eigen expressions to TensorBase
  template <typename OtherDerived>
  TensorBase & operator=(const Eigen::MatrixBase<OtherDerived> & other) {
    EigenParent::operator=(other);
    return *this;
  }

  TensorBase(std::initializer_list<std::initializer_list<T>> list) {
    if (is_proxy)
      AKANTU_EXCEPTION("Cannot create a proxy class from a initializer_list");
    static_assert(std::is_trivially_copyable<T>{},
                  "Cannot create a tensor on non trivial types");
    idx_t m = list.size();
    idx_t n = 0;
    for (auto row : list) {
      n = std::max(n, row.size());
    }

    this->resize(m, n);

    idx_t i = 0, j = 0;
    for (auto & row : list) {
      for (auto & val : row) {
        this->operator()(i, j++) = val;
      }
      ++i;
      j = 0;
    }
  }

public:
  inline T & operator()(idx_t i, idx_t j) {
    return EigenParent::operator()(i, j);
  }

  inline const T & operator()(idx_t i, idx_t j) const {
    return EigenParent::operator()(i, j);
  }

  inline auto operator()(idx_t s) {
    return TensorBase<T, 1, true>(this->data() + s * this->cols(),
                                  this->cols());
  }

  inline auto operator()(idx_t s) const {
    return TensorBase<const T, 1, true>(this->data() + s * this->cols(),
                                        this->cols());
  }

  // clang-format off
  [[deprecated("use data instead to be stl compatible")]]
  T * storage() {
    return this->data();
  }

  [[deprecated("use data instead to be stl compatible")]]
  const T * storage() const {
    return this->data();
  }
  // clang-format on

  auto size() const { return this->cols() * this->rows(); }
  auto size(idx_t i) const {
    AKANTU_DEBUG_ASSERT(i < 2, "This tensor has only "
                                   << 2 << " dimensions, not " << (i + 1));
    return i == 0 ? this->cols() : this->rows();
  };

  inline void set(const T & t) { this->fill(t); };
  inline void clear() { set(T()); };

public:
  void printself(std::ostream & stream) const {
    stream << "[";
    for (idx_t i = 0; i < this->rows(); ++i) {
      if (i != 0)
        stream << ",";
      stream << this->operator()(i);
    }
    stream << "]";
  };
};

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
namespace types {
  namespace details {
    template <typename reference_> class vector_iterator {
    public:
      using difference_type = std::ptrdiff_t;
      using value_type = std::decay_t<reference_>;
      using pointer = value_type *;
      using reference = reference_;
      using iterator_category = std::input_iterator_tag;

      vector_iterator(pointer ptr) : ptr(ptr) {}

      // input iterator ++it
      vector_iterator & operator++() {
        ++ptr;
        return *this;
      }
      // input iterator it++
      vector_iterator operator++(int) {
        auto cpy = *this;
        ++ptr;
        return cpy;
      }
      vector_iterator & operator+=(int n) {
        ptr += n;
        return *this;
      }

      vector_iterator operator+(int n) {
        vector_iterator cpy(*this);
        cpy += n;
        return cpy;
      }

      // input iterator it != other_it
      bool operator!=(const vector_iterator & other) const {
        return ptr != other.ptr;
      }
      bool operator==(const vector_iterator & other) const {
        return ptr == other.ptr;
      }

      difference_type operator-(const vector_iterator & other) const {
        return this->ptr - other.ptr;
      }

      // input iterator dereference *it
      reference operator*() { return *ptr; }
      pointer operator->() { return ptr; }

    private:
      pointer ptr;
    };
  } // namespace details
} // namespace types

// /* --------------------------------------------------------------------------
// */
// /* --------------------------------------------------------------------------
// */ template <typename T, bool is_proxy> class VectorBase : public
// TensorBase<T, 1, is_proxy> {
//   using Parent = TensorBase<T, 1, is_proxy>;
//   using EigenView = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>>;

// protected:
//   using idx_t = typename Parent::idx_t;

// public:
//   template <typename... Args>
//   VectorBase(Args &&... args) : Parent(std::forward<Args>(args)...) {}

//   operator EigenView() { return EigenView(this->values, this->_size); }

// public:
//   using iterator = types::details::vector_iterator<T &>;
//   using const_iterator = types::details::vector_iterator<const T &>;

//   iterator begin() { return iterator(this->data()); }
//   iterator end() { return iterator(this->data() + this->size()); }

//   const_iterator begin() const { return const_iterator(this->data()); }
//   const_iterator end() const {
//     return const_iterator(this->data() + this->size());
//   }

// public:
//   /* ------------------------------------------------------------------------
//   */ template <bool ip> inline T dot(const VectorBase<T, ip> & vect) const {
//     return Math::vectorDot(this->values, vect.data(), this->_size);
//   }

//   /* ------------------------------------------------------------------------
//   */ inline T mean() const {
//     return accumulate(T(), [](auto && init, auto && a) { return init + a; })
//     /
//            this->_size;
//   }

//   template <bool tr, bool ip1, bool ip2>
//   inline void mul(const MatrixBase<T, ip1> & A, const VectorBase<T, ip2> & x,
//                   T alpha = 1);

//   /* ------------------------------------------------------------------------
//   */ template <bool ip1, bool ip2> inline VectorBase & crossProduct(const
//   VectorBase<T, ip1> & v1,
//                                    const VectorBase<T, ip2> & v2) {
//     AKANTU_DEBUG_ASSERT(this->size() == 3,
//                         "crossProduct is only defined in 3D (n=" <<
//                         this->size()
//                                                                  << ")");
//     AKANTU_DEBUG_ASSERT(
//         this->size() == v1.size() && this->size() == v2.size(),
//         "crossProduct is not a valid operation non matching size vectors");
//     Math::vectorProduct3(v1.data(), v2.data(), this->values);
//     return *this;
//   }

//   template <bool ip> inline auto crossProduct(const VectorBase<T, ip> & v) {
//     VectorBase<T, false> tmp(this->size());
//     tmp.crossProduct(*this, v);
//     return tmp;
//   }

//   /* ------------------------------------------------------------------------
//   */ template <bool ip1, bool ip2> inline void solve(const MatrixBase<T, ip1>
//   & A,
//                     const VectorBase<T, ip2> & b) {
//     AKANTU_DEBUG_ASSERT(
//         this->size() == A.rows() && this->_size == A.cols(),
//         "The size of the solution vector mismatches the size of the matrix");
//     AKANTU_DEBUG_ASSERT(
//         this->_size == b._size,
//         "The rhs vector has a mismatch in size with the matrix");
//     Math::solve(this->_size, A.data(), this->values, b.data());
//   }

//   /* ------------------------------------------------------------------------
//   */ inline VectorBase & normalize() {
//     T n = this->norm();
//     this->operator/=(n);
//     return *this;
//   }

//   /* ------------------------------------------------------------------------
//   */
//   /// norm of (*this - x)
//   template <bool ip> inline T distance(const VectorBase<T, ip> & y) const {
//     return std::sqrt(transform_reduce(
//         y, T(), [](auto && init, auto && a) { return init + a; },
//         [](auto && a, auto && b) { return (a - b) * (a - b); }));
//   }

//   /* ------------------------------------------------------------------------
//   */ template <bool ip, typename R = T,
//             std::enable_if_t<std::is_floating_point<R>::value, int> = 0>
//   inline bool equal(const VectorBase<R, ip> & v,
//                     R tolerance = Math::getTolerance()) const {
//     T * a = this->data();
//     T * b = v.data();
//     idx_t i = 0;
//     while (i < this->_size && (std::abs(*(a++) - *(b++)) < tolerance))
//       ++i;
//     return i == this->_size;
//   }

//   /* ------------------------------------------------------------------------
//   */ template <bool ip, typename R = T,
//             std::enable_if_t<std::is_floating_point<R>::value, int> = 0>
//   inline short compare(const VectorBase<R, ip> & v,
//                        Real tolerance = Math::getTolerance()) const {
//     T * a = this->data();
//     T * b = v.data();
//     for (UInt i(0); i < this->_size; ++i, ++a, ++b) {
//       if (std::abs(*a - *b) > tolerance)
//         return (((*a - *b) > tolerance) ? 1 : -1);
//     }
//     return 0;
//   }

//   template <bool ip, typename R = T,
//             std::enable_if_t<not std::is_floating_point<R>::value, int> = 0>
//   inline bool equal(const VectorBase<R, ip> & v) const {
//     return std::equal(this->values, this->values + this->_size, v.data());
//   }

//   /* ------------------------------------------------------------------------
//   */ template <bool ip, typename R = T,
//             std::enable_if_t<not std::is_floating_point<R>::value, int> = 0>
//   inline short compare(const VectorBase<R, ip> & v) const {
//     T * a = this->data();
//     T * b = v.data();
//     for (idx_t i(0); i < this->_size; ++i, ++a, ++b) {
//       if (*a > *b)
//         return 1;
//       else if (*a < *b)
//         return -1;
//     }
//     return 0;
//   }

//   /* ------------------------------------------------------------------------
//   */ template <bool ip> inline VectorBase & operator*=(const VectorBase<T,
//   ip> & other) const {
//     return transform(other, [](auto && a, auto && b) { return a * b; });
//   }

//   /* ------------------------------------------------------------------------
//   */ template <bool ip> inline bool operator==(const VectorBase<T, ip> & v)
//   const {
//     return equal(v);
//   }
//   template <bool ip> inline bool operator!=(const VectorBase<T, ip> & v)
//   const {
//     return !operator==(v);
//   }
//   template <bool ip> inline bool operator<(const VectorBase<T, ip> & v) const
//   {
//     return compare(v) == -1;
//   }
//   template <bool ip> inline bool operator>(const VectorBase<T, ip> & v) const
//   {
//     return compare(v) == 1;
//   }

//   template <bool ip> inline bool operator<=(const VectorBase<T, ip> & v)
//   const {
//     return compare(v) <= 0;
//   }

//   template <bool ip> inline bool operator>=(const VectorBase<T, ip> & v)
//   const {
//     return compare(v) >= 0;
//   }
// };

// /* --------------------------------------------------------------------------
// */ template <typename T> class Vector : public VectorBase<T, false> {
//   using Parent = VectorBase<T, false>;

// protected:
//   using idx_t = typename Parent::idx_t;

// public:
//   Vector(std::initializer_list<T> list) : Parent(list.size()) {
//     idx_t i = 0;
//     for (auto val : list) {
//       this->operator()(i++) = val;
//     }
//   }
//   using Parent::Parent;

//   template <bool ip>
//   Vector(const TensorBase<T, 1, ip> & other) : Parent(other) {}

//   static inline Vector<T> zeros(idx_t n) {
//     Vector<T> tmp(n);
//     tmp.clear();
//     return tmp;
//   }
// };

// template <typename T> using VectorProxy = VectorBase<T, true>;

// /* --------------------------------------------------------------------------
// */
// /* --------------------------------------------------------------------------
// */ template <typename T, bool is_proxy> class MatrixBase : public
// TensorBase<T, 2, is_proxy> {
//   using Parent = TensorBase<T, 2, is_proxy>;
//   using EigenView = Eigen::Map<
//       Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>;

// protected:
//   using idx_t = typename Parent::idx_t;

// public:
//   using Parent::Parent;

//   /* ----------------------------------------------------------------------
//   */ idx_t rows() const { return this->n[0]; } idx_t cols() const { return
//   this->n[1]; }

//   operator EigenView() {
//     return EigenView(this->values, this->n[0], this->n[1]);
//   }

//   /* ----------------------------------------------------------------------
//   */ template <bool ip> inline MatrixBase<T, false> operator*(const
//   MatrixBase<T, ip> & B) {
//     MatrixBase<T, false> C(this->rows(), B.cols());
//     C.mul<false, false>(*this, B);
//     return C;
//   }

//   /* -----------------------------------------------------------------------
//   */ template <bool ip> inline MatrixBase & operator*=(const T & x) {
//     return Parent::operator*=(x);
//   }

//   template <bool ip> inline auto & operator*=(const MatrixBase<T, ip> & B) {
//     MatrixBase<T, false> C(*this);
//     this->mul<false, false>(C, B);
//     return *this;
//   }

//   /* ----------------------------------------------------------------------
//   */ template <bool tr_A, bool tr_B, bool ip1, bool ip2> inline void
//   mul(const MatrixBase<T, ip1> & A, const MatrixBase<T, ip2> & B,
//                   T alpha = 1.) {
//     auto k = A.cols();
//     if (tr_A)
//       k = A.rows();
// #ifndef AKANTU_NDEBUG
//     if (tr_B) {
//       AKANTU_DEBUG_ASSERT(k == B.cols(),
//                           "matrices to multiply have no fit dimensions");
//       AKANTU_DEBUG_ASSERT(this->cols() == B.rows(),
//                           "matrices to multiply have no fit dimensions");
//     } else {
//       AKANTU_DEBUG_ASSERT(k == B.rows(),
//                           "matrices to multiply have no fit dimensions");
//       AKANTU_DEBUG_ASSERT(this->cols() == B.cols(),
//                           "matrices to multiply have no fit dimensions");
//     }
//     if (tr_A) {
//       AKANTU_DEBUG_ASSERT(this->rows() == A.cols(),
//                           "matrices to multiply have no fit dimensions");
//     } else {
//       AKANTU_DEBUG_ASSERT(this->rows() == A.rows(),
//                           "matrices to multiply have no fit dimensions");
//     }
// #endif // AKANTU_NDEBUG
//     Math::matMul<tr_A, tr_B>(this->rows(), this->cols(), k, alpha, A.data(),
//                              B.data(), 0., this->data());
//   }

//   /* ----------------------------------------------------------------------
//   */ template <bool ip1, bool ip2> inline void outerProduct(const
//   VectorBase<T, ip1> & A,
//                            const VectorBase<T, ip2> & B) {
//     AKANTU_DEBUG_ASSERT(
//         A.size() == this->rows() && B.size() == this->cols(),
//         "A and B are not compatible with the size of the matrix");
//     for (idx_t i = 0; i < this->rows(); ++i) {
//       for (idx_t j = 0; j < this->cols(); ++j) {
//         this->values[i + j * this->rows()] += A[i] * B[j];
//       }
//     }
//   }

// private:
//   class EigenSorter {
//   public:
//     EigenSorter(const Vector<T> & eigs) : eigs(eigs) {}

//     bool operator()(const UInt & a, const UInt & b) const {
//       return (eigs(a) > eigs(b));
//     }

//   private:
//     const Vector<T> & eigs;
//   };

// public:
//   /* ----------------------------------------------------------------------
//   */ template <bool ip1, bool ip2> inline void eig(VectorBase<T, ip1> &
//   eigenvalues,
//                   MatrixBase<T, ip2> & eigenvectors) const {
//     AKANTU_DEBUG_ASSERT(this->cols() == this->rows(),
//                         "eig is not a valid operation on a rectangular
//                         matrix");
//     AKANTU_DEBUG_ASSERT(eigenvalues.size() == this->cols(),
//                         "eigenvalues should be of size " << this->cols()
//                                                          << ".");
//     AKANTU_DEBUG_ASSERT(eigenvectors.data() != nullptr and
//                             (eigenvectors.rows() == eigenvectors.cols()) and
//                             (eigenvectors.rows() == this->cols()),
//                         "Eigenvectors needs to be a square matrix of size "
//                             << this->cols() << " x " << this->cols() << ".");

//     MatrixBase<T, false> tmp = *this;
//     Vector<T> tmp_eigs(eigenvalues.size());
//     MatrixBase<T, false> tmp_eig_vects(eigenvectors.rows(),
//                                        eigenvectors.cols());

//     if (tmp_eig_vects.rows() == 0 || tmp_eig_vects.cols() == 0)
//       Math::matrixEig(tmp.cols(), tmp.data(), tmp_eigs.data());
//     else
//       Math::matrixEig(tmp.cols(), tmp.data(), tmp_eigs.data(),
//                       tmp_eig_vects.data());

//     Vector<UInt> perm(eigenvalues.size());
//     for (idx_t i = 0; i < perm.size(); ++i)
//       perm(i) = i;

//     std::sort(perm.begin(), perm.end(), EigenSorter(tmp_eigs));

//     for (UInt i = 0; i < perm.size(); ++i)
//       eigenvalues(i) = tmp_eigs(perm(i));

//     if (tmp_eig_vects.rows() != 0 && tmp_eig_vects.cols() != 0)
//       for (UInt i = 0; i < perm.size(); ++i) {
//         for (UInt j = 0; j < eigenvectors.rows(); ++j) {
//           eigenvectors(j, i) = tmp_eig_vects(j, perm(i));
//         }
//       }
//   }

//   /* ----------------------------------------------------------------------
//   */ template <bool ip> inline void eig(VectorBase<T, ip> & eigenvalues)
//   const {
//     MatrixBase<T, false> empty;
//     eig(eigenvalues, empty);
//   }

//   /* ----------------------------------------------------------------------
//   */ inline void eye(T alpha = 1.) {
//     AKANTU_DEBUG_ASSERT(this->cols() == this->rows(),
//                         "eye is not a valid operation on a rectangular
//                         matrix");
//     this->clear();
//     for (idx_t i = 0; i < this->cols(); ++i) {
//       this->values[i + i * this->rows()] = alpha;
//     }
//   }

//   /* ----------------------------------------------------------------------
//   */ inline T trace() const {
//     AKANTU_DEBUG_ASSERT(
//         this->cols() == this->rows(),
//         "trace is not a valid operation on a rectangular matrix");
//     T trace = 0.;
//     for (idx_t i = 0; i < this->rows(); ++i) {
//       trace += this->values[i + i * this->rows()];
//     }
//     return trace;
//   }

//   /* ----------------------------------------------------------------------
//   */ inline auto transpose() const {
//     MatrixBase<T, false> tmp(this->cols(), this->rows());
//     for (UInt i = 0; i < this->rows(); ++i) {
//       for (UInt j = 0; j < this->cols(); ++j) {
//         tmp(j, i) = this->operator()(i, j);
//       }
//     }
//     return tmp;
//   }

//   /* ----------------------------------------------------------------------
//   */ template <bool ip> inline void inverse(const MatrixBase<T, ip> & A) {
//     AKANTU_DEBUG_ASSERT(A.cols() == A.rows(),
//                         "inv is not a valid operation on a rectangular
//                         matrix");
//     AKANTU_DEBUG_ASSERT(this->cols() == A.cols(),
//                         "the matrix should have the same size as its
//                         inverse");

//     if (this->cols() == 1)
//       *this->values = 1. / *A.data();
//     else if (this->cols() == 2)
//       Math::inv2(A.data(), this->values);
//     else if (this->cols() == 3)
//       Math::inv3(A.data(), this->values);
//     else
//       Math::inv(this->cols(), A.data(), this->values);
//   }

//   inline auto inverse() {
//     MatrixBase<T, false> inv(this->rows(), this->cols());
//     inv.inverse(*this);
//     return inv;
//   }

//   /* --------------------------------------------------------------------- */
//   inline T det() const {
//     AKANTU_DEBUG_ASSERT(this->cols() == this->rows(),
//                         "inv is not a valid operation on a rectangular
//                         matrix");
//     if (this->cols() == 1)
//       return *(this->values);
//     else if (this->cols() == 2)
//       return Math::det2(this->values);
//     else if (this->cols() == 3)
//       return Math::det3(this->values);
//     else
//       return Math::det(this->cols(), this->values);
//   }

//   /* --------------------------------------------------------------------- */
//   template <bool ip> inline T doubleDot(const MatrixBase<T, ip> & other)
//   const {
//     AKANTU_DEBUG_ASSERT(
//         this->cols() == this->rows(),
//         "doubleDot is not a valid operation on a rectangular matrix");
//     if (this->cols() == 1)
//       return *(this->values) * *(other.data());
//     else if (this->cols() == 2)
//       return Math::matrixDoubleDot22(this->values, other.data());
//     else if (this->cols() == 3)
//       return Math::matrixDoubleDot33(this->values, other.data());
//     else
//       AKANTU_ERROR("doubleDot is not defined for other spatial dimensions"
//                    << " than 1, 2 or 3.");
//     return T();
//   }

//   /* ------------------------------------------------------------------------
//   */ template <bool ip> inline void block(const MatrixBase<T, ip> & block,
//   UInt pos_i, UInt pos_j) {
//     AKANTU_DEBUG_ASSERT(pos_i + block.rows() <= rows(),
//                         "The block size or position are not correct");
//     AKANTU_DEBUG_ASSERT(pos_i + block.cols() <= cols(),
//                         "The block size or position are not correct");
//     for (idx_t i = 0; i < block.rows(); ++i)
//       for (idx_t j = 0; j < block.cols(); ++j)
//         (*this)(i + pos_i, j + pos_j) = block(i, j);
//   }

//   inline auto block(UInt pos_i, UInt pos_j, UInt block_rows,
//                     UInt block_cols) const {
//     AKANTU_DEBUG_ASSERT(pos_i + block_rows <= rows(),
//                         "The block size or position are not correct");
//     AKANTU_DEBUG_ASSERT(pos_i + block_cols <= cols(),
//                         "The block size or position are not correct");
//     MatrixBase<T, false> block(block_rows, block_cols);
//     for (UInt i = 0; i < block_rows; ++i)
//       for (UInt j = 0; j < block_cols; ++j)
//         block(i, j) = this->at(i + pos_i, j + pos_j);
//     return block;
//   }
// };

// /* --------------------------------------------------------------------------
// */ template <typename T> class Matrix : public MatrixBase<T, false> {
// protected:
//   using Parent = MatrixBase<T, false>;

// protected:
//   using idx_t = typename Parent::idx_t;

// public:
//   Matrix(std::initializer_list<std::initializer_list<T>> list) : Parent(0, 0)
//   {
//     static_assert(std::is_trivially_copyable<T>{},
//                   "Cannot create a tensor on non trivial types");
//     this->n[0] = list.size();
//     this->n[1] = 0;
//     for (auto row : list) {
//       this->n[1] = std::max(this->n[1], row.size());
//     }

//     this->_size = this->n[0] * this->n[1];
//     delete[] this->values; // delete of size 0 to avoid LSan to complain
//     this->values = new T[this->_size];
//     this->wrapped = false;

//     idx_t i = 0, j = 0;
//     for (auto & row : list) {
//       for (auto & val : row) {
//         this->operator()(i, j++) = val;
//       }
//       ++i;
//       j = 0;
//     }
//   }

//   using Parent::Parent;

//   template <bool ip>
//   Matrix(const TensorBase<T, 1, ip> & other) : Parent(other) {}

//   /* ----------------------------------------------------------------------
//   */ static inline auto eye(idx_t m, T alpha = 1.) {
//     Matrix<T> tmp(m, m);
//     tmp.eye(alpha);
//     return tmp;
//   }
// };

// template <typename T> using MatrixProxy = MatrixBase<T, true>;
// /* --------------------------------------------------------------------------
// */

// template <typename T, std::size_t N, bool is_proxy>
// std::ostream & operator<<(std::ostream & stream,
//                           const TensorBase<T, N, is_proxy> & tensor) {
//   tensor.printself(stream);
//   return stream;
// }

// /* ------------------------------------------------------------------------
// */ template <typename T, bool ip> template <bool tr_A, bool ip1, bool ip2>
// inline void VectorBase<T, ip>::mul(const MatrixBase<T, ip1> & A,
//                                    const VectorBase<T, ip2> & x, T alpha) {
// #ifndef AKANTU_NDEBUG
//   auto n = x.size();
//   if (tr_A) {
//     AKANTU_DEBUG_ASSERT(n == A.rows(),
//                         "matrix and vector to multiply have no fit
//                         dimensions");
//     AKANTU_DEBUG_ASSERT(this->size() == A.cols(),
//                         "matrix and vector to multiply have no fit
//                         dimensions");
//   } else {
//     AKANTU_DEBUG_ASSERT(n == A.cols(),
//                         "matrix and vector to multiply have no fit
//                         dimensions");
//     AKANTU_DEBUG_ASSERT(this->size() == A.rows(),
//                         "matrix and vector to multiply have no fit
//                         dimensions");
//   }
// #endif
//   Math::matVectMul<tr_A>(A.rows(), A.cols(), alpha, A.data(), x.data(),
//                          0., this->data());
// }

template <typename T> using Vector = TensorBase<T, 1, false>;
template <typename T> using VectorProxy = TensorBase<T, 1, true>;

template <typename T> using Matrix = TensorBase<T, 2, false>;
template <typename T> using MatrixProxy = TensorBase<T, 2, true>;

template <typename T, size_t n> using Tensor = TensorBase<T, n, false>;
template <typename T, size_t n> using TensorProxy = TensorBase<T, n, true>;

/* ------------------------------------------------------------------------ */
/* Tensor3                                                                  */
/* ------------------------------------------------------------------------ */
template <typename T> using Tensor3 = TensorBase<T, 3, false>;
template <typename T> using Tensor3Proxy = TensorBase<T, 3, true>;

/* -------------------------------------------------------------------------- */
// support operations for the creation of other vectors
/* -------------------------------------------------------------------------- */
// template <typename T, bool ip>
// auto operator*(const T & scalar, const VectorBase<T, ip> & a) {
//   Vector<T> r(a);
//   r *= scalar;
//   return r;
// }

// template <typename T, bool ip>
// auto operator*(const VectorBase<T, ip> & a, const T & scalar) {
//   Vector<T> r(a);
//   r *= scalar;
//   return r;
// }

// template <typename T, bool ip>
// auto operator/(const VectorBase<T, ip> & a, const T & scalar) {
//   Vector<T> r(a);
//   r /= scalar;
//   return r;
// }

// template <typename T, bool ip1, bool ip2>
// auto operator*(const VectorBase<T, ip1> & a, const VectorBase<T, ip2> & b) {
//   Vector<T> r(a);
//   r *= b;
//   return r;
// }

// template <typename T, bool ip1, bool ip2>
// auto operator+(const VectorBase<T, ip1> & a, const VectorBase<T, ip2> & b) {
//   Vector<T> r(a);
//   r += b;
//   return r;
// }

// template <typename T, bool ip1, bool ip2>
// auto operator-(const VectorBase<T, ip1> & a, const VectorBase<T, ip2> & b) {
//   Vector<T> r(a);
//   r -= b;
//   return r;
// }

// template <typename T, bool ip1, bool ip2>
// auto operator*(const MatrixBase<T, ip1> & A, const VectorBase<T, ip2> & b) {
//   Vector<T> r(b.size());
//   r.mul(A, b);
//   return r;
// }

// /* --------------------------------------------------------------------------
// */ template <typename T, bool ip> auto operator*(const T & scalar, const
// MatrixBase<T, ip> & a) {
//   Matrix<T> r(a);
//   r *= scalar;
//   return r;
// }

// template <typename T, bool ip>
// auto operator*(const MatrixBase<T, ip> & a, const T & scalar) {
//   Matrix<T> r(a);
//   r *= scalar;
//   return r;
// }

// template <typename T, bool ip>
// auto operator/(const MatrixBase<T, ip> & a, const T & scalar) {
//   Matrix<T> r(a);
//   r /= scalar;
//   return r;
// }

// template <typename T, bool ip1, bool ip2>
// auto operator+(const MatrixBase<T, ip1> & a, const MatrixBase<T, ip2> & b) {
//   Matrix<T> r(a);
//   r += b;
//   return r;
// }

// template <typename T, bool ip1, bool ip2>
// auto operator-(const MatrixBase<T, ip1> & a, const MatrixBase<T, ip2> & b) {
//   Matrix<T> r(a);
//   r -= b;
//   return r;
// }
/* -------------------------------------------------------------------------- */
namespace {
  template <std::size_t dim, typename T> struct ViewIteratorHelper {};

  template <typename T> struct ViewIteratorHelper<0, T> { using type = T; };
  template <typename T> struct ViewIteratorHelper<1, T> {
    using type = Vector<T>;
  };
  template <typename T> struct ViewIteratorHelper<2, T> {
    using type = Matrix<T>;
  };
  template <typename T> struct ViewIteratorHelper<3, T> {
    using type = Tensor3<T>;
  };

  template <std::size_t dim, typename T>
  using ViewIteratorHelper_t = typename ViewIteratorHelper<dim, T>::type;
} // namespace

} // namespace akantu

#include <iterator>

namespace std {
template <typename R>
struct iterator_traits<::akantu::types::details::vector_iterator<R>> {
protected:
  using iterator = ::akantu::types::details::vector_iterator<R>;

public:
  using iterator_category = typename iterator::iterator_category;
  using value_type = typename iterator::value_type;
  using difference_type = typename iterator::difference_type;
  using pointer = typename iterator::pointer;
  using reference = typename iterator::reference;
};

template <typename Mat>
struct iterator_traits<::akantu::types::details::column_iterator<Mat>> {
protected:
  using iterator = ::akantu::types::details::column_iterator<Mat>;

public:
  using iterator_category = typename iterator::iterator_category;
  using value_type = typename iterator::value_type;
  using difference_type = typename iterator::difference_type;
  using pointer = typename iterator::pointer;
  using reference = typename iterator::reference;
};

} // namespace std

#endif /* AKANTU_AKA_TYPES_HH_ */
