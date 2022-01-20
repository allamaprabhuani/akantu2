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
#include "aka_compatibilty_with_cpp_standard.hh"
#include "aka_error.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <array>
#include <initializer_list>
#include <numeric>
#include <type_traits>

#ifndef AKANTU_AKA_TYPES_HH
#define AKANTU_AKA_TYPES_HH

/* -------------------------------------------------------------------------- */
namespace aka {
template <typename T> struct is_eigen_map : public std::false_type {};
} // namespace aka

#define EIGEN_DEFAULT_DENSE_INDEX_TYPE akantu::Idx
#define EIGEN_DEFAULT_MATRIX_STORAGE_ORDER_OPTION Eigen::ColMajor
/* -------------------------------------------------------------------------- */
#define EIGEN_MATRIXBASE_PLUGIN "aka_types_eigen_matrix_base_plugin.hh"
#define EIGEN_MATRIX_PLUGIN "aka_types_eigen_matrix_plugin.hh"
#define EIGEN_PLAINOBJECTBASE_PLUGIN                                           \
  "aka_types_eigen_plain_object_base_plugin.hh"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
/* -------------------------------------------------------------------------- */

namespace akantu {

using Eigen::Ref;

template <typename T, Eigen::Index n = Eigen::Dynamic>
using Vector = Eigen::Matrix<T, n, 1>;

template <typename T, Eigen::Index m = Eigen::Dynamic,
          Eigen::Index n = Eigen::Dynamic>
using Matrix = Eigen::Matrix<T, m, n>;

template <typename T, Eigen::Index n = Eigen::Dynamic>
using VectorProxy =
    Eigen::Map<std::conditional_t<std::is_const<T>::value,
                                  const Vector<std::remove_const_t<T>, n>,
                                  Vector<std::remove_const_t<T>, n>>>;

template <typename T, Eigen::Index m = Eigen::Dynamic,
          Eigen::Index n = Eigen::Dynamic>
using MatrixProxy =
    Eigen::Map<std::conditional_t<std::is_const<T>::value,
                                  const Matrix<std::remove_const_t<T>, m, n>,
                                  Matrix<std::remove_const_t<T>, m, n>>>;

using VectorXr = Vector<Real>;
using MatrixXr = Matrix<Real>;

enum NormType : int8_t { L_1 = 1, L_2 = 2, L_inf = -1 };

struct TensorTraitBase {};

template <size_t n> struct TensorTrait : public TensorTraitBase {};

} // namespace akantu

namespace aka {
template <typename Derived>
using is_matrix = aka::bool_constant<not Derived::IsVectorAtCompileTime>;

template <typename Derived>
using is_vector = aka::bool_constant<Derived::IsVectorAtCompileTime>;

template <typename... Ds>
using are_vectors = aka::conjunction<is_vector<Ds>...>;

template <typename... Ds>
using are_matrices = aka::conjunction<is_matrix<Ds>...>;

template <typename... Ds>
using enable_if_matrices_t = std::enable_if_t<are_matrices<Ds...>::value>;

template <typename... Ds>
using enable_if_vectors_t = std::enable_if_t<are_vectors<Ds...>::value>;

// template <typename T> struct is_eigen_map : public std::false_type {};

template <typename PlainObjectType, int MapOptions, typename StrideType>
struct is_eigen_map<Eigen::Map<PlainObjectType, MapOptions, StrideType>>
    : public std::true_type {};
} // namespace aka

/* -------------------------------------------------------------------------- */
namespace aka {

template <typename T>
struct is_tensor : public std::is_base_of<akantu::TensorTraitBase, T> {};

template <typename PlainObjectType, int MapOptions, typename StrideType>
struct is_tensor<Eigen::Map<PlainObjectType, MapOptions, StrideType>>
    : public std::true_type {};

template <typename T, Eigen::Index m, Eigen::Index n>
struct is_tensor<Eigen::Matrix<T, m, n>> : public std::true_type {};

template <typename T, size_t n>
using is_tensor_n = std::is_base_of<akantu::TensorTrait<n>, T>;

//   template <typename T> using is_vector = is_tensor_n<T, 1>;

//   template <typename T> using is_matrix = is_tensor_n<T, 2>;

template <std::size_t n, typename T = void, typename... Ts>
using enable_if_tensors_n = std::enable_if<
    aka::conjunction<
        is_tensor_n<Ts, n>...,
        std::is_same<
            std::common_type_t<std::decay_t<typename Ts::value_type>...>,
            std::decay_t<typename Ts::value_type>>...>::value,
    T>;

template <std::size_t n, typename T = void, typename... Ts>
using enable_if_tensors_n_t = typename enable_if_tensors_n<n, T, Ts...>::type;

//   template <typename T = void, typename... Ms>
//   using enable_if_vectors_t = typename enable_if_tensors_n<1, T,
//   Ms...>::type; template <typename T = void, typename... Ms>

//   using enable_if_matricies_t = typename enable_if_tensors_n<2, T,
//   Ms...>::type;
} // namespace aka

namespace akantu {
template <typename T, Int ndim> class TensorBase;
template <typename T, Int ndim> class TensorProxy;
template <typename T, Int ndim> class Tensor;
} // namespace akantu

#include "aka_view_iterators.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename T, Int ndim> class TensorBase : public TensorTrait<ndim> {
  using RetType = TensorBase<T, ndim>;
  static_assert(ndim > 2, "TensorBase cannot by used for dimensions < 3");

protected:
  using idx_t = Int;

  template <typename... Args>
  using valid_args_t = typename std::enable_if<
      aka::conjunction<aka::disjunction<std::is_integral<Args>,
                                        std::is_enum<Args>>...>::value and
          ndim == sizeof...(Args),
      int>::type;

public:
  using proxy = TensorBase<T, ndim>;
  using size_type = Idx;
  template <Int _ndim = ndim,
            std::enable_if_t<_ndim == 1 or _ndim == 2, int> = 0>
  TensorBase() {
    n.fill(0);
  }

  TensorBase() { n.fill(0); }

  template <typename... Args, valid_args_t<Args...> = 0>
  constexpr TensorBase(Args... args)
      : n{idx_t(args)...}, _size(detail::product_all(args...)) {}

  constexpr TensorBase(const TensorBase & other)
      : n(other.n), _size(other._size), values(other.values) {}

  constexpr TensorBase(TensorBase && other) noexcept
      : n(std::move(other.n)), _size(std::exchange(other._size, 0)),
        values(std::exchange(other.values, nullptr)) {}

protected:
  template <typename Array, Int... I>
  constexpr auto check_indices(const Array & idx,
                               std::index_sequence<I...>) const {
    bool result = true;
    (void)std::initializer_list<int>{(result &= idx[I] < n[I], 0)...};
    return result;
  }

  template <typename... Args> constexpr auto compute_index(Args... args) const {
    std::array<idx_t, sizeof...(Args)> idx{idx_t(args)...};
    static_assert(
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

  template <typename S, int... I>
  constexpr auto get_slice(idx_t s, std::index_sequence<I...>) {
    return S(this->values + s * detail::product_all(n[I]...), n[I]...);
  }

  template <typename S, std::size_t... I>
  constexpr auto get_slice(idx_t s, std::index_sequence<I...>) const {
    return S(this->values + s * detail::product_all(n[I]...), n[I]...);
  }

public:
  template <typename... Args, valid_args_t<Args...> = 0>
  inline auto operator()(Args... args) -> T & {
    return *(this->values + compute_index(std::move(args)...));
  }

  template <typename... Args, valid_args_t<Args...> = 0>
  inline auto operator()(Args... args) const -> const T & {
    return *(this->values + compute_index(std::move(args)...));
  }

  template <
      class R = T, idx_t _ndim = ndim,
      std::enable_if_t<(_ndim > 3) and not std::is_const<R>::value> * = nullptr>
  inline auto operator()(idx_t s) {
    return get_slice<TensorProxy<T, ndim - 1>>(
        s, std::make_index_sequence<ndim - 1>());
  }

  template <idx_t _ndim = ndim, std::enable_if_t<(_ndim > 3)> * = nullptr>
  inline auto operator()(idx_t s) const {
    return get_slice<TensorProxy<T, ndim - 1>>(
        s, std::make_index_sequence<ndim - 1>());
  }

  template <class R = T, idx_t _ndim = ndim,
            std::enable_if_t<(_ndim == 3) and not std::is_const<R>::value> * =
                nullptr>
  inline auto operator()(idx_t s) {
    return get_slice<
        Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>>(
        s, std::make_index_sequence<ndim - 1>());
  }

  template <idx_t _ndim = ndim, std::enable_if_t<_ndim == 3> * = nullptr>
  inline auto operator()(idx_t s) const {
    return get_slice<Eigen::Map<const Eigen::Matrix<
        std::remove_const_t<T>, Eigen::Dynamic, Eigen::Dynamic>>>(
        s, std::make_index_sequence<ndim - 1>());
  }

protected:
  template <class Operator> auto transform(Operator && op) -> RetType & {
    std::transform(this->values, this->values + this->_size, this->values,
                   std::forward<Operator>(op));
    return *(static_cast<RetType *>(this));
  }

  template <class Other, class Operator>
  auto transform(Other && other, Operator && op) -> RetType & {
    AKANTU_DEBUG_ASSERT(_size == other.size(),
                        "The two tensors do not have the same size "
                            << this->_size << " != " << other._size);

    std::transform(this->values, this->values + this->_size, other.values,
                   this->values, std::forward<Operator>(op));
    return *(static_cast<RetType *>(this));
  }

  template <class Operator> auto accumulate(T init, Operator && op) -> T {
    return std::accumulate(this->values, this->values + this->_size,
                           std::move(init), std::forward<Operator>(op));
  }

  template <class Other, class Init, class Accumulate, class Operator>
  auto transform_reduce(Other && other, T init, Accumulate && acc,
                        Operator && op) -> T {
    return std::inner_product(
        this->values, this->values + this->_size, other.data(), std::move(init),
        std::forward<Accumulate>(acc), std::forward<Operator>(op));
  }

  // arithmetic operators ------------------------------------------------------
  /* ------------------------------------------------------------------------ */
public:
  inline decltype(auto) operator+=(const TensorBase & other) {
    return transform(other, [](auto && a, auto && b) { return a + b; });
  }

  /* ------------------------------------------------------------------------ */
  inline auto operator-=(const TensorBase & other) -> TensorBase & {
    return transform(other, [](auto && a, auto && b) { return a - b; });
  }

  /* ------------------------------------------------------------------------ */
  inline auto operator+=(const T & x) -> TensorBase & {
    return transform([&x](auto && a) { return a + x; });
  }

  /* ------------------------------------------------------------------------ */
  inline auto operator-=(const T & x) -> TensorBase & {
    return transform([&x](auto && a) { return a - x; });
  }

  /* ------------------------------------------------------------------------ */
  inline auto operator*=(const T & x) -> TensorBase & {
    return transform([&x](auto && a) { return a * x; });
  }

  /* ---------------------------------------------------------------------- */
  inline auto operator/=(const T & x) -> TensorBase & {
    return transform([&x](auto && a) { return a / x; });
  }

  /// Y = \alpha X + Y
  inline auto aXplusY(const TensorBase & other, const T alpha = 1.)
      -> TensorBase & {
    return transform(other,
                     [&alpha](auto && a, auto && b) { return alpha * a + b; });
  }

  /* ------------------------------------------------------------------------ */
  auto data() -> T * { return values; }
  auto data() const -> const T * { return values; }

  // clang-format off
  [[deprecated("use data instead to be stl compatible")]]
  auto storage() -> T*{
    return values;
  }

  [[deprecated("use data instead to be stl compatible")]]
  auto storage() const -> const T * {
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

  auto norm() const -> T { return norm<L_2>(); }

protected:
  template <Int N, typename... Args,
            std::enable_if_t<(sizeof...(Args) == ndim), int> = 0>
  void serialize(std::ostream & stream, Args... args) const {
    stream << this->operator()(std::move(args)...);
  }

  template <Int N, typename... Args,
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
  template <std::size_t... I>
  constexpr decltype(auto) begin(std::index_sequence<I...>) {
    return view_iterator<ViewIteratorHelper_t<sizeof...(I), T>>(values,
                                                                n[I]...);
  }

  template <std::size_t... I>
  constexpr decltype(auto) end(std::index_sequence<I...>) {
    return view_iterator<ViewIteratorHelper_t<sizeof...(I), T>>(values + _size,
                                                                n[I]...);
  }

  template <std::size_t... I>
  constexpr decltype(auto) begin(std::index_sequence<I...>) const {
    return const_view_iterator<ViewIteratorHelper_t<sizeof...(I), T>>(values,
                                                                      n[I]...);
  }

  template <std::size_t... I>
  constexpr decltype(auto) end(std::index_sequence<I...>) const {
    return const_view_iterator<ViewIteratorHelper_t<sizeof...(I), T>>(
        values + _size, n[I]...);
  }

public:
  decltype(auto) begin() { return begin(std::make_index_sequence<ndim - 1>{}); }
  decltype(auto) end() { return end(std::make_index_sequence<ndim - 1>{}); }

  decltype(auto) begin() const {
    return begin(std::make_index_sequence<ndim - 1>{});
  }
  decltype(auto) end() const {
    return end(std::make_index_sequence<ndim - 1>{});
  }

protected:
  // size per dimension
  std::array<idx_t, ndim> n;

  // total storage size
  idx_t _size{0};

  // actual data location
  T * values{nullptr};
};

/* -------------------------------------------------------------------------- */
template <typename T, Int ndim> class TensorProxy : public TensorBase<T, ndim> {
private:
  using parent = TensorBase<T, ndim>;

public:
  // proxy constructor
  template <typename... Args>
  constexpr TensorProxy(T * data, Args... args) : parent(args...) {
    this->values = data;
  }

  constexpr TensorProxy(const TensorProxy<T, ndim> & other) : parent(other) {
    this->values = other.values;
  }

  constexpr TensorProxy(const Tensor<T, ndim> & other) : parent(other) {
    this->values = other.values;
  }

  // move constructors ---------------------------------------------------------
  // proxy -> proxy
  TensorProxy(TensorProxy && other) noexcept : parent(other) {}

  auto operator=(const TensorBase<T, ndim> & other) -> TensorProxy & {
    AKANTU_DEBUG_ASSERT(
        other.size() == this->size(),
        "You are trying to copy too a tensors proxy with the wrong size "
            << this->_size << " != " << other._size);

    static_assert(std::is_trivially_copyable<T>{},
                  "Cannot copy a tensor on non trivial types");

    std::copy(other.values, other.values + this->_size, this->values);
    return *this;
  }
};

/* -------------------------------------------------------------------------- */
template <typename T, Int ndim> class Tensor : public TensorBase<T, ndim> {
private:
  using parent = TensorBase<T, ndim>;

public:
  template <typename... Args> constexpr Tensor(Args... args) : parent(args...) {
    static_assert(
        std::is_trivially_constructible<T>{},
        "Cannot create a tensor on non trivially constructible types");
    this->values = new T[this->_size];
  }

  /* ------------------------------------------------------------------------ */
  virtual ~Tensor() { delete[] this->values; }

  // copy constructors ---------------------------------------------------------
  constexpr Tensor(const Tensor & other) : parent(other) {
    this->values = new T[this->_size];
    std::copy(other.values, other.values + this->_size, this->values);
  }

  constexpr explicit Tensor(const TensorProxy<T, ndim> & other)
      : parent(other) {
    //    static_assert(false, "Copying data are you sure");
    this->values = new T[this->_size];
    std::copy(other.values, other.values + this->_size, this->values);
  }

  // move constructors ---------------------------------------------------------
  // proxy -> proxy, non proxy -> non proxy
  Tensor(Tensor && other) noexcept : parent(other) {}

  // copy operator -------------------------------------------------------------
  /// operator= copy-and-swap
  auto operator=(const TensorBase<T, ndim> & other) -> Tensor & {
    if (&other == this)
      return *this;

    std::cout << "Warning: operator= delete data" << std::endl;
    delete[] this->values;
    this->n = other.n;
    this->_size = other._size;

    static_assert(
        std::is_trivially_constructible<T>{},
        "Cannot create a tensor on non trivially constructible types");

    this->values = new T[this->_size];

    static_assert(std::is_trivially_copyable<T>{},
                  "Cannot copy a tensor on non trivial types");

    std::copy(other.values, other.values + this->_size, this->values);

    return *this;
  }
};

/* -------------------------------------------------------------------------- */
template <typename T> using Tensor3 = Tensor<T, 3>;
template <typename T> using Tensor3Proxy = TensorProxy<T, 3>;
template <typename T> using Tensor3Base = TensorBase<T, 3>;

class ArrayBase;

/* -------------------------------------------------------------------------- */
namespace details {
  template <typename T> struct MapPlainObjectType { using type = T; };

  template <typename PlainObjectType, int MapOptions, typename StrideType>
  struct MapPlainObjectType<
      Eigen::Map<PlainObjectType, MapOptions, StrideType>> {
    using type = PlainObjectType;
  };

  template <typename T>
  using MapPlainObjectType_t = typename MapPlainObjectType<T>::type;

  template <typename Scalar, Idx...> struct EigenMatrixViewHelper {};

  template <typename Scalar, Idx RowsAtCompileTime>
  struct EigenMatrixViewHelper<Scalar, RowsAtCompileTime> {
    using type = Eigen::Matrix<Scalar, RowsAtCompileTime, 1>;
  };

  template <typename Scalar, Idx RowsAtCompileTime, Idx ColsAtCompileTime>
  struct EigenMatrixViewHelper<Scalar, RowsAtCompileTime, ColsAtCompileTime> {
    using type = Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime>;
  };

  template <typename Scalar, Idx... sizes>
  using EigenMatrixViewHelper_t =
      typename EigenMatrixViewHelper<Scalar, sizes...>::type;

  template <typename Array, Idx... sizes> class EigenView {
    static_assert(sizeof...(sizes) == 1 or sizeof...(sizes) == 2,
                  "Eigen only supports Vector and Matrices");

  public:
    using size_type = typename std::decay_t<Array>::size_type;
    using value_type = typename std::decay_t<Array>::value_type;

    EigenView(Array && array, decltype(sizes)... sizes_)
        : array(array), sizes_(sizes_...) {}

    EigenView(Array && array) : array(array), sizes_(sizes...) {}

    EigenView(const EigenView & other) = default;
    EigenView(EigenView && other) noexcept = default;

    auto operator=(const EigenView & other) -> EigenView & = default;
    auto operator=(EigenView && other) noexcept -> EigenView & = default;

    template <typename T = std::remove_reference_t<
                  decltype(*std::declval<Array>().data())>,
              std::enable_if_t<not std::is_const<T>::value> * = nullptr>
    decltype(auto) begin() {
      return aka::make_from_tuple<::akantu::view_iterator<
          Eigen::Map<EigenMatrixViewHelper_t<value_type, sizes...>>>>(
          std::tuple_cat(std::make_tuple(array.get().data()), sizes_));
    }

    template <typename T = std::remove_reference_t<
                  decltype(*std::declval<Array>().data())>,
              std::enable_if_t<not std::is_const<T>::value> * = nullptr>
    decltype(auto) end() {
      return aka::make_from_tuple<::akantu::view_iterator<
          Eigen::Map<EigenMatrixViewHelper_t<value_type, sizes...>>>>(
          std::tuple_cat(std::make_tuple(array.get().data() + array_size()),
                         sizes_));
    }

    decltype(auto) begin() const {
      return aka::make_from_tuple<::akantu::view_iterator<
          Eigen::Map<const EigenMatrixViewHelper_t<value_type, sizes...>>>>(
          std::tuple_cat(std::make_tuple(array.get().data()), sizes_));
    }
    decltype(auto) end() const {
      return aka::make_from_tuple<::akantu::view_iterator<
          Eigen::Map<const EigenMatrixViewHelper_t<value_type, sizes...>>>>(
          std::tuple_cat(std::make_tuple(array.get().data() + array_size()),
                         sizes_));
    }

  private:
    template <
        class A = Array,
        std::enable_if_t<std::is_base_of<ArrayBase, std::decay_t<A>>::value> * =
            nullptr>
    auto array_size() {
      return array.get().size() * array.get().getNbComponent();
    }

    template <class A = Array,
              std::enable_if_t<not std::is_base_of<
                  ArrayBase, std::decay_t<A>>::value> * = nullptr>
    auto array_size() {
      return array.get().size();
    }

  private:
    std::reference_wrapper<std::remove_reference_t<Array>> array;
    std::tuple<decltype(sizes)...> sizes_;
  };

} // namespace details

template <Idx RowsAtCompileTime, typename Array>
decltype(auto) make_view(Array && array, Idx rows = RowsAtCompileTime) {
  return details::EigenView<Array, RowsAtCompileTime>(
      std::forward<Array>(array), rows);
}

template <Idx RowsAtCompileTime, Idx ColsAtCompileTime, typename Array>
decltype(auto) make_view(Array && array, Idx rows = RowsAtCompileTime,
                         Idx cols = ColsAtCompileTime) {
  return details::EigenView<Array, RowsAtCompileTime, ColsAtCompileTime>(
      std::forward<Array>(array), rows, cols);
}

} // namespace akantu

namespace Eigen {

template <typename Derived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void MatrixBase<Derived>::zero() {
  return this->fill(0.);
}

/* -------------------------------------------------------------------------- */
/* Vector                                                                     */
/* -------------------------------------------------------------------------- */
template <typename Derived>
template <typename ED, typename T,
          std::enable_if_t<not std::is_const<T>::value and
                           ED::IsVectorAtCompileTime> *>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto)
MatrixBase<Derived>::begin() {
  return ::akantu::view_iterator<typename Derived::Scalar>(
      this->derived().data());
}

template <typename Derived>
template <typename ED, typename T,
          std::enable_if_t<not std::is_const<T>::value and
                           ED::IsVectorAtCompileTime> *>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto)
MatrixBase<Derived>::end() {
  return ::akantu::view_iterator<typename Derived::Scalar>(
      this->derived().data() + this->size());
}

template <typename Derived>
template <typename ED, std::enable_if_t<ED::IsVectorAtCompileTime> *>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto)
MatrixBase<Derived>::begin() const {
  using Scalar = typename Derived::Scalar;
  return ::akantu::const_view_iterator<Scalar>(this->derived().data());
}

template <typename Derived>
template <typename ED, std::enable_if_t<ED::IsVectorAtCompileTime> *>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto)
MatrixBase<Derived>::end() const {
  using Scalar = typename Derived::Scalar;
  return ::akantu::const_view_iterator<Scalar>(this->derived().data() +
                                               this->size());
}

/* -------------------------------------------------------------------------- */
/* Matrix                                                                     */
/* -------------------------------------------------------------------------- */
template <typename Derived>
template <typename ED, typename T,
          std::enable_if_t<not std::is_const<T>::value and
                           not ED::IsVectorAtCompileTime> *>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto)
MatrixBase<Derived>::begin() {
  return ::akantu::view_iterator<
      Map<Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, 1>>>(
      this->derived().data(), this->rows());
}

template <typename Derived>
template <typename ED, typename T,
          std::enable_if_t<not std::is_const<T>::value and
                           not ED::IsVectorAtCompileTime> *>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto)
MatrixBase<Derived>::end() {
  return ::akantu::view_iterator<
      Map<Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, 1>>>(
      this->derived().data() + this->size(), this->rows());
}

template <typename Derived>
template <typename ED, std::enable_if_t<not ED::IsVectorAtCompileTime> *>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto)
MatrixBase<Derived>::begin() const {
  using Scalar = typename Derived::Scalar;
  return ::akantu::const_view_iterator<
      Map<const Matrix<Scalar, Derived::RowsAtCompileTime, 1>>>(
      const_cast<Scalar *>(this->derived().data()), this->rows());
}

template <typename Derived>
template <typename ED, std::enable_if_t<not ED::IsVectorAtCompileTime> *>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE decltype(auto)
MatrixBase<Derived>::end() const {
  using Scalar = typename Derived::Scalar;
  return ::akantu::const_view_iterator<
      Map<const Matrix<Scalar, Derived::RowsAtCompileTime, 1>>>(
      const_cast<Scalar *>(this->derived().data()) + this->size(),
      this->rows());
}

template <typename Derived>
template <typename OtherDerived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void
MatrixBase<Derived>::eig(const MatrixBase<OtherDerived> & values_) const {
  EigenSolver<akantu::details::MapPlainObjectType_t<std::decay_t<Derived>>>
      solver(*this, false);
  Eigen::MatrixBase<OtherDerived> & values =
      const_cast<Eigen::MatrixBase<OtherDerived> &>(
          values_); // as advised by the Eigen developers

  values = solver.eigenvalues();
}

template <typename Derived>
template <typename D1, typename D2>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void
MatrixBase<Derived>::eig(const MatrixBase<D1> & values_,
                         const MatrixBase<D2> & vectors_) const {
  EigenSolver<akantu::details::MapPlainObjectType_t<std::decay_t<Derived>>>
      solver(*this, false);

  Eigen::MatrixBase<D1> & values = const_cast<Eigen::MatrixBase<D1> &>(
      values_); // as advised by the Eigen developers
  Eigen::MatrixBase<D2> & vectors = const_cast<Eigen::MatrixBase<D2> &>(
      vectors_); // as advised by the Eigen developers

  values = solver.eigenvalues();
  vectors = solver.eigenvectors();
}

template <typename Derived>
template <typename OtherDerived>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void
MatrixBase<Derived>::eigh(const MatrixBase<OtherDerived> & values_) const {
  SelfAdjointEigenSolver<
      akantu::details::MapPlainObjectType_t<std::decay_t<Derived>>>
      solver(*this, false);
  Eigen::MatrixBase<OtherDerived> & values =
      const_cast<Eigen::MatrixBase<OtherDerived> &>(
          values_); // as advised by the Eigen developers
  values = solver.eigenvalues();
}

template <typename Derived>
template <typename D1, typename D2>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void
MatrixBase<Derived>::eigh(const MatrixBase<D1> & values_,
                          const MatrixBase<D2> & vectors_) const {
  SelfAdjointEigenSolver<
      akantu::details::MapPlainObjectType_t<std::decay_t<Derived>>>
      solver(*this, false);

  Eigen::MatrixBase<D1> & values = const_cast<Eigen::MatrixBase<D1> &>(
      values_); // as advised by the Eigen developers
  Eigen::MatrixBase<D2> & vectors = const_cast<Eigen::MatrixBase<D2> &>(
      vectors_); // as advised by the Eigen developers

  values = solver.eigenvalues();
  vectors = solver.eigenvectors();
}

} // namespace Eigen

namespace std {
template <typename POT1, typename POT2, int MapOptions, typename StrideType>
struct is_convertible<Eigen::Map<POT1, MapOptions, StrideType>,
                      Eigen::Map<POT2, MapOptions, StrideType>>
    : aka::bool_constant<is_convertible<POT1, POT2>::value> {};
} // namespace std

#endif /* AKANTU_AKA_TYPES_HH */
