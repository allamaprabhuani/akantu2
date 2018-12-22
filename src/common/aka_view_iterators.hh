/**
 * @file   aka_view_iterators.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Thu Nov 15 2018
 *
 * @brief View iterators
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_types.hh"
/* -------------------------------------------------------------------------- */
#include <memory>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_AKA_VIEW_ITERATORS_HH__
#define __AKANTU_AKA_VIEW_ITERATORS_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
/* Iterators                                                                  */
/* -------------------------------------------------------------------------- */
namespace detail {
  template <class R> struct IteratorHelper { static constexpr Int dim = 0; };

  template <class Derived> struct IteratorHelper<Eigen::MatrixBase<Derived>> {
  private:
    using T = typename Derived::Scalar;
    static constexpr Int m = Derived::RowsAtCompileTime;
    static constexpr Int n = Derived::ColsAtCompileTime;

  public:
    static constexpr Int dim =
        Eigen::MatrixBase<Derived>::IsVectorAtCompileTime ? 1 : 2;
    using pointer = T *;
    using proxy = Eigen::Map<Eigen::Matrix<T, m, n>>;
    using const_proxy = Eigen::Map<const Eigen::Matrix<T, m, n>>;
  };

  template <class Derived> struct IteratorHelper<Eigen::Map<Derived>> {
  private:
    using T = typename Derived::Scalar;
    static constexpr Int m = Derived::RowsAtCompileTime;
    static constexpr Int n = Derived::ColsAtCompileTime;

  public:
    static constexpr Int dim = Derived::IsVectorAtCompileTime ? 1 : 2;
    using pointer = T *;
    using proxy = Eigen::Map<Derived>;
    using const_proxy = Eigen::Map<const Derived>;
  };

  template <typename T, std::size_t _dim, bool is_proxy>
  struct IteratorHelper<TensorBase<T, _dim, is_proxy>> {
    static constexpr Int dim = _dim;
    using pointer = T *;
    using proxy = TensorBase<T, _dim, true>;
    using const_proxy = TensorBase<const T, _dim, true>;
  };

  /* --------------------------------------------------------------------------
   */
  template <class R, class daughter, class IR = R,
            Int dim = IteratorHelper<std::decay_t<R>>::dim>
  class internal_view_iterator {
  protected:
    using helper = IteratorHelper<std::decay_t<R>>;
    using internal_value_type = IR;
    using internal_pointer = IR *;

    using scalar_pointer = typename helper::pointer;
    using proxy_t = typename helper::proxy;
    using const_proxy_t = typename helper::const_proxy;

  public:
    using pointer = proxy_t *;
    using value_type = proxy_t;
    using reference = proxy_t &;
    using const_reference = const_proxy_t &;
    using difference_type = Int;
    using iterator_category = std::random_access_iterator_tag;

  private:
    template <std::size_t... I>
    constexpr auto get_new_const_proxy(scalar_pointer data,
                                 std::index_sequence<I...>) const {
      return const_proxy_t(data, dims[I]...);
    }

    constexpr auto get_new_const_proxy(scalar_pointer data) const {
      return get_new_const_proxy(data, std::make_index_sequence<dim>());
    }
    template <std::size_t... I>
    constexpr auto get_new_proxy(scalar_pointer data,
                                 std::index_sequence<I...>) {
      return proxy_t(data, dims[I]...);
    }

    constexpr auto get_new_proxy(scalar_pointer data) {
      return get_new_proxy(data, std::make_index_sequence<dim>());
    }

    template <typename T, std::size_t... I>
    constexpr void reset_proxy(T & t, scalar_pointer data,
                               std::index_sequence<I...>) const {
      new (&t) T(data, dims[I]...);
    }

    template <typename T> constexpr auto reset_proxy(T & t) const {
      return reset_proxy(t, this->ret_ptr, std::make_index_sequence<dim>());
    }

  public:
    template <typename... Ns>
    internal_view_iterator(scalar_pointer data, Ns... ns)
        : dims({Int(ns)...}), _offset(detail::product_all(std::forward<Ns>(ns)...)),
          initial(data), ret_ptr(data), proxy(data, ns...),
          const_proxy(data, ns...) {}

    internal_view_iterator() : internal_view_iterator(nullptr) {}

    internal_view_iterator(const internal_view_iterator & it)
        : proxy(get_new_proxy(it.ret_ptr)),
          const_proxy(get_new_const_proxy(it.ret_ptr)) {
      if (this != &it) {
        this->dims = it.dims;
        this->_offset = it._offset;
        this->initial = it.initial;
        this->ret_ptr = it.ret_ptr;
      }
    }

    internal_view_iterator(internal_view_iterator && it) = default;

    virtual ~internal_view_iterator() = default;

    inline internal_view_iterator &
    operator=(const internal_view_iterator & it) {
      if (this != &it) {
        this->_offset = it._offset;
        this->initial = it.initial;
        this->ret_ptr = it.ret_ptr;
        reset_proxy(this->proxy);
        reset_proxy(this->const_proxy);
      }
      return *this;
    }

  public:
    UInt getCurrentIndex() {
      return (this->ret_ptr - this->initial) / this->_offset;
    }

    inline reference operator*() {
      reset_proxy(proxy);
      return proxy;
    }

    inline const_reference operator*() const {
      reset_proxy(const_proxy);
      return proxy;
    }

    inline daughter & operator++() {
      ret_ptr += _offset;
      return static_cast<daughter &>(*this);
    }

    inline daughter & operator--() {
      ret_ptr -= _offset;
      return static_cast<daughter &>(*this);
    }

    inline daughter & operator+=(UInt n) {
      ret_ptr += _offset * n;
      return static_cast<daughter &>(*this);
    }

    inline daughter & operator-=(UInt n) {
      ret_ptr -= _offset * n;
      return static_cast<daughter &>(*this);
    }

    inline auto operator[](UInt n) {
      return get_new_proxy(ret_ptr + n * _offset);
    }

    inline auto operator[](UInt n) const {
      return get_new_proxy(ret_ptr + n * _offset);
    }

    inline bool operator==(const internal_view_iterator & other) const {
      return this->ret_ptr == other.ret_ptr;
    }

    inline bool operator!=(const internal_view_iterator & other) const {
      return this->ret_ptr != other.ret_ptr;
    }

    inline bool operator<(const internal_view_iterator & other) const {
      return this->ret_ptr < other.ret_ptr;
    }

    inline bool operator<=(const internal_view_iterator & other) const {
      return this->ret_ptr <= other.ret_ptr;
    }

    inline bool operator>(const internal_view_iterator & other) const {
      return this->ret_ptr > other.ret_ptr;
    }

    inline bool operator>=(const internal_view_iterator & other) const {
      return this->ret_ptr >= other.ret_ptr;
    }

    inline daughter operator+(difference_type n) {
      daughter tmp(static_cast<daughter &>(*this));
      tmp += n;
      return tmp;
    }

    inline daughter operator-(difference_type n) {
      daughter tmp(static_cast<daughter &>(*this));
      tmp -= n;
      return tmp;
    }

    inline difference_type operator-(const internal_view_iterator & b) {
      return (this->ret_ptr - b.ret_ptr) / _offset;
    }

    inline scalar_pointer data() const { return ret_ptr; }
    inline difference_type offset() const { return _offset; }

  protected:
    std::array<Int, dim> dims;
    difference_type _offset{0};
    scalar_pointer initial{nullptr};
    scalar_pointer ret_ptr{nullptr};
    proxy_t proxy;
    const_proxy_t const_proxy;
  };

  /* --------------------------------------------------------------------------
   */
  /**
   * Specialization for scalar types
   */
  template <class R, class daughter, class IR>
  class internal_view_iterator<R, daughter, IR, 0> {
  public:
    using value_type = R;
    using pointer = R *;
    using reference = R &;
    using const_reference = const R &;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::random_access_iterator_tag;

  protected:
    using internal_value_type = IR;
    using internal_pointer = IR *;

  public:
    internal_view_iterator(pointer data = nullptr) : ret(data), initial(data){};
    internal_view_iterator(const internal_view_iterator & it) = default;
    internal_view_iterator(internal_view_iterator && it) = default;

    virtual ~internal_view_iterator() = default;

    inline internal_view_iterator &
    operator=(const internal_view_iterator & it) = default;

    UInt getCurrentIndex() { return (this->ret - this->initial); };

    inline reference operator*() { return *ret; }
    inline const_reference operator*() const { return *ret; }
    inline pointer operator->() { return ret; };
    inline daughter & operator++() {
      ++ret;
      return static_cast<daughter &>(*this);
    }
    inline daughter & operator--() {
      --ret;
      return static_cast<daughter &>(*this);
    }

    inline daughter & operator+=(const UInt n) {
      ret += n;
      return static_cast<daughter &>(*this);
    }
    inline daughter & operator-=(const UInt n) {
      ret -= n;
      return static_cast<daughter &>(*this);
    }

    inline reference operator[](const UInt n) { return ret[n]; }

    inline bool operator==(const internal_view_iterator & other) const {
      return ret == other.ret;
    }
    inline bool operator!=(const internal_view_iterator & other) const {
      return ret != other.ret;
    }
    inline bool operator<(const internal_view_iterator & other) const {
      return ret < other.ret;
    }
    inline bool operator<=(const internal_view_iterator & other) const {
      return ret <= other.ret;
    }
    inline bool operator>(const internal_view_iterator & other) const {
      return ret > other.ret;
    }
    inline bool operator>=(const internal_view_iterator & other) const {
      return ret >= other.ret;
    }

    inline daughter operator-(difference_type n) { return daughter(ret - n); }
    inline daughter operator+(difference_type n) { return daughter(ret + n); }

    inline difference_type operator-(const internal_view_iterator & b) {
      return ret - b.ret;
    }

    inline pointer data() const { return ret; }

  protected:
    pointer ret{nullptr};
    pointer initial{nullptr};
  };
} // namespace detail

/* -------------------------------------------------------------------------- */
template <typename R> class view_iterator;

template <typename R>
class const_view_iterator
    : public detail::internal_view_iterator<const R, const_view_iterator<R>,
                                             R> {
public:
  using parent =
      detail::internal_view_iterator<const R, const_view_iterator, R>;
  using value_type = typename parent::value_type;
  using pointer = typename parent::pointer;
  using reference = typename parent::reference;
  using difference_type = typename parent::difference_type;
  using iterator_category = typename parent::iterator_category;

public:
  const_view_iterator() : parent(){};
  const_view_iterator(const const_view_iterator & it) = default;
  const_view_iterator(const_view_iterator && it) = default;

  template <typename P, typename... Ns>
  const_view_iterator(P * data, Ns... ns) : parent(data, ns...) {}

  const_view_iterator & operator=(const const_view_iterator & it) = default;
};

template <class R, bool is_tensor_ = is_tensor<R>::value>
struct ConstConverterIteratorHelper {
  static inline auto convert(const view_iterator<R> & it) {
    return const_view_iterator<R>(std::unique_ptr<R>(new R(*it, false)));
  }
};

template <class R> struct ConstConverterIteratorHelper<R, false> {
  static inline auto convert(const view_iterator<R> & it) {
    return const_view_iterator<R>(it.data());
  }
};

template <typename R>
class view_iterator
    : public detail::internal_view_iterator<R, view_iterator<R>> {
public:
  using parent = detail::internal_view_iterator<R, view_iterator>;
  using value_type = typename parent::value_type;
  using pointer = typename parent::pointer;
  using reference = typename parent::reference;
  using difference_type = typename parent::difference_type;
  using iterator_category = typename parent::iterator_category;

public:
  view_iterator() : parent(){};
  view_iterator(const view_iterator & it) = default;
  view_iterator(view_iterator && it) = default;

  template <typename P, typename... Ns>
  view_iterator(P * data, Ns... ns) : parent(data, ns...) {}

  view_iterator & operator=(const view_iterator & it) = default;

  operator const_view_iterator<R>() {
    return ConstConverterIteratorHelper<R>::convert(*this);
  }
};

namespace {
  template <std::size_t dim, typename T> struct ViewIteratorHelper {
    using type = TensorBase<T, dim, true>;
  };

  template <typename T> struct ViewIteratorHelper<0, T> { using type = T; };
  template <typename T> struct ViewIteratorHelper<1, T> {
    using type = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>>;
  };
  template <typename T> struct ViewIteratorHelper<1, const T> {
    using type = Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>>;
  };

  template <typename T> struct ViewIteratorHelper<2, T> {
    using type = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>;
  };

  template <typename T> struct ViewIteratorHelper<2, const T> {
    using type =
        Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>;
  };

  template <std::size_t dim, typename T>
  using ViewIteratorHelper_t = typename ViewIteratorHelper<dim, T>::type;
} // namespace

// #include <iterator>

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

} // namespace akantu

#endif /* __AKANTU_AKA_VIEW_ITERATORS_HH__ */
