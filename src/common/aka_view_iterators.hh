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
//#include "aka_types.hh"
/* -------------------------------------------------------------------------- */
#include <memory>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_AKA_VIEW_ITERATORS_HH__
#define __AKANTU_AKA_VIEW_ITERATORS_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
/* Iterators                                                                  */
/* -------------------------------------------------------------------------- */
namespace details {
  template <class R, class daughter, class IR = R,
            bool is_tensor = is_tensor<R>::value>
  class internal_view_iterator {
  public:
    using value_type = R;
    using pointer = R *;
    using reference = R &;
    using proxy = typename R::proxy;
    using const_proxy = const typename R::proxy;
    using const_reference = const R &;
    using internal_value_type = IR;
    using internal_pointer = IR *;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::random_access_iterator_tag;

  public:
    internal_view_iterator() = default;

    internal_view_iterator(pointer data, UInt _offset)
        : _offset(_offset), initial(data), ret(nullptr), ret_ptr(data) {
      AKANTU_ERROR(
          "The constructor should never be called it is just an ugly trick...");
    }

    internal_view_iterator(std::unique_ptr<internal_value_type> && wrapped)
        : _offset(wrapped->size()), initial(wrapped->storage()),
          ret(std::move(wrapped)), ret_ptr(ret->storage()) {}

    internal_view_iterator(const internal_view_iterator & it) {
      if (this != &it) {
        this->_offset = it._offset;
        this->initial = it.initial;
        this->ret_ptr = it.ret_ptr;
        this->ret = std::make_unique<internal_value_type>(*it.ret, false);
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
        if (this->ret)
          this->ret->shallowCopy(*it.ret);
        else
          this->ret = std::make_unique<internal_value_type>(*it.ret, false);
      }
      return *this;
    }

    UInt getCurrentIndex() {
      return (this->ret_ptr - this->initial) / this->_offset;
    };

    inline reference operator*() {
      ret->values = ret_ptr;
      return *ret;
    };
    inline const_reference operator*() const {
      ret->values = ret_ptr;
      return *ret;
    };
    inline pointer operator->() {
      ret->values = ret_ptr;
      return ret.get();
    };
    inline daughter & operator++() {
      ret_ptr += _offset;
      return static_cast<daughter &>(*this);
    };
    inline daughter & operator--() {
      ret_ptr -= _offset;
      return static_cast<daughter &>(*this);
    };

    inline daughter & operator+=(const UInt n) {
      ret_ptr += _offset * n;
      return static_cast<daughter &>(*this);
    }
    inline daughter & operator-=(const UInt n) {
      ret_ptr -= _offset * n;
      return static_cast<daughter &>(*this);
    }

    inline proxy operator[](const UInt n) {
      ret->values = ret_ptr + n * _offset;
      return proxy(*ret);
    }
    inline const_proxy operator[](const UInt n) const {
      ret->values = ret_ptr + n * _offset;
      return const_proxy(*ret);
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

    inline pointer data() const { return ret_ptr; }
    inline difference_type offset() const { return _offset; }

  protected:
    UInt _offset{0};
    pointer initial{nullptr};
    std::unique_ptr<internal_value_type> ret{nullptr};
    pointer ret_ptr{nullptr};
  };

  /* --------------------------------------------------------------------------
   */
  /**
   * Specialization for scalar types
   */
  template <class R, class daughter, class IR>
  class internal_view_iterator<R, daughter, IR, false> {
  public:
    using value_type = R;
    using pointer = R *;
    using reference = R &;
    using const_reference = const R &;
    using internal_value_type = IR;
    using internal_pointer = IR *;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::random_access_iterator_tag;

  public:
    internal_view_iterator(pointer data = nullptr) : ret(data), initial(data){};
    internal_view_iterator(const internal_view_iterator & it) = default;
    internal_view_iterator(internal_view_iterator && it) = default;

    virtual ~internal_view_iterator() = default;

    inline internal_view_iterator &
    operator=(const internal_view_iterator & it) = default;

    UInt getCurrentIndex() { return (this->ret - this->initial); };

    inline reference operator*() { return *ret; };
    inline const_reference operator*() const { return *ret; };
    inline pointer operator->() { return ret; };
    inline daughter & operator++() {
      ++ret;
      return static_cast<daughter &>(*this);
    };
    inline daughter & operator--() {
      --ret;
      return static_cast<daughter &>(*this);
    };

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
} // namespace details

/* -------------------------------------------------------------------------- */
template <typename R>
class view_iterator;

template <typename R>
class const_view_iterator
    : public details::internal_view_iterator<const R, const_view_iterator<R>,
                                             R> {
public:
  using parent =
      details::internal_view_iterator<const R, const_view_iterator, R>;
  using value_type = typename parent::value_type;
  using pointer = typename parent::pointer;
  using reference = typename parent::reference;
  using difference_type = typename parent::difference_type;
  using iterator_category = typename parent::iterator_category;

public:
  const_view_iterator() : parent(){};
  // const_view_iterator(pointer data, UInt offset) : parent(data, offset)
  // {} const_view_iterator(pointer warped) : parent(warped) {}
  // const_view_iterator(const parent & it) : parent(it) {}

  const_view_iterator(const const_view_iterator & it) = default;
  const_view_iterator(const_view_iterator && it) = default;

  template <typename P, typename = std::enable_if_t<not is_tensor<P>::value>>
  const_view_iterator(P * data) : parent(data) {}

  template <typename UP_P, typename = std::enable_if_t<
                               is_tensor<typename UP_P::element_type>::value>>
  const_view_iterator(UP_P && tensor) : parent(std::forward<UP_P>(tensor)) {}

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
class view_iterator : public details::internal_view_iterator<R, view_iterator<R>> {
public:
  using parent = details::internal_view_iterator<R, view_iterator>;
  using value_type = typename parent::value_type;
  using pointer = typename parent::pointer;
  using reference = typename parent::reference;
  using difference_type = typename parent::difference_type;
  using iterator_category = typename parent::iterator_category;

public:
  view_iterator() : parent(){};
  view_iterator(const view_iterator & it) = default;
  view_iterator(view_iterator && it) = default;

  template <typename P, typename = std::enable_if_t<not is_tensor<P>::value>>
  view_iterator(P * data) : parent(data) {}

  template <typename UP_P, typename = std::enable_if_t<
                               is_tensor<typename UP_P::element_type>::value>>
  view_iterator(UP_P && tensor) : parent(std::forward<UP_P>(tensor)) {}

  view_iterator & operator=(const view_iterator & it) = default;

  operator const_view_iterator<R>() {
    return ConstConverterIteratorHelper<R>::convert(*this);
  }
};

} // namespace akantu

#endif /* __AKANTU_AKA_VIEW_ITERATORS_HH__ */
