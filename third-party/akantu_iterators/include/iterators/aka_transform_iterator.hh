/**
 * @file   aka_transform_iterator.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  jeu déc 12 2019
 *
 * @brief transform adaptors
 *
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * akantu-iterators is free  software: you can redistribute it and/or  modify it
 * under the terms  of the  GNU Lesser  General Public  License as  published by
 * the Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * akantu-iterators is  distributed in the  hope that it  will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public
 * License  for more details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with akantu-iterators. If not, see <http://www.gnu.org/licenses/>.
 *
 */
/* -------------------------------------------------------------------------- */
#include "aka_iterator_tools.hh"
/* -------------------------------------------------------------------------- */
#include <utility>
/* -------------------------------------------------------------------------- */

#ifndef AKA_TRANSFORM_ITERATOR_H
#define AKA_TRANSFORM_ITERATOR_H

namespace AKANTU_ITERATORS_NAMESPACE {

namespace iterators AKA_ITERATOR_EXPORT_NAMESPACE {
  template <class iterator_t, class operator_t>
  class transform_adaptor_iterator {
  public:
    using value_type = decltype(std::declval<operator_t>()(
        std::declval<typename iterator_t::value_type>()));
    using difference_type = typename iterator_t::difference_type;
    using pointer = std::decay_t<value_type> *;
    using reference = value_type &;
    using iterator_category = typename iterator_t::iterator_category;

    constexpr transform_adaptor_iterator(iterator_t it, operator_t op)
        : it(std::move(it)), op(op) {}
    constexpr transform_adaptor_iterator(const transform_adaptor_iterator &) =
        default;

    constexpr transform_adaptor_iterator & operator++() {
      ++it;
      return *this;
    }

    constexpr decltype(auto) operator*() {
      return op(std::forward<decltype(*it)>(*it));
    }

    constexpr bool operator==(const transform_adaptor_iterator & other) const {
      return (it == other.it);
    }

    constexpr bool operator!=(const transform_adaptor_iterator & other) const {
      return not operator==(other);
    }

    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::random_access_iterator_tag>::value> * = nullptr>
    constexpr difference_type
    operator-(const transform_adaptor_iterator & other) {
      return other - *this;
    }

  private:
    iterator_t it;
    operator_t op;
  };

  template <class iterator_t, class operator_t>
  constexpr auto make_transform_adaptor_iterator(iterator_t it, operator_t op)
      -> decltype(auto) {
    return transform_adaptor_iterator<iterator_t, operator_t>(
        it, std::forward<operator_t>(op));
  }

} // namespace AKA_ITERATOR_EXPORT_NAMESPACE

namespace containers AKA_ITERATOR_EXPORT_NAMESPACE {
  template <class container_t, class operator_t>
  class TransformIteratorAdaptor {
  public:
    // using const_iterator = typename
    // std::decay_t<container_t>::const_iterator; using iterator = typename
    // std::decay_t<container_t>::iterator;
    using size_type = typename std::decay_t<container_t>::size_type;

    constexpr TransformIteratorAdaptor(container_t && cont, operator_t op)
        : cont(std::forward<container_t>(cont)),
          op(std::forward<operator_t>(op)) {}

    constexpr auto begin() const -> decltype(auto) {
      return iterators::make_transform_adaptor_iterator(cont.begin(), op);
    }
    constexpr auto begin() -> decltype(auto) {
      return iterators::make_transform_adaptor_iterator(cont.begin(), op);
    }

    constexpr auto end() const -> decltype(auto) {
      return iterators::make_transform_adaptor_iterator(cont.end(), op);
    }
    constexpr auto end() -> decltype(auto) {
      return iterators::make_transform_adaptor_iterator(cont.end(), op);
    }

  private:
    container_t cont;
    operator_t op;
  };
} // namespace AKA_ITERATOR_EXPORT_NAMESPACE

template <class container_t, class operator_t>
constexpr auto make_transform_adaptor(container_t && cont, operator_t && op)
    -> decltype(auto) {
  return containers::TransformIteratorAdaptor<container_t, operator_t>(
      std::forward<container_t>(cont), std::forward<operator_t>(op));
}

template <class container_t>
constexpr auto make_keys_adaptor(container_t && cont) -> decltype(auto) {
  return make_transform_adaptor(
      std::forward<container_t>(cont),
      [](auto && pair) -> const auto & { return pair.first; });
}

template <class container_t>
constexpr auto make_values_adaptor(container_t && cont) -> decltype(auto) {
  return make_transform_adaptor(
      std::forward<container_t>(cont),
      [](auto && pair) -> auto & { return pair.second; });
}

template <class container_t>
constexpr auto make_dereference_adaptor(container_t && cont) -> decltype(auto) {
  return make_transform_adaptor(
      std::forward<container_t>(cont),
      [](auto && value) -> decltype(*value) { return *value; });
}

namespace details {
  template <typename T> struct BroadcastHelper {
    BroadcastHelper(T && t) : t(std::forward<T>(t)) {}
    constexpr auto operator()() -> std::remove_reference_t<T> & { return t; }
    T t;
  };
} // namespace details

template <typename Type, typename Size>
constexpr auto broadcast(Type && data, Size size) -> decltype(auto) {
  auto && accessor = details::BroadcastHelper<Type>(std::forward<Type>(data));
  return make_transform_adaptor(
      arange(size), [accessor](auto && /*value*/) constexpr mutable
                            ->std::remove_reference_t<Type> &
                        { return accessor(); });
}

} // namespace AKANTU_ITERATORS_NAMESPACE

#endif // AKA_TRANSFORM_ITERATOR_H
