/**
 * @file   aka_enumerate_iterator.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  jeu déc 12 2019
 *
 * @brief implementation of enumerate.
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
#include "iterators/aka_zip_iterator.hh"
/* -------------------------------------------------------------------------- */
#include <iterator>
#include <utility>
/* -------------------------------------------------------------------------- */

#ifndef AKA_ENUMERATE_ITERATOR_HH
#define AKA_ENUMERATE_ITERATOR_HH

#ifndef AKANTU_ITERATORS_NAMESPACE
#define AKANTU_ITERATORS_NAMESPACE akantu
#endif

namespace AKANTU_ITERATORS_NAMESPACE {

/* -------------------------------------------------------------------------- */
namespace iterators {
  template <class Iterator, class size_type_> class EnumerateIterator {
  public:
    using value_type =
        std::tuple<std::size_t,
                   typename std::iterator_traits<Iterator>::value_type>;
    using difference_type = size_type_;
    using pointer =
        std::tuple<size_type_,
                   typename std::iterator_traits<Iterator>::pointer>;
    using reference =
        std::tuple<size_type_,
                   typename std::iterator_traits<Iterator>::reference>;
    using iterator_category = std::input_iterator_tag;

  public:
    explicit EnumerateIterator(Iterator && iterator) : iterator(iterator) {}

    // input iterator ++it
    EnumerateIterator & operator++() {
      ++iterator;
      ++index;
      return *this;
    }

    // input iterator it++
    EnumerateIterator operator++(int) {
      auto cpy = *this;
      this->operator++();
      return cpy;
    }

    // input iterator it != other_it
    bool operator!=(const EnumerateIterator & other) const {
      return iterator != other.iterator;
    }

    // input iterator dereference *it
    decltype(auto) operator*() {
      return std::tuple_cat(std::make_tuple(index), *iterator);
    }

    bool operator==(const EnumerateIterator & other) const {
      return not this->operator!=(other);
    }

  private:
    Iterator iterator;
    size_type_ index{0};
  };

  template <class Iterator, class size_type_ = std::size_t>
  inline constexpr decltype(auto) enumerate(Iterator && iterator,
                                            size_type_ /*size*/) {
    return EnumerateIterator<Iterator, size_type_>(
        std::forward<Iterator>(iterator));
  }

} // namespace iterators

namespace containers {
  template <class... Containers> class EnumerateContainer {
    using ZipContainer_t = ZipContainer<Containers...>;
    using size_type = typename ZipContainer_t::size_type;

  public:
    explicit EnumerateContainer(Containers &&... containers)
        : zip_container(std::forward<Containers>(containers)...) {}

    decltype(auto) begin() {
      return iterators::enumerate(zip_container.begin(), size_type{});
    }

    decltype(auto) begin() const {
      return iterators::enumerate(zip_container.begin(), size_type{});
    }

    decltype(auto) end() {
      return iterators::enumerate(zip_container.end(), size_type{});
    }

    decltype(auto) end() const {
      return iterators::enumerate(zip_container.end(), size_type{});
    }

  private:
    ZipContainer_t zip_container;
  };
} // namespace containers

template <class... Container>
inline constexpr decltype(auto) enumerate(Container &&... container) {
  return containers::EnumerateContainer<Container...>(
      std::forward<Container>(container)...);
}

} // namespace AKANTU_ITERATORS_NAMESPACE

namespace std {
template <class Iterator, class size_type>
struct iterator_traits<::AKANTU_ITERATORS_NAMESPACE::iterators::
                           EnumerateIterator<Iterator, size_type>> {
private:
  using iterator_type =
      typename ::AKANTU_ITERATORS_NAMESPACE::iterators::EnumerateIterator<
          Iterator, size_type>;

public:
  using iterator_category = typename iterator_type::iterator_category;
  using value_type = typename iterator_type::value_type;
  using difference_type = typename iterator_type::difference_type;
  using pointer = typename iterator_type::pointer;
  using reference = typename iterator_type::reference;
};
} // namespace std

#endif /* AKA_ENUMERATE_ITERATOR_HH */
