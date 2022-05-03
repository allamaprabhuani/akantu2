/**
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
 */

/* -------------------------------------------------------------------------- */
#include "iterators/aka_iterator_tools.hh"
/* -------------------------------------------------------------------------- */
#include <cmath>
#include <utility>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_AKA_REPEAT_ITERATOR_HH
#define AKANTU_AKA_REPEAT_ITERATOR_HH
namespace AKANTU_ITERATORS_NAMESPACE {

namespace AKA_ITERATOR_EXPORT_NAMESPACE iterators {

  template <class Iterator, typename N>
  class RepeatNIterator
      : public details::CopyAssignmentEnabler<
            aka::conjunction<std::is_copy_assignable<Iterator>,
                             std::is_copy_constructible<Iterator>>::value>,
        public details::MoveAssignmentEnabler<
            aka::conjunction<std::is_move_assignable<Iterator>,
                             std::is_move_constructible<Iterator>>::value> {

  private:
    using it_traits = typename std::iterator_traits<Iterator>;

  public:
    using value_type = typename it_traits::value_type;
    using difference_type = typename it_traits::difference_type;
    using pointer = typename it_traits::pointer;
    using reference = typename it_traits::reference;
    using iterator_category = typename it_traits::iterator_category;

    explicit RepeatNIterator(Iterator iterator, N n)
        : iterator(std::move(iterator)), n(std::move(n)) {}

    // input iterator ++it
    auto operator++() -> RepeatNIterator & {
      ++i;
      if (i == n) {
        ++iterator;
        i = 0;
      }
      return *this;
    }

    // input iterator it++
    auto operator++(int) -> RepeatNIterator {
      auto cpy = *this;
      this->operator++();
      return cpy;
    }

    // input iterator it != other_it
    auto operator!=(const RepeatNIterator & other) const -> bool {
      return iterator != other.iterator or i != other.i;
    }

    // input iterator dereference *it
    auto operator*() -> decltype(auto) { return *iterator; }

    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::bidirectional_iterator_tag>::value> * = nullptr>
    auto operator--() -> RepeatNIterator & {
      if (i == 0) {
        --iterator;
        i = n - 1;
      } else {
        --i;
      }
      return *this;
    }

    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::bidirectional_iterator_tag>::value> * = nullptr>
    auto operator--(int) -> RepeatNIterator {
      auto cpy = *this;
      this->operator--();
      return cpy;
    }

    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::random_access_iterator_tag>::value> * = nullptr>
    auto operator-(const RepeatNIterator & other) const -> difference_type {
      return (iterator - other.iterator) + (i - other.i);
    }

    // random iterator it[idx]
    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::random_access_iterator_tag>::value> * = nullptr>
    auto operator[](const N & idx) -> decltype(auto) {
      return iterator[(i + idx) / n];
    }

    // random iterator it + n
    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::random_access_iterator_tag>::value> * = nullptr>
    auto operator+(const N & n_) const -> decltype(auto) {
      auto it = RepeatNIterator(iterator + n_ / n, n);
      it.i = (i + n_) % n;
      return it;
    }

    // random iterator it - n
    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::random_access_iterator_tag>::value> * = nullptr>
    auto operator-(const N & n_) const -> decltype(auto) {
      auto it = RepeatNIterator(iterator + N(std::floor(double(-n_) / n)), n);
      it.i = (i + n - n_) % n;
      return it;
    }

    template <
        class iterator_category_ = iterator_category,
        std::enable_if_t<aka::is_iterator_category_at_least<
            iterator_category_, std::forward_iterator_tag>::value> * = nullptr>
    auto operator==(const RepeatNIterator & other) const -> bool {
      return not this->operator!=(other);
    }

  private:
    Iterator iterator;
    const N n;
    N i{0};
  };

  template <class Iterator, typename N>
  inline constexpr auto repeat_n(Iterator && iterator, N n) -> decltype(auto) {
    return RepeatNIterator<Iterator, N>(std::forward<Iterator>(iterator),
                                        std::move(n));
  }

} // namespace iterators

namespace AKA_ITERATOR_EXPORT_NAMESPACE containers {
  template <class Container, typename N> class RepeatNContainer {
  public:
    using size_type = typename std::decay_t<Container>::size_type;

    explicit RepeatNContainer(Container && container, N n)
        : container(std::forward<Container>(container)), n(std::move(n)) {}
    auto begin() -> decltype(auto) {
      return iterators::repeat_n(container.begin(), n);
    }

    auto begin() const -> decltype(auto) {
      return iterators::repeat_n(container.begin(), n);
    }

    auto end() -> decltype(auto) {
      return iterators::repeat_n(container.end(), n);
    }

    auto end() const -> decltype(auto) {
      return iterators::repeat_n(container.end(), n);
    }

  private:
    Container container;
    N n;
  };
} // namespace containers

template <class Container, typename N>
inline constexpr auto repeat_n(Container && container, N n) -> decltype(auto) {
  return containers::RepeatNContainer<Container, N>(
      std::forward<Container>(container), std::move(n));
}

} // namespace AKANTU_ITERATORS_NAMESPACE

namespace std {
template <class Iterator, typename N>
struct iterator_traits<
    ::AKANTU_ITERATORS_NAMESPACE::iterators::RepeatNIterator<Iterator, N>> {
private:
  using iterator_type =
      typename ::AKANTU_ITERATORS_NAMESPACE::iterators::RepeatNIterator<
          Iterator, N>;

public:
  using iterator_category = typename iterator_type::iterator_category;
  using value_type = typename iterator_type::value_type;
  using difference_type = typename iterator_type::difference_type;
  using pointer = typename iterator_type::pointer;
  using reference = typename iterator_type::reference;
};

} // namespace std

#endif /* AKANTU_AKA_REPEAT_ITERATOR_HH */
