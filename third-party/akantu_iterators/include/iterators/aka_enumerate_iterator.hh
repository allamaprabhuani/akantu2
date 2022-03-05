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
#include "iterators/aka_iterator_tools.hh"
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
namespace AKA_ITERATOR_EXPORT_NAMESPACE iterators {
  template <class Iterator, class size_type_>
  class EnumerateIterator
      : public details::CopyAssignmentEnabler<
            aka::conjunction<std::is_copy_assignable<Iterator>,
                             std::is_copy_constructible<Iterator>>::value>,
        public details::MoveAssignmentEnabler<
            aka::conjunction<std::is_move_assignable<Iterator>,
                             std::is_move_constructible<Iterator>>::value> {

  private:
    using it_traits = typename std::iterator_traits<Iterator>;

  public:
    using value_type = std::tuple<size_type_, typename it_traits::value_type>;
    using difference_type = size_type_;
    using pointer = std::tuple<size_type_, typename it_traits::pointer>;
    using reference = std::tuple<size_type_, typename it_traits::reference>;
    using iterator_category = typename it_traits::iterator_category;

    explicit EnumerateIterator(Iterator iterator)
        : iterator(std::move(iterator)) {}

    // input iterator ++it
    auto operator++() -> EnumerateIterator & {
      ++iterator;
      ++index;
      return *this;
    }

    // input iterator it++
    auto operator++(int) -> EnumerateIterator {
      auto cpy = *this;
      this->operator++();
      return cpy;
    }

    // input iterator it != other_it
    auto operator!=(const EnumerateIterator & other) const -> bool {
      return iterator != other.iterator;
    }

    // input iterator dereference *it
    auto operator*() -> decltype(auto) {
      return std::tuple_cat(std::make_tuple(index), *iterator);
    }

    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::bidirectional_iterator_tag>::value> * = nullptr>
    auto operator--() -> EnumerateIterator & {
      --iterator;
      --index;
      return *this;
    }

    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::bidirectional_iterator_tag>::value> * = nullptr>
    auto operator--(int) -> EnumerateIterator {
      auto cpy = *this;
      this->operator--();
      return cpy;
    }

    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::random_access_iterator_tag>::value> * = nullptr>
    auto operator-(const EnumerateIterator & other) -> difference_type {
      return index - other.index;
    }

    // random iterator it[idx]
    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::random_access_iterator_tag>::value> * = nullptr>
    auto operator[](size_type_ idx) -> decltype(auto) {
      return std::tuple_cat(index + idx, iterator[idx]);
    }

    // random iterator it + n
    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::random_access_iterator_tag>::value> * = nullptr>
    auto operator+(size_type_ n) -> decltype(auto) {
      auto && it = EnumerateIterator(iterator + n);
      it.index = this->index + n;
      return it;
    }

    // random iterator it - n
    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::random_access_iterator_tag>::value> * = nullptr>
    auto operator-(size_type_ n) -> decltype(auto) {
      auto && it = EnumerateIterator(iterator - n);
      it.index = this->index - n;
      return it;
    }

    template <
        class iterator_category_ = iterator_category,
        std::enable_if_t<aka::is_iterator_category_at_least<
            iterator_category_, std::forward_iterator_tag>::value> * = nullptr>
    auto operator==(const EnumerateIterator & other) const -> bool {
      return not this->operator!=(other);
    }

  private:
    Iterator iterator;
    size_type_ index{0};
  };

  template <class Iterator, class size_type_ = std::size_t>
  inline constexpr auto enumerate(Iterator && iterator, size_type_ /*size*/)
      -> decltype(auto) {
    return EnumerateIterator<Iterator, size_type_>(
        std::forward<Iterator>(iterator));
  }

} // namespace iterators

namespace containers {
  template <class... Containers>
  class [[gnu::visibility("hidden")]] EnumerateContainer {
    using ZipContainer_t = ZipContainer<Containers...>;
    using size_type = typename ZipContainer_t::size_type;

  public:
    explicit EnumerateContainer(Containers && ... containers)
        : zip_container(std::forward<Containers>(containers)...) {}

    auto begin()->decltype(auto) {
      return iterators::enumerate(zip_container.begin(), size_type{});
    }

    auto begin() const->decltype(auto) {
      return iterators::enumerate(zip_container.begin(), size_type{});
    }

    auto end()->decltype(auto) {
      return iterators::enumerate(zip_container.end(), size_type{});
    }

    auto end() const->decltype(auto) {
      return iterators::enumerate(zip_container.end(), size_type{});
    }

  private:
    ZipContainer_t zip_container;
  };
} // namespace containers

template <class... Container>
inline constexpr auto enumerate(Container &&... container) -> decltype(auto) {
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
