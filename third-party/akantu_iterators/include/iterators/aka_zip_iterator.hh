/**
 * @file   aka_zip_iterator.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  jeu déc 12 2019
 *
 * @brief A Documented file.
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
#include "aka_compatibilty_with_cpp_standard.hh"
#include "aka_iterator_tools.hh"
#include "aka_tuple_tools.hh"
/* -------------------------------------------------------------------------- */
#include <iterator>
#include <utility>
/* -------------------------------------------------------------------------- */

#ifndef AKA_ZIP_ITERATOR_HH
#define AKA_ZIP_ITERATOR_HH

#ifndef AKANTU_ITERATORS_NAMESPACE
#define AKANTU_ITERATORS_NAMESPACE akantu
#endif

namespace AKANTU_ITERATORS_NAMESPACE {

/* -------------------------------------------------------------------------- */
namespace iterators {
  /* ------------------------------------------------------------------------ */
  template <template <class...> class Tuple, class... Iterators>
  class ZipIterator_
      : public details::CopyAssignmentEnabler<
            aka::conjunction<std::is_copy_assignable<Iterators>...,
                             std::is_copy_constructible<Iterators>...>::value>,
        public details::MoveAssignmentEnabler<
            aka::conjunction<std::is_move_assignable<Iterators>...,
                             std::is_move_constructible<Iterators>...>::value> {
  private:
    using tuple_t = Tuple<Iterators...>;

  public:
    using value_type =
        Tuple<typename std::iterator_traits<Iterators>::value_type...>;
    using difference_type = std::common_type_t<
        typename std::iterator_traits<Iterators>::difference_type...>;
    using pointer = Tuple<typename std::iterator_traits<Iterators>::pointer...>;
    using reference =
        Tuple<typename std::iterator_traits<Iterators>::reference...>;
    using iterator_category = // std::input_iterator_tag;
        std::common_type_t<
            typename std::iterator_traits<Iterators>::iterator_category...>;

    //    using nb_iterators = sizeof...(Iterators);

  public:
    explicit ZipIterator_(tuple_t iterators)
        : iterators(std::move(iterators)) {}

    /* ---------------------------------------------------------------------- */
    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::bidirectional_iterator_tag>::value> * = nullptr>
    ZipIterator_ & operator--() {
      tuple::foreach ([](auto && it) { --it; }, iterators);
      return *this;
    }

    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::bidirectional_iterator_tag>::value> * = nullptr>
    ZipIterator_ operator--(int) {
      auto cpy = *this;
      this->operator--();
      return cpy;
    }

    // input iterator ++it
    ZipIterator_ & operator++() {
      tuple::foreach ([](auto && it) { ++it; }, iterators);
      return *this;
    }

    // input iterator it++
    ZipIterator_ operator++(int) {
      auto cpy = *this;
      this->operator++();
      return cpy;
    }

    // input iterator it != other_it
    bool operator!=(const ZipIterator_ & other) const {
      // return tuple::are_not_equal(iterators, other.iterators);
      return std::get<0>(iterators) !=
             std::get<0>(other.iterators); // helps the compiler to optimize
    }

    // input iterator dereference *it
    decltype(auto) operator*() {
      return tuple::transform([](auto && it) -> decltype(auto) { return *it; },
                              iterators);
    }

    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::random_access_iterator_tag>::value> * = nullptr>
    difference_type operator-(const ZipIterator_ & other) {
      return std::get<0>(this->iterators) - std::get<0>(other.iterators);
    }

    // random iterator it[idx]
    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::random_access_iterator_tag>::value> * = nullptr>
    decltype(auto) operator[](std::size_t idx) {
      return tuple::transform(
          [idx](auto && it) -> decltype(auto) { return it[idx]; }, iterators);
    }

    // random iterator it + n
    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::random_access_iterator_tag>::value> * = nullptr>
    decltype(auto) operator+(std::size_t n) {
      return ZipIterator_(std::forward<tuple_t>(tuple::transform(
          [n](auto && it) -> decltype(auto) { return it + n; }, iterators)));
    }

    // random iterator it - n
    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::random_access_iterator_tag>::value> * = nullptr>
    decltype(auto) operator-(std::size_t n) {
      return ZipIterator_(std::forward<tuple_t>(tuple::transform(
          [n](auto && it) -> decltype(auto) { return it - n; }, iterators)));
    }

    template <
        class iterator_category_ = iterator_category,
        std::enable_if_t<aka::is_iterator_category_at_least<
            iterator_category_, std::forward_iterator_tag>::value> * = nullptr>
    bool operator==(const ZipIterator_ & other) const {
      return not tuple::are_not_equal(iterators, other.iterators);
    }

  private:
    tuple_t iterators;
  };

  template <class... Iterators>
  using ZipIterator = ZipIterator_<std::tuple, Iterators...>;

  template <class... Iterators>
  using NamedZipIterator = ZipIterator_<tuple::named_tuple, Iterators...>;
} // namespace iterators

/* -------------------------------------------------------------------------- */
template <class... Iterators>
decltype(auto) zip_iterator(std::tuple<Iterators...> && iterators_tuple) {
  auto zip = iterators::ZipIterator<Iterators...>(
      std::forward<decltype(iterators_tuple)>(iterators_tuple));
  return zip;
}

template <class... Iterators>
decltype(auto)
zip_iterator(tuple::named_tuple<Iterators...> && iterators_tuple) {
  auto zip = iterators::NamedZipIterator<Iterators...>(
      std::forward<decltype(iterators_tuple)>(iterators_tuple));
  return zip;
}

/* -------------------------------------------------------------------------- */
namespace containers {
  template <template <class...> class Tuple, class... Containers>
  class ZipContainer_ {
  public:
    using containers_t = Tuple<Containers...>;
    using size_type = std::common_type_t<aka::size_type_t<Containers>...>;

    explicit ZipContainer_(Containers &&... containers)
        : containers(std::forward<Containers>(containers)...) {}

    decltype(auto) begin() const {
      return zip_iterator(
          tuple::transform([](auto && c) { return c.begin(); },
                           std::forward<containers_t>(containers)));
    }

    decltype(auto) end() const {
      return zip_iterator(
          tuple::transform([](auto && c) { return c.end(); },
                           std::forward<containers_t>(containers)));
    }

    decltype(auto) begin() {
      return zip_iterator(
          tuple::transform([](auto && c) { return c.begin(); },
                           std::forward<containers_t>(containers)));
    }

    decltype(auto) end() {
      return zip_iterator(
          tuple::transform([](auto && c) { return c.end(); },
                           std::forward<containers_t>(containers)));
    }

  protected:
    template <class... OtherContainers, std::size_t... Is>
    decltype(auto) append(std::index_sequence<Is...> && /*unused*/,
                          OtherContainers &&... other_containers) {
      return ZipContainer_<Tuple, Containers..., OtherContainers...>(
          std::get<Is>(std::forward<containers_t>(containers))...,
          std::forward<OtherContainers>(
              std::forward<OtherContainers>(other_containers))...);
    }

    template <class OtherContainers, std::size_t... Is>
    decltype(auto) extend(std::index_sequence<Is...> && /*unused*/,
                          OtherContainers && other_containers) {
      return append(
          std::make_index_sequence<sizeof...(Containers)>{},
          std::get<Is>(std::forward<containers_t>(other_containers))...);
    }

    // template <std::size_t nth, std::size_t... Is_before,
    //           std::size_t... Is_after>
    // decltype(auto) remove(std::index_sequence<Is_before...> && /*unused*/,
    //                       std::index_sequence<Is_after...> && /*unused*/) {
    //   using tuple::tuple_element_t;
    //   return ZipContainer_<
    //       Tuple, tuple_element_t<Is_before, containers_t>...,
    //       tuple_element_t<Is_after + nth + 1, containers_t>...>(
    //       std::forward<tuple_element_t<Is_before, containers_t>>(
    //           std::get<Is_before>(containers))...,
    //       std::forward<tuple_element_t<Is_after + nth + 1, containers_t>>(
    //           std::get<Is_after + nth + 1>(containers))...);
    // }

  public:
    template <class... OtherContainers>
    decltype(auto) append(OtherContainers &&... other_containers) {
      static_assert(
          tuple::is_named_tuple<containers_t>::value ==
              aka::conjunction<
                  tuple::is_named_tag<std::decay_t<OtherContainers>>...>::value,
          "Cannot append named tag in non named zip or vice "
          "versa");
      return append(std::make_index_sequence<sizeof...(Containers)>{},
                    std::forward<OtherContainers>(other_containers)...);
    }

    template <class... OtherContainers>
    decltype(auto)
    extend(const ZipContainer_<Tuple, OtherContainers...> & other) {
      return extend(std::make_index_sequence<sizeof...(OtherContainers)>{},
                    other.containers);
    }

    // template <size_t Tag, class Containers_t = containers_t,
    //           std::enable_if_t<tuple::is_named_tuple<Containers_t>::value> *
    //           =
    //               nullptr>
    // decltype(auto) remove() {
    //   constexpr auto nth =
    //       containers.template get_element_index(tuple::get<Tag>());
    //   return remove<nth>(
    //       std::make_index_sequence<nth>{},
    //       std::make_index_sequence<sizeof...(Containers) - nth - 1>{});
    // }

    template <size_t Tag> decltype(auto) remove();

    template <size_t Tag, class Container>
    decltype(auto) replace(Container && cont);

  private:
    template <template <class...> class OtherTuple, class... OtherContainers>
    friend class ZipContainer_;

    containers_t containers;
  };

  template <class... Containers>
  using ZipContainer = ZipContainer_<std::tuple, Containers...>;

  template <class... Containers>
  using NamedZipContainer = ZipContainer_<tuple::named_tuple, Containers...>;
} // namespace containers

/* -------------------------------------------------------------------------- */
template <class... Containers,
          std::enable_if_t<not aka::conjunction<tuple::is_named_tag<
              std::decay_t<Containers>>...>::value> * = nullptr>
decltype(auto) zip(Containers &&... conts) {
  return containers::ZipContainer<Containers...>(
      std::forward<Containers>(conts)...);
}

template <class... Containers,
          std::enable_if_t<aka::conjunction<tuple::is_named_tag<
              std::decay_t<Containers>>...>::value> * = nullptr>
decltype(auto) zip(Containers &&... conts) {
  return containers::NamedZipContainer<Containers...>(
      std::forward<Containers>(conts)...);
}

/* -------------------------------------------------------------------------- */
template <class zip_container_1_t, class zip_container_2_t>
decltype(auto) zip_cat(zip_container_1_t && cont1, zip_container_2_t && cont2) {
  return std::forward<zip_container_1_t>(cont1).extend(
      std::forward<zip_container_2_t>(cont2));
}

template <class zip_container_t, class... container_t>
decltype(auto) zip_append(zip_container_t && zip_container,
                          container_t &&... cont) {
  return std::forward<zip_container_t>(zip_container)
      .append(std::forward<container_t>(cont)...);
}

template <size_t Tag, class zip_container_t, class container_t>
decltype(auto) zip_replace(zip_container_t && zip_container,
                           container_t && cont) {
  return std::forward<zip_container_t>(zip_container)
      .template replace<Tag>(std::forward<container_t>(cont));
}

template <size_t Tag, class zip_container_t>
decltype(auto) zip_remove(zip_container_t && zip_container) {
  return std::forward<zip_container_t>(zip_container).template remove<Tag>();
}

namespace details {
  template <class Tuple, std::size_t... Is,
            std::enable_if_t<
                tuple::is_named_tuple<std::decay_t<Tuple>>::value> * = nullptr>
  decltype(auto)
  make_zip_from_tuple_impl(std::index_sequence<Is...> && /*unused*/,
                           Tuple && tuple) {
    return zip(
        std::get<tuple::tuple_name_tag<Is, std::decay_t<Tuple>>>() =
            std::forward<tuple::tuple_element_t<Is, std::decay_t<Tuple>>>(
                std::get<Is>(tuple))...);
  }
} // namespace details

template <class Tuple, std::enable_if_t<tuple::is_named_tuple<
                           std::decay_t<Tuple>>::value> * = nullptr>
decltype(auto) make_zip_from_tuple(Tuple && tuple) {
  return details::make_zip_from_tuple_impl(
      std::make_index_sequence<
          std::tuple_size<typename std::decay_t<Tuple>::parent>::value>{},
      std::forward<Tuple>(tuple));
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
namespace containers {
  template <template <class...> class Tuple, class... Containers>
  template <std::size_t nth>
  decltype(auto) ZipContainer_<Tuple, Containers...>::remove() {
    return make_zip_from_tuple(
        tuple::remove<nth>(std::forward<containers_t>(containers)));
  }

  template <template <class...> class Tuple, class... Containers>
  template <std::size_t nth, class Container>
  decltype(auto)
  ZipContainer_<Tuple, Containers...>::replace(Container && cont) {
    return make_zip_from_tuple(tuple::replace<nth>(
        std::forward<containers_t>(containers), std::forward<Container>(cont)));
  }

} // namespace containers
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

} // namespace AKANTU_ITERATORS_NAMESPACE

namespace std {
template <template <class...> class Tuple, typename... Its>
struct iterator_traits<
    ::AKANTU_ITERATORS_NAMESPACE::iterators::ZipIterator_<Tuple, Its...>> {
private:
  using iterator_type =
      typename ::AKANTU_ITERATORS_NAMESPACE::iterators::ZipIterator_<Tuple,
                                                                     Its...>;

public:
  using iterator_category = typename iterator_type::iterator_category;
  using value_type = typename iterator_type::value_type;
  using difference_type = typename iterator_type::difference_type;
  using pointer = typename iterator_type::pointer;
  using reference = typename iterator_type::reference;
};
} // namespace std

#endif /* AKA_ZIP_ITERATOR_HH */
