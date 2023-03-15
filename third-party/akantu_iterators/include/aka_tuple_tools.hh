/**
 * @file   aka_tuple_tools.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Aug 11 2017
 * @date last modification: Mon Jan 29 2018
 *
 * @brief  iterator interfaces
 *
 *
 * Copyright 2019 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_compatibilty_with_cpp_standard.hh"
#include "aka_named_tuple.hh"
/* -------------------------------------------------------------------------- */
#include <tuple>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_AKA_TUPLE_TOOLS_HH
#define AKANTU_AKA_TUPLE_TOOLS_HH

#ifndef AKANTU_ITERATORS_NAMESPACE
#define AKANTU_ITERATORS_NAMESPACE akantu
#endif

#define AKA_ITERATOR_EXPORT_NAMESPACE __attribute__((visibility("hidden")))
#define FWD(x) std::forward<decltype(x)>(x)

namespace AKANTU_ITERATORS_NAMESPACE {

namespace tuple {
  template <typename T> struct is_tuple : public std::false_type {};
  template <typename... Params>
  struct is_tuple<std::tuple<Params...>> : public std::true_type {};

  /* ------------------------------------------------------------------------ */
  namespace details AKA_ITERATOR_EXPORT_NAMESPACE {
    template <std::size_t N> struct Foreach {
      template <class Tuple>

      static constexpr inline auto not_equal(Tuple && a, Tuple && b) -> bool {
        if (std::get<N - 1>(std::forward<Tuple>(a)) ==
            std::get<N - 1>(std::forward<Tuple>(b))) {
          return false;
        }
        return Foreach<N - 1>::not_equal(std::forward<Tuple>(a),
                                         std::forward<Tuple>(b));
      }

      template <class Tuple, class V>
      static constexpr inline auto find(Tuple && tuple, V && value)
          -> std::size_t {
        constexpr auto size = std::tuple_size<std::decay_t<Tuple>>::value;
        if (std::get<size - N>(std::forward<Tuple>(tuple)) == value) {
          return size - N;
        }

        return Foreach<N - 1>::find(std::forward<Tuple>(tuple),
                                    std::forward<V>(value));
      }
    };

    /* ---------------------------------------------------------------------- */
    template <> struct Foreach<0> {
      template <class Tuple>
      static constexpr inline auto not_equal(Tuple && a, Tuple && b) -> bool {
        return std::get<0>(std::forward<Tuple>(a)) !=
               std::get<0>(std::forward<Tuple>(b));
      }

      template <class Tuple, class V>
      static constexpr inline auto find(Tuple && tuple, V && value)
          -> std::size_t {
        constexpr auto size = std::tuple_size<std::decay_t<Tuple>>::value;
        if (std::get<size - 1>(std::forward<Tuple>(tuple)) == value) {
          return size - 1;
        }

        return size;
      }
    };

    template <typename... Ts>
    auto make_tuple_no_decay(Ts &&... args) -> decltype(auto) {
      return std::tuple<Ts...>(std::forward<Ts>(args)...);
    }

    template <typename... Names, typename... Ts>
    auto make_named_tuple_no_decay(std::tuple<Names...> /*unused*/,
                                   Ts &&... args) -> decltype(auto) {
      return named_tuple<named_tag<Names, Ts>...>(std::forward<Ts>(args)...);
    }

    template <class F, class Tuple, std::size_t... Is>
    constexpr void foreach_impl(F && func, Tuple && tuple,
                                std::index_sequence<Is...> && /*unused*/) {
      (void)std::initializer_list<int>{
          (std::forward<F>(func)(std::get<Is>(std::forward<Tuple>(tuple))),
           0)...};
    }

    template <class F, class Tuple, std::size_t... Is>
    constexpr auto transform_impl(F && func, Tuple && tuple,
                                  std::index_sequence<Is...> && /*unused*/)
        -> decltype(auto) {
      return make_tuple_no_decay(
          std::forward<F>(func)(std::get<Is>(std::forward<Tuple>(tuple)))...);
    }

    template <class F, class Tuple, std::size_t... Is>
    constexpr auto transform_impl(F && func, Tuple && tuple1, Tuple && tuple2,
                                  std::index_sequence<Is...> && /*unused*/)
        -> decltype(auto) {
      return make_tuple_no_decay(
          std::forward<F>(func)(std::get<Is>(std::forward<Tuple>(tuple1)),
                                std::get<Is>(std::forward<Tuple>(tuple2)))...);
    }

    template <class F, class Tuple, std::size_t... Is>
    constexpr auto
    transform_named_impl(F && func, Tuple && tuple,
                         std::index_sequence<Is...> && /*unused*/)
        -> decltype(auto) {
      return make_named_tuple_no_decay(
          typename std::decay_t<Tuple>::Names_t{},
          std::forward<F>(func)(std::get<Is>(std::forward<Tuple>(tuple)))...);
    }

    template <class Tuple, std::size_t... Is>
    auto flatten(Tuple && tuples, std::index_sequence<Is...> /*unused*/)
        -> decltype(auto) {
      return std::tuple_cat(std::get<Is>(tuples)...);
    }

    template <class Tuple, class... Vals, std::size_t... Is>
    auto append_impl(std::index_sequence<Is...> && /*unused*/, Tuple && tuple,
                     Vals &&... other_vals) -> decltype(auto) {
      return make_tuple_no_decay(
          FWD(std::get<Is>(std::forward<Tuple>(tuple)))...,
          std::forward<Vals>(other_vals)...);
    }

    template <class Tuple, class... Vals, std::size_t... Is>
    auto append_named_impl(std::index_sequence<Is...> && /*unused*/,
                           Tuple && tuple, Vals &&... other_vals)
        -> decltype(auto) {
      return make_named_tuple_no_decay(
          std::tuple<tuple_name_tag_t<Is, Tuple>...,
                     typename std::decay_t<Vals>::_tag...>(),
          FWD(std::get<Is>(std::forward<Tuple>(tuple)))...,
          FWD(other_vals._value)...);
    }

    template <std::size_t nth, class Tuple, std::size_t... Is_before,
              std::size_t... Is_after>
    auto remove_impl(std::index_sequence<Is_before...> && /*unused*/,
                     std::index_sequence<Is_after...> && /*unused*/,
                     Tuple && tuple) -> decltype(auto) {
      return make_tuple_no_decay(
          FWD(std::get<Is_before>(std::forward<Tuple>(tuple)))...,
          FWD(std::get<Is_after + nth + 1>(std::forward<Tuple>(tuple)))...);
    }

    template <std::size_t nth, class Tuple, std::size_t... Is_before,
              std::size_t... Is_after>
    auto remove_named_impl(std::index_sequence<Is_before...> && /*unused*/,
                           std::index_sequence<Is_after...> && /*unused*/,
                           Tuple && tuple) -> decltype(auto) {
      return make_named_tuple_no_decay(
          std::tuple<tuple::tuple_name_tag_t<Is_before, std::decay_t<Tuple>>...,
                     tuple::tuple_name_tag_t<Is_after + nth + 1,
                                             std::decay_t<Tuple>>...>{},
          std::get<Is_before>(std::forward<Tuple>(tuple))...,
          std::get<Is_after + nth + 1>(std::forward<Tuple>(tuple))...);
    }

    template <std::size_t nth, class Tuple, class Value,
              std::size_t... Is_before, std::size_t... Is_after>
    auto replace_impl(std::index_sequence<Is_before...> && /*unused*/,
                      std::index_sequence<Is_after...> && /*unused*/,
                      Tuple && tuple, Value && value) -> decltype(auto) {
      return make_tuple_no_decay(
          std::get<Is_before>(std::forward<Tuple>(tuple))...,
          std::forward<Value>(value),
          std::get<Is_after + nth + 1>(std::forward<Tuple>(tuple))...);
    }

    template <std::size_t nth, typename Tag, class Tuple, class Value,
              std::size_t... Is_before, std::size_t... Is_after>
    auto replace_named_impl(std::index_sequence<Is_before...> && /*unused*/,
                            std::index_sequence<Is_after...> && /*unused*/,
                            Tuple && tuple, Value && value) -> decltype(auto) {
      using tuple_type = std::decay_t<Tuple>;
      return make_named_tuple_no_decay(
          std::tuple<
              tuple::tuple_name_tag_t<Is_before, tuple_type>..., Tag,
              tuple::tuple_name_tag_t<Is_after + nth + 1, tuple_type>...>{},
          std::get<Is_before>(std::forward<Tuple>(tuple))...,
          std::forward<Value>(value),
          std::get<Is_after + nth + 1>(std::forward<Tuple>(tuple))...);
    }
  } // namespace AKA_ITERATOR_EXPORT_NAMESPACE

  /* ------------------------------------------------------------------------
   */
  template <class Tuple,
            std::enable_if_t<is_tuple<std::decay_t<Tuple>>::value> * = nullptr>
  constexpr auto are_not_equal(Tuple && a, Tuple && b) -> bool {
    return details::Foreach<std::tuple_size<std::decay_t<Tuple>>::value>::
        not_equal(std::forward<Tuple>(a), std::forward<Tuple>(b));
  }

  template <
      class Tuple,
      std::enable_if_t<is_named_tuple<std::decay_t<Tuple>>::value> * = nullptr>
  constexpr auto are_not_equal(Tuple && a, Tuple && b) -> bool {
    return details::Foreach<
        std::tuple_size<typename std::decay_t<Tuple>::parent>::value>::
        not_equal(std::forward<Tuple>(a), std::forward<Tuple>(b));
  }

  template <class Tuple, class V>
  constexpr decltype(auto) find(Tuple && tuple, V && value) {
    return details::Foreach<std::tuple_size<std::decay_t<Tuple>>::value>::find(
        std::forward<Tuple>(tuple), std::forward<V>(value));
  }

  template <class F, class Tuple,
            std::enable_if_t<is_tuple<std::decay_t<Tuple>>::value> * = nullptr>
  constexpr void foreach (F && func, Tuple && tuple) {
    return details::foreach_impl(
        std::forward<F>(func), std::forward<Tuple>(tuple),
        std::make_index_sequence<
            std::tuple_size<std::decay_t<Tuple>>::value>{});
  }

  template <
      class F, class Tuple,
      std::enable_if_t<is_named_tuple<std::decay_t<Tuple>>::value> * = nullptr>
  constexpr void foreach (F && func, Tuple && tuple) {
    return details::foreach_impl(
        std::forward<F>(func), std::forward<Tuple>(tuple),
        std::make_index_sequence<
            std::tuple_size<typename std::decay_t<Tuple>::parent>::value>{});
  }

  /* ------------------------------------------------------------------------ */
  template <class F, class Tuple,
            std::enable_if_t<is_tuple<std::decay_t<Tuple>>::value> * = nullptr>
  constexpr auto transform(F && func, Tuple && tuple) -> decltype(auto) {
    return details::transform_impl(
        std::forward<F>(func), std::forward<Tuple>(tuple),
        std::make_index_sequence<
            std::tuple_size<std::decay_t<Tuple>>::value>{});
  }

  template <class F, class Tuple,
            std::enable_if_t<is_tuple<std::decay_t<Tuple>>::value> * = nullptr>
  constexpr auto transform(F && func, Tuple && tuple1, Tuple && tuple2)
      -> decltype(auto) {
    return details::transform_impl(
        std::forward<F>(func), std::forward<Tuple>(tuple1),
        std::forward<Tuple>(tuple2),
        std::make_index_sequence<
            std::tuple_size<std::decay_t<Tuple>>::value>{});
  }

  template <
      class F, class Tuple,
      std::enable_if_t<is_named_tuple<std::decay_t<Tuple>>::value> * = nullptr>
  constexpr auto transform(F && func, Tuple && tuple) -> decltype(auto) {
    return details::transform_named_impl(
        std::forward<F>(func), std::forward<Tuple>(tuple),
        std::make_index_sequence<
            std::tuple_size<typename std::decay_t<Tuple>::parent>::value>{});
  }

  /* ------------------------------------------------------------------------ */
  template <class Tuple> auto flatten(Tuple && tuples) -> decltype(auto) {
    return details::flatten(std::forward<Tuple>(tuples),
                            std::make_index_sequence<
                                std::tuple_size<std::decay_t<Tuple>>::value>());
  }

  namespace details {
    template <size_t n, class Tuple>
    auto dynamic_get_impl(size_t i, Tuple && tuple) -> decltype(auto) {
      constexpr auto size = std::tuple_size<std::decay_t<Tuple>>::value;
      if (i == n) {
        return std::get<n>(tuple);
      }

      if (n == size - 1) {
        throw std::out_of_range("Tuple element out of range.");
      }

      return dynamic_get_impl<(n < size - 1 ? n + 1 : 0)>(i, tuple);
    }
  } // namespace details

  /* ------------------------------------------------------------------------ */
  template <typename... Ts> struct cat {
    using type = decltype(std::tuple_cat(std::declval<Ts>()...));
  };

  template <typename... T> using cat_t = typename cat<T...>::type;

  /* ------------------------------------------------------------------------ */
  template <template <typename> class Pred, typename... Ts> struct filter {};

  template <template <typename> class Pred, typename T>
  struct filter<Pred, std::tuple<T>> {
    using type =
        std::conditional_t<Pred<T>::value, std::tuple<T>, std::tuple<>>;
  };

  template <template <typename> class Pred, typename T, typename... Ts>
  struct filter<Pred, std::tuple<T, Ts...>> {
    using type = cat_t<typename filter<Pred, std::tuple<T>>::type,
                       typename filter<Pred, std::tuple<Ts...>>::type>;
  };

  template <template <typename> class Pred, typename... Ts>
  using filter_t = typename filter<Pred, Ts...>::type;

  /* ------------------------------------------------------------------------ */
  template <class Tuple>
  constexpr auto dynamic_get(std::size_t i, Tuple && tuple) -> decltype(auto) {
    return details::dynamic_get_impl<0>(i, std::forward<Tuple>(tuple));
  }

  /* ------------------------------------------------------------------------ */
  template <class Tuple, class... Vals,
            std::enable_if_t<is_tuple<std::decay_t<Tuple>>::value> * = nullptr>
  auto append(Tuple && tuple, Vals &&... vals) -> decltype(auto) {
    return details::append_impl(
        std::make_index_sequence<std::tuple_size<std::decay_t<Tuple>>::value>{},
        FWD(tuple), std::forward<Vals>(vals)...);
  }

  template <
      class Tuple, class... Vals,
      std::enable_if_t<is_named_tuple<std::decay_t<Tuple>>::value> * = nullptr>
  auto append(Tuple && tuple, Vals &&... vals) -> decltype(auto) {
    return details::append_named_impl(
        std::make_index_sequence<
            std::tuple_size<typename std::decay_t<Tuple>::parent>::value>{},
        FWD(tuple), std::forward<Vals>(vals)...);
  }

  /* ------------------------------------------------------------------------ */
  template <std::size_t nth, class Tuple,
            std::enable_if_t<is_tuple<std::decay_t<Tuple>>::value> * = nullptr>
  auto remove(Tuple && tuple) -> decltype(auto) {
    static_assert(nth < std::tuple_size<Tuple>(),
                  "You are trying to remove a non existing entry");
    return details::remove_impl<nth>(
        std::make_index_sequence<nth>{},
        std::make_index_sequence<std::tuple_size<Tuple>::value - nth - 1>{},
        std::forward<Tuple>(tuple));
  }

  template <
      typename Tag, class Tuple,
      std::enable_if_t<is_named_tuple<std::decay_t<Tuple>>::value> * = nullptr>
  auto remove(Tuple && tuple) -> decltype(auto) {
    static_assert(tuple.has(Tag{}),
                  "You are trying to remove a non existing entry");
    constexpr auto nth = tuple.get_element_index(Tag{});
    return details::remove_named_impl<nth>(
        std::make_index_sequence<nth>{},
        std::make_index_sequence<
            std::tuple_size<typename std::decay_t<Tuple>::parent>::value - nth -
            1>{},
        std::forward<Tuple>(tuple));
  }

  template <std::size_t nth, class Tuple, class Value,
            std::enable_if_t<is_tuple<std::decay_t<Tuple>>::value> * = nullptr>
  auto replace(Tuple && tuple, Value && value) -> decltype(auto) {
    static_assert(nth < std::tuple_size<Tuple>(),
                  "You are trying to replace a non existing entry");
    return details::replace_impl<nth>(
        std::make_index_sequence<nth>{},
        std::make_index_sequence<std::tuple_size<Tuple>::value - nth - 1>{},
        std::forward<Tuple>(tuple), std::forward<Value>(value));
  }

  template <typename Tag, class Tuple, class Value,
            std::enable_if_t<is_named_tuple<std::decay_t<Tuple>>::value and
                             not is_named_tag_v<Value>> * = nullptr>
  auto replace(Tuple && tuple, Value && value) -> decltype(auto) {
    static_assert(tuple::has<Tag, Tuple>(),
                  "You are trying to replace a non existing entry");
    constexpr std::size_t nth =
        std::decay_t<Tuple>::template get_element_index<Tag>();
    return details::replace_named_impl<nth, typename Tag::_tag>(
        std::make_index_sequence<nth>{},
        std::make_index_sequence<
            std::tuple_size<typename std::decay_t<Tuple>::parent>::value - nth -
            1>{},
        std::forward<Tuple>(tuple), std::forward<Value>(value));
  }

  template <class Tuple, class Value,
            std::enable_if_t<is_named_tuple_v<std::decay_t<Tuple>> and
                             is_named_tag_v<Value>> * = nullptr>
  auto replace(Tuple && tuple, Value && value) -> decltype(auto) {
    using Tag = decltype(make_named_tag<typename std::decay_t<Value>::_tag>());
    return replace<Tag>(std::forward<Tuple>(tuple),
                        std::forward<decltype(value._value)>(value._value));
  }

} // namespace tuple

} // namespace AKANTU_ITERATORS_NAMESPACE

#undef FWD

#endif /* AKANTU_AKA_TUPLE_TOOLS_HH */
