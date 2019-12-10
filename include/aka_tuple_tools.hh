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
 * @section LICENSE
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
#include "aka_str_hash.hh"
/* -------------------------------------------------------------------------- */
#include <tuple>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_AKA_TUPLE_TOOLS_HH
#define AKANTU_AKA_TUPLE_TOOLS_HH

#ifndef AKANTU_ITERATORS_NAMESPACE
#define AKANTU_ITERATORS_NAMESPACE akantu
#endif

namespace AKANTU_ITERATORS_NAMESPACE {

namespace tuple {
  /* ------------------------------------------------------------------------ */
  namespace details {
    template <std::size_t N> struct Foreach {
      template <class Tuple>
      static inline bool not_equal(Tuple && a, Tuple && b) {
        if (std::get<N - 1>(std::forward<Tuple>(a)) ==
            std::get<N - 1>(std::forward<Tuple>(b)))
          return false;
        return Foreach<N - 1>::not_equal(std::forward<Tuple>(a),
                                         std::forward<Tuple>(b));
      }
    };

    /* ---------------------------------------------------------------------- */
    template <> struct Foreach<0> {
      template <class Tuple>
      static inline bool not_equal(Tuple && a, Tuple && b) {
        return std::get<0>(std::forward<Tuple>(a)) !=
               std::get<0>(std::forward<Tuple>(b));
      }
    };

    template <typename... Ts>
    decltype(auto) make_tuple_no_decay(Ts &&... args) {
      return std::tuple<Ts...>(std::forward<Ts>(args)...);
    }

    template <class F, class Tuple, std::size_t... Is>
    void foreach_impl(F && func, Tuple && tuple,
                      std::index_sequence<Is...> &&) {
      (void)std::initializer_list<int>{
          (std::forward<F>(func)(std::get<Is>(std::forward<Tuple>(tuple))),
           0)...};
    }

    template <class F, class Tuple, std::size_t... Is>
    decltype(auto) transform_impl(F && func, Tuple && tuple,
                                  std::index_sequence<Is...> &&) {
      return make_tuple_no_decay(
          std::forward<F>(func)(std::get<Is>(std::forward<Tuple>(tuple)))...);
    }
  } // namespace details

  /* ------------------------------------------------------------------------ */
  template <class Tuple> bool are_not_equal(Tuple && a, Tuple && b) {
    return details::Foreach<std::tuple_size<std::decay_t<Tuple>>::value>::
        not_equal(std::forward<Tuple>(a), std::forward<Tuple>(b));
  }

  template <class F, class Tuple> void foreach (F && func, Tuple && tuple) {
    return details::foreach_impl(
        std::forward<F>(func), std::forward<Tuple>(tuple),
        std::make_index_sequence<
            std::tuple_size<std::decay_t<Tuple>>::value>{});
  }

  template <class F, class Tuple>
  decltype(auto) transform(F && func, Tuple && tuple) {
    return details::transform_impl(
        std::forward<F>(func), std::forward<Tuple>(tuple),
        std::make_index_sequence<
            std::tuple_size<std::decay_t<Tuple>>::value>{});
  }

  namespace details {
    template <class Tuple, std::size_t... Is>
    decltype(auto) flatten(Tuple && tuples, std::index_sequence<Is...>) {
      return std::tuple_cat(std::get<Is>(tuples)...);
    }
  } // namespace details

  template <class Tuple> decltype(auto) flatten(Tuple && tuples) {
    return details::flatten(std::forward<Tuple>(tuples),
                            std::make_index_sequence<
                                std::tuple_size<std::decay_t<Tuple>>::value>());
  }

  template <template <class> class traits, typename... Ts>
  struct conjunction {};

  template <template <class> class traits, typename... Ts>
  struct conjunction<traits, std::tuple<Ts...>>
      : public aka::conjunction<typename traits<Ts>::value...> {};

} // namespace tuple

namespace tuple {
  namespace details {
    /* ---------------------------------------------------------------------- */
    template <typename tag, typename type> struct named_tag {
      using _tag = tag;   ///< key
      using _type = type; ///< value type

      template <typename T>
      explicit named_tag(T && value) : _value(std::forward<T>(value)) {}

      type _value;
    };

/* ---------------------------------------------------------------------- */
#if (defined(__GNUC__) || defined(__GNUG__))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif
    template <typename tag> struct named_tag_proxy {
      using _tag = tag;

      template <typename T> decltype(auto) operator=(T && value) {
        return named_tag<tag, T>{std::forward<T>(value)};
      }
    };
#if (defined(__GNUC__) || defined(__GNUG__))
#pragma GCC diagnostic pop
#endif

    template <typename T> struct is_named_tag : public std::false_type {};
    template <typename tag>
    struct is_named_tag<named_tag_proxy<tag>> : public std::true_type {};

    /* ---------------------------------------------------------------------- */
    /// Named tuple is just tuple of named params
    template <typename... Params>
    struct named_tuple : public std::tuple<Params...> {
      template <typename... Args>
      named_tuple(Args &&... args)
          : std::tuple<Args...>(std::forward<Args>(args)...) {}

    private:
      template <typename tag, std::size_t Idx,
                std::enable_if_t<Idx == sizeof...(Params)> * = nullptr>
      static constexpr std::size_t get_element_index() noexcept {
        return -1;
      }

      template <typename tag, std::size_t Idx,
                std::enable_if_t<(Idx < sizeof...(Params))> * = nullptr>
      static constexpr std::size_t get_element_index() noexcept {
        using element_t = std::tuple_element_t<Idx, std::tuple<Params...>>;
        return (std::is_same<typename element_t::_tag, tag>::value)
                   ? Idx
                   : get_element_index<tag, Idx + 1>();
      }

    public:
      template <typename NT,
                std::enable_if_t<is_named_tag<NT>::value> * = nullptr>
      constexpr decltype(auto) get(NT &&) noexcept {
        const auto index = get_element_index<typename NT::_tag, 0>();
        static_assert((index != -1), "wrong named_tag");
        return (std::get<index>(*this)._value);
      }

      template <typename NT,
                std::enable_if_t<is_named_tag<NT>::value> * = nullptr>
      constexpr decltype(auto) get(NT &&) const noexcept {
        const auto index = get_element_index<typename NT::_tag, 0>();
        static_assert((index != -1), "wrong named_tag");

        return std::get<index>(*this)._value;
      }
    };

  } // namespace details

  template <typename... Params>
  constexpr decltype(auto) make_named_tuple(Params &&... params) noexcept {
    return details::named_tuple<Params...>(std::forward<Params>(params)...);
  }

  template <typename tag> constexpr decltype(auto) make_named_tag() noexcept {
    return details::named_tag_proxy<tag>{};
  }

  template <size_t HashCode> constexpr decltype(auto) get() {
    return make_named_tag<std::integral_constant<size_t, HashCode>>();
  }

  template <size_t HashCode, class Tuple>
  constexpr decltype(auto) get(Tuple && tuple) {
    return tuple.get(get<HashCode>());
  }

  template <typename Param, typename Tuple,
            std::enable_if_t<details::is_named_tag<Param>::value> * = nullptr>
  constexpr decltype(auto) get(Tuple && tuple) noexcept {
    return tuple.template get<typename Param::hash>();
  }


#if defined(__INTEL_COMPILER)
// intel warnings here
#elif defined(__clang__)
// clang warnings here
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wgnu-string-literal-operator-template"
#elif (defined(__GNUC__) || defined(__GNUG__))
// gcc warnings here
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#endif

  /// this is a GNU exstension
  /// http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3599.html
  template <class CharT, CharT... chars>
  constexpr decltype(auto) operator"" _n() {
    return make_named_tag<string_literal<CharT, chars...>>();
  }

#if defined(__INTEL_COMPILER)
#elif defined(__clang__)
#pragma clang diagnostic pop
#elif (defined(__GNUC__) || defined(__GNUG__))
#pragma GCC diagnostic pop
#endif

} // namespace tuple
} // namespace AKANTU_ITERATORS_NAMESPACE

#endif /* AKANTU_AKA_TUPLE_TOOLS_HH */
