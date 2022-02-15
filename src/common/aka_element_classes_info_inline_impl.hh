/**
 * @file   aka_element_classes_info_inline_impl.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Jun 18 2015
 * @date last modification: Tue Sep 29 2020
 *
 * @brief  Implementation of the streaming fonction for the element classes
 * enums
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_config.hh"
#include "aka_tuple_tools.hh"
/* -------------------------------------------------------------------------- */
#include <tuple>
#include <unordered_map>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_AKA_ELEMENT_CLASSES_INFO_INLINE_IMPL_HH_
#define AKANTU_AKA_ELEMENT_CLASSES_INFO_INLINE_IMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
// BOOST PART: TOUCH ONLY IF YOU KNOW WHAT YOU ARE DOING
#define AKANTU_BOOST_CASE_MACRO(r, macro, _type)                               \
  case _type: {                                                                \
    macro(_type);                                                              \
    break;                                                                     \
  }

#define AKANTU_BOOST_LIST_SWITCH(macro1, list1, var)                           \
  do {                                                                         \
    switch (var) {                                                             \
      BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_CASE_MACRO, macro1, list1)            \
    default: {                                                                 \
      AKANTU_ERROR("Type (" << var << ") not handled by this function");       \
    }                                                                          \
    }                                                                          \
  } while (0)

#define AKANTU_BOOST_LIST_SWITCH_NO_DEFAULT(macro1, list1, var)                \
  do {                                                                         \
    switch (var) {                                                             \
      BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_CASE_MACRO, macro1, list1)            \
    case _not_defined:                                                         \
      break;                                                                   \
    case _max_element_type:                                                    \
      break;                                                                   \
    }                                                                          \
  } while (0)

#define AKANTU_BOOST_LIST_SWITCH_CONSTEXPR(macro1, list1, var)                 \
  do {                                                                         \
    switch (var) {                                                             \
      BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_CASE_MACRO, macro1, list1)            \
    default: {                                                                 \
      macro1(_not_defined);                                                    \
    }                                                                          \
    }                                                                          \
  } while (0)

#define AKANTU_BOOST_ELEMENT_SWITCH(macro1, list1)                             \
  AKANTU_BOOST_LIST_SWITCH(macro1, list1, type)

#define AKANTU_BOOST_ELEMENT_SWITCH_NO_DEFAULT(macro1, list1)                  \
  AKANTU_BOOST_LIST_SWITCH_NO_DEFAULT(macro1, list1, type)

#define AKANTU_BOOST_ELEMENT_SWITCH_CONSTEXPR(macro1, list1)                   \
  AKANTU_BOOST_LIST_SWITCH_CONSTEXPR(macro1, list1, type)

#define AKANTU_BOOST_ALL_ELEMENT_SWITCH(macro)                                 \
  AKANTU_BOOST_ELEMENT_SWITCH(macro, AKANTU_ALL_ELEMENT_TYPE)

#define AKANTU_BOOST_ALL_ELEMENT_SWITCH_NO_DEFAULT(macro)                      \
  AKANTU_BOOST_ELEMENT_SWITCH_NO_DEFAULT(macro, AKANTU_ALL_ELEMENT_TYPE)

#define AKANTU_BOOST_ALL_ELEMENT_SWITCH_CONSTEXPR(macro)                       \
  AKANTU_BOOST_ELEMENT_SWITCH_CONSTEXPR(macro, AKANTU_ALL_ELEMENT_TYPE)

#define AKANTU_BOOST_LIST_MACRO(r, macro, type) macro(type)

#define AKANTU_BOOST_APPLY_ON_LIST(macro, list)                                \
  BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_LIST_MACRO, macro, list)

#define AKANTU_BOOST_ALL_ELEMENT_LIST(macro)                                   \
  AKANTU_BOOST_APPLY_ON_LIST(macro, AKANTU_ALL_ELEMENT_TYPE)

#define AKANTU_GET_ELEMENT_LIST(kind) AKANTU##kind##_ELEMENT_TYPE

#define AKANTU_BOOST_KIND_ELEMENT_SWITCH(macro, kind)                          \
  AKANTU_BOOST_ELEMENT_SWITCH(macro, AKANTU_GET_ELEMENT_LIST(kind))

// BOOST_PP_SEQ_TO_LIST does not exists in Boost < 1.49
#define AKANTU_GENERATE_KIND_LIST(seq)                                         \
  BOOST_PP_TUPLE_TO_LIST(BOOST_PP_SEQ_SIZE(seq), BOOST_PP_SEQ_TO_TUPLE(seq))

#define AKANTU_ELEMENT_KIND_BOOST_LIST                                         \
  AKANTU_GENERATE_KIND_LIST(AKANTU_ELEMENT_KIND)

#define AKANTU_BOOST_ALL_KIND_LIST(macro, list)                                \
  BOOST_PP_LIST_FOR_EACH(AKANTU_BOOST_LIST_MACRO, macro, list)

#define AKANTU_BOOST_ALL_KIND(macro)                                           \
  AKANTU_BOOST_ALL_KIND_LIST(macro, AKANTU_ELEMENT_KIND_BOOST_LIST)

#define AKANTU_BOOST_ALL_KIND_SWITCH(macro)                                    \
  AKANTU_BOOST_LIST_SWITCH(macro, AKANTU_ELEMENT_KIND, kind)

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

AKANTU_ENUM_OUTPUT_STREAM(
    ElementType, AKANTU_ALL_ELEMENT_TYPE(_not_defined)(_max_element_type))
AKANTU_ENUM_INPUT_STREAM(ElementType, AKANTU_ALL_ELEMENT_TYPE)

AKANTU_ENUM_OUTPUT_STREAM(InterpolationType, AKANTU_INTERPOLATION_TYPES)
AKANTU_ENUM_INPUT_STREAM(InterpolationType, AKANTU_INTERPOLATION_TYPES)

AKANTU_ENUM_OUTPUT_STREAM(ElementKind, AKANTU_ELEMENT_KIND)
AKANTU_ENUM_INPUT_STREAM(ElementKind, AKANTU_ELEMENT_KIND)

template <::akantu::ElementType t>
using element_type_t = std::integral_constant<::akantu::ElementType, t>;

#define OP_CAT(s, data, elem) BOOST_PP_CAT(_element_type, elem)

// creating a type instead of a using helps to debug
#define AKANTU_DECLARE_ELEMENT_TYPE_STRUCT(r, data, elem)                      \
  struct BOOST_PP_CAT(_element_type, elem)                                     \
      : public element_type_t<::akantu::elem> {};

BOOST_PP_SEQ_FOR_EACH(AKANTU_DECLARE_ELEMENT_TYPE_STRUCT, _,
                      AKANTU_ALL_ELEMENT_TYPE)

#undef AKANTU_DECLARE_ELEMENT_TYPE_STRUCT

template <ElementKind kind> struct ElementTypes;

template <> struct ElementTypes<_ek_regular> {
  using type = std::tuple<BOOST_PP_SEQ_ENUM(
      BOOST_PP_SEQ_TRANSFORM(OP_CAT, _, AKANTU_ek_regular_ELEMENT_TYPE))>;
};

#if defined(AKANTU_COHESIVE_ELEMENT)
template <> struct ElementTypes<_ek_cohesive> {
  using type = std::tuple<BOOST_PP_SEQ_ENUM(
      BOOST_PP_SEQ_TRANSFORM(OP_CAT, _, AKANTU_ek_cohesive_ELEMENT_TYPE))>;
};
#endif

#if defined(AKANTU_STRUCTURAL_MECHANICS)
template <> struct ElementTypes<_ek_structural> {
  using type = std::tuple<BOOST_PP_SEQ_ENUM(
      BOOST_PP_SEQ_TRANSFORM(OP_CAT, _, AKANTU_ek_structural_ELEMENT_TYPE))>;
};
#endif

#undef OP_CAT

template <ElementKind kind>
using ElementTypes_t = typename ElementTypes<kind>::type;

#define OP_CAT(s, data, elem) ElementTypes_t<elem>
using AllElementTypes = tuple::cat_t<BOOST_PP_SEQ_ENUM(
    BOOST_PP_SEQ_TRANSFORM(OP_CAT, _, AKANTU_ELEMENT_KIND))>;
#undef OP_CAT

namespace details {
  // Author Jason Turner C++ Weekly ep 233
  template <typename Key, typename Value, std::size_t Size>
  struct ConstexprMap {
    std::array<std::pair<Key, Value>, Size> data;
    [[nodiscard]] constexpr Value at(const Key & key) const {
      const auto it =
          std::find_if(data.begin(), data.end(),
                       [&key](const auto & val) { return val.first == key; });

      if (it != data.end()) {
        return it->second;
      } else {
        throw std::range_error("Key not in map");
      }
    }
  };

  // magic_switch from
  // https://stackoverflow.com/questions/39915986/solutions-for-dynamic-dispatch-on-unrelated-types
  template <class Function, class DynamicType, class Tuple, std::size_t... Is>
  decltype(auto) static_switch_dispatch(const Tuple &, Function && function,
                                        const DynamicType & type,
                                        std::index_sequence<Is...> /*is*/) {
    auto * function_pointer = std::addressof(function);
    using FunctionPointer = decltype(function_pointer);
    using Ret = decltype(function(std::tuple_element_t<0, Tuple>{}));
    using TableEntry = Ret (*)(FunctionPointer);

    // needed since gcc 5.4.1 does not understand the expansion other why
    auto unpack = [&](auto I) {
      using type = std::tuple_element_t<decltype(I)::value, Tuple>;
      return std::make_pair(type::value,
                            [](FunctionPointer function_pointer) -> Ret {
                              return (*function_pointer)(type{});
                            });
    };

    static constexpr std::array<std::pair<DynamicType, TableEntry>,
                                sizeof...(Is)>
        data{unpack(std::integral_constant<std::size_t, Is>{})...};

    static constexpr auto map =
        ConstexprMap<DynamicType, TableEntry, data.size()>{{data}};

    return map.at(type)(function_pointer);
  }
} // namespace details

template <class Tuple, class Function, class DynamicType>
constexpr decltype(auto) tuple_dispatch(Function && function,
                                        const DynamicType & type) {
  return details::static_switch_dispatch(
      Tuple{}, std::forward<Function>(function), type,
      std::make_index_sequence<std::tuple_size<Tuple>::value>{});
}

} // namespace akantu

#endif /* AKANTU_AKA_ELEMENT_CLASSES_INFO_INLINE_IMPL_HH_ */
