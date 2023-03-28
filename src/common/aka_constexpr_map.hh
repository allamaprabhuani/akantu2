/**
 * Copyright (©) 2022-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
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
 */

#include <array>
#include <tuple>
#include <type_traits>

#ifndef AKANTU_CONSTEXPR_MAP_HH
#define AKANTU_CONSTEXPR_MAP_HH

namespace akantu {
namespace details {
  template <class InputIt, class UnaryPredicate>
  constexpr InputIt my_find_if(InputIt first, InputIt last,
                               UnaryPredicate && p) {
    for (; first != last; ++first) {
      if (std::forward<UnaryPredicate>(p)(*first)) {
        return first;
      }
    }
    return last;
  }

  // Author Jason Turner C++ Weekly ep 233
  template <typename Key, typename Value, std::size_t Size>
  struct ConstexprMap {
    std::array<std::pair<Key, Value>, Size> data;

    [[nodiscard]] constexpr Value at(const Key & key) const {
      const auto it =
          my_find_if(data.begin(), data.end(),
                     [&key](const auto & val) { return val.first == key; });

      if (it != data.end()) {
        return it->second;
      } else {
        throw std::range_error("Key out of range");
      }
    }

    [[nodiscard]] constexpr auto find(const Key & key) const {
      const auto it =
          my_find_if(data.begin(), data.end(),
                     [&key](const auto & val) { return val.first == key; });

      return it;
    }

    [[nodiscard]] constexpr auto begin() const { return data.begin(); }
    [[nodiscard]] constexpr auto end() const { return data.end(); }
  };

  // magic_switch from
  // https://stackoverflow.com/questions/39915986/solutions-for-dynamic-dispatch-on-unrelated-types
  template <class Function, class DynamicType, class Tuple,
            class DefaultFunction, std::size_t... Is>
  [[gnu::visibility("hidden")]] constexpr decltype(auto) static_switch_dispatch(
      const Tuple &, Function && function, const DynamicType & type,
      DefaultFunction && default_function, std::index_sequence<Is...> /*is*/) {
    auto * function_pointer = std::addressof(function);
    using FunctionPointer = decltype(function_pointer);
    using Ret = decltype(function(std::tuple_element_t<0, Tuple>{}));
    using TableEntry = Ret (*)(FunctionPointer);

    constexpr std::array<std::pair<DynamicType, std::size_t>, sizeof...(Is)>
        data{{{std::tuple_element_t<Is, Tuple>::value, Is}...}};

    constexpr auto map =
        ConstexprMap<DynamicType, std::size_t, data.size()>{{data}};

    constexpr TableEntry table[] = {
        [](FunctionPointer function_pointer) -> Ret {
          return (*function_pointer)(std::tuple_element_t<Is, Tuple>{});
        }...};

    auto it = map.find(type);
    if (it != map.end()) {
      return table[it->second](function_pointer);
    } else {
      return default_function(type);
    }
  }
} // namespace details
} // namespace akantu

#endif // AKANTU_CONSTEXPR_MAP_HH
