#include <algorithm>
#include <exception>
#include <iostream>
#include <tuple>
#include <type_traits>
#include <utility>

namespace details {
// Author Jason Turner C++ Weekly ep 233
template <typename Key, typename Value, std::size_t Size> struct ConstexprMap {
  std::array<std::pair<Key, Value>, Size> data;
  [[nodiscard]] constexpr Value at(const Key & key) const {
    const auto it =
        std::find_if(data.begin(), data.end(),
                     [&key](const auto & val) { return val.first == key; });

    if (it != data.end()) {
      return it->second;
    } else {
      throw std::range_error("Key out of range");
    }
  }

  [[nodiscard]] constexpr auto find(const Key & key) const {
    const auto it =
        std::find_if(data.begin(), data.end(),
                     [&key](const auto & val) { return val.first == key; });

    return it;
  }

  [[nodiscard]] constexpr auto begin() const { return data.begin(); }
  [[nodiscard]] constexpr auto end() const { return data.end(); }
};

template <class Function, class DynamicType, class Tuple, std::size_t... Is>
constexpr decltype(auto)
static_switch_dispatch(Function && function, const DynamicType & type,
                       const Tuple &, std::index_sequence<Is...> /*is*/) {
  auto * function_pointer = std::addressof(function);
  using FunctionPointer = decltype(function_pointer);
  using Ret = decltype(function(std::tuple_element_t<0, Tuple>{}));
  using TableEntry = Ret (*)(FunctionPointer);

  constexpr std::array<std::pair<DynamicType, TableEntry>, sizeof...(Is)> data{
      {{std::tuple_element_t<Is, Tuple>::value,
        [](FunctionPointer function_pointer) -> Ret {
          return (*function_pointer)(std::tuple_element_t<Is, Tuple>{});
        }}...}};

  constexpr auto map =
      ConstexprMap<DynamicType, TableEntry, data.size()>{{data}};

  auto it = map.find(type);
  if (it != map.end()) {
    return it->second(function_pointer);
  }
}

} // namespace details

int main() {
  details::static_switch_dispatch(
      [](auto && type) {
        std::cout << std::decay_t<decltype(type)>::value << std::endl;
      },
      2,
      std::tuple<std::integral_constant<int, 1>,
                 std::integral_constant<int, 2>>{},
      std::make_index_sequence<2>{});
}
