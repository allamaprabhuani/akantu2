#include <iostream>

#include "aka_constexpr_map.hh"

int main() {
  akantu::details::static_switch_dispatch(
      std::tuple<std::integral_constant<int, 1>,
                 std::integral_constant<int, 2>>{},

      [](auto && type) {
        std::cout << std::decay_t<decltype(type)>::value << std::endl;
      },
      2, [](auto && type) { std::cout << "Default" << std::endl; },
      std::make_index_sequence<2>{});
}
