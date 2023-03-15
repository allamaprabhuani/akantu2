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

#include <iostream>

#include "aka_constexpr_map.hh"

int main() {
  akantu::details::static_switch_dispatch(
      std::tuple<std::integral_constant<int, 1>,
                 std::integral_constant<int, 2>>{},

      [](auto && type) {
        std::cout << std::decay_t<decltype(type)>::value << std::endl;
      },
      2, [](auto && /*type*/) { std::cout << "Default" << std::endl; },
      std::make_index_sequence<2>{});
}
