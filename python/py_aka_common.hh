/**
 * @file   py_aka_common.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Fri May 03 2019
 * @date last modification: Tue Sep 29 2020
 *
 * @brief  pybind11 interface to aka_common.hh
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include <pybind11/pybind11.h>

#ifndef AKANTU_PY_AKA_COMMON_HH_
#define AKANTU_PY_AKA_COMMON_HH_

namespace akantu {

void register_enums(pybind11::module & mod);
void register_initialize(pybind11::module & mod);
void register_functions(pybind11::module & mod);

} // namespace akantu

#endif
