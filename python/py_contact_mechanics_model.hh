/**
 * @file   py_contact_mechanics_model.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Thu Jun 20 2019
 * @date last modification: Sat Dec 12 2020
 *
 * @brief  Contact mechanics python binding
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

#ifndef __AKANTU_PY_CONTACT_MECHANICS_MODEL_HH__
#define __AKANTU_PY_CONTACT_MECHANICS_MODEL_HH__

namespace akantu {
void register_contact_mechanics_model(pybind11::module & mod);
} // namespace akantu

#endif //  __AKANTU_PY_CONTACT_MECHANICS_MODEL_HH__
