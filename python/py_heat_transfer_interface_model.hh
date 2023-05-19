/**
 * @file   py_heat_transfer_interface_model.hh
 *
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 *
 * @date creation: Tue May 16 2023
 * @date last modification: Tue May 16 2023
 *
 * @brief  pybind11 interface to HeatTransferInterfaceModel
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

/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>

#ifndef AKANTU_PY_HEAT_TRANSFER_INTERFACE_MODEL_HH_
#define AKANTU_PY_HEAT_TRANSFER_INTERFACE_MODEL_HH_

namespace akantu {

void register_heat_transfer_interface_model(pybind11::module & mod);

} // namespace akantu

#endif // AKANTU_PY_HEAT_TRANSFER_INTERFACE_MODEL_HH_
