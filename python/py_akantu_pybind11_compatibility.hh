/**
 * Copyright (©) 2021-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef PY_AKANTU_PYBIND11_COMPATIBILITY_HH_
#define PY_AKANTU_PYBIND11_COMPATIBILITY_HH_

#if not defined(PYBIND11_OVERRIDE)
#define PYBIND11_OVERRIDE PYBIND11_OVERLOAD
#define PYBIND11_OVERRIDE_NAME PYBIND11_OVERLOAD_NAME
#define PYBIND11_OVERRIDE_PURE PYBIND11_OVERLOAD_PURE
#define PYBIND11_OVERRIDE_PURE_NAME PYBIND11_OVERLOAD_PURE_NAME
#endif

#endif // PY_AKANTU_PYBIND11_COMPATIBILITY_HH_
